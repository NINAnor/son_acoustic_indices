library(lubridate)
library(readxl)
library(ranger)
library(vip)
library(patchwork)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)  # for diversity calculation
library(broom)

###################
# GET THE DATASET #
###################

# Run the to_duckdb.R

####################################################
# VARIABLE COMBINING ALL THE VARIABLES OF INTEREST #
####################################################

ACOUSTIC_VAR <- c("EVNsp", "MFC", "ROItotal", "ROIcover", "ACI", "BI", "roughness")
CONFIDENCE = 0.8

############################################################
# How do the acoustic indicators correlate with each other #
############################################################

df_annotated %>% 
  select(all_of(ACOUSTIC_VAR)) %>% 
  drop_na() %>% 
  cor(.)

###################################
####### DATA VIZ ##################
###################################

# indices at locations
df_annotated %>% 
  select(all_of(ACOUSTIC_VAR), site) %>% 
  group_by(site) %>% 
  summarise_all(.funs="mean") %>% 
  pivot_longer(!site, names_to = "indice", values_to = "value") %>% 
  ggplot(., aes(x=site, y=value)) +
  geom_point() +
  facet_wrap(~indice, scales="free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Indices vs time
df_annotated %>% 
  select(all_of(ACOUSTIC_VAR), week) %>% 
  group_by(week) %>% 
  summarise_all(.funs="mean") %>% 
  pivot_longer(!week, names_to = "indice", values_to = "value") %>% 
  ggplot(., aes(x=week, y=value)) +
  geom_point() +
  facet_wrap(~indice, scales="free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#################################################################################
# Make the df for number of calls, species richness and number of func. species #
#################################################################################

# DF that summarises the number of calls
df_processed <- df_annotated %>% 
  filter(tags != "None") %>%
  filter(confidence > CONFIDENCE) %>% 
  select(!filename, site) %>% 
  mutate(tags_number = ifelse(is.na(tags), 0, 1)) %>% 
  group_by(week) %>% 
  summarise(across(all_of(ACOUSTIC_VAR), mean, na.rm=TRUE),
            n = sum(tags_number)) %>% 
  pivot_longer(!c(n, week), names_to = "indice", values_to = "value")

# DF that summarises species richness
df_richness <- df_annotated %>% 
  filter(tags != "None") %>%
  filter(confidence > CONFIDENCE) %>% 
  select(!filename) %>% 
  group_by(week) %>% 
  summarise(across(all_of(ACOUSTIC_VAR), mean, na.rm=TRUE),
            n_species = n_distinct(tags, na.rm = TRUE)) %>% 
  pivot_longer(!c(n_species, week), names_to = "indice", values_to = "value")

# DF that summarises functional number of species using Hill number
df_functional <- df_annotated %>% 
  filter(tags != "None") %>%
  filter(confidence > CONFIDENCE) %>% 
  select(week, all_of(ACOUSTIC_VAR), tags) %>% 
  group_by(week) %>% 
  summarise(across(all_of(ACOUSTIC_VAR), mean, na.rm=TRUE),
            hill_number = exp(diversity(table(tags), index = "shannon")))  # Hill number for q = 1

#################################
##########  PLOTS ##############
################################

# Plot for number of calls per week
p1 <- df_processed %>%
  ggplot(aes(x = week, y = n)) + 
  geom_point(stat = "identity", fill = "blue", alpha = 0.5) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  ylab("N. Calls") 
p1

# Plot for species richness per week
p2 <- df_richness %>%
  ggplot(aes(x = week, y = n_species)) + 
  geom_point(stat = "identity", fill = "blue", alpha = 0.5) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  ylab("Sp. Rich") 
p2

# Plot for functional number of species per week
p4 <- df_functional %>%
  ggplot(aes(x = week, y = hill_number)) + 
  geom_point(stat = "identity", fill = "blue", alpha = 0.5) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  ylab("Hill Nb") 
p4

# Plot for indices
p3 <- df_processed %>%
  ggplot(aes(y = value, x=week)) + 
  geom_point(aes(color = indice, fill = indice), alpha = 0.7) +
  geom_line(aes(color = indice, fill = indice, group=indice), alpha = 0.7) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  xlab("Week") +
  ylab("Value of Indice") +
  facet_grid(rows = vars(indice), scales = "free_y")

# Combine the plots
(p1 / p2 / p4 / p3) + plot_layout(heights = c(1, 1, 1, 6))

# Combine data for modeling
df_model <- df_functional %>%
  select(week, hill_number, all_of(ACOUSTIC_VAR))

##########################################
# Build the Statistical Model            #
##########################################

# Fit the linear regression model
model <- lm(hill_number ~ ., data = df_model %>% select(-week))

# Summarize the model
summary(model)

# Visualize the model diagnostics
par(mfrow = c(2, 2))
plot(model)

# Extract model coefficients
coefficients <- tidy(model)
print(coefficients)

# Predict and visualize
df_model <- df_model %>%
  mutate(predicted_hill_number = predict(model, newdata = df_model))

# Plot observed vs predicted
ggplot(df_model, aes(x = hill_number, y = predicted_hill_number)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Observed Hill Number", y = "Predicted Hill Number") +
  theme_minimal()

# Plot residuals
ggplot(df_model, aes(x = week, y = model$residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Week", y = "Residuals") +
  theme_minimal()
