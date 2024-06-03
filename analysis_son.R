# Install and load required packages
library(lubridate)
library(readxl)
library(patchwork)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)  # for diversity calculation
library(glmmTMB)
library(broom.mixed)
library(caret)


###################
# GET THE DATASET #
###################

# Run the data_wrangling.R

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
  group_by(week, site) %>% 
  summarise(across(all_of(ACOUSTIC_VAR), mean, na.rm=TRUE),
            n_calls = sum(tags_number),
            n_species = n_distinct(tags, na.rm = TRUE),
            hill_number = exp(diversity(table(tags), index = "shannon")))  # Hill number for q = 1

#################################
##########  PLOTS ##############
################################

# Overall plots
# Plot for number of calls per week
p1 <- df_processed %>%
  ggplot(aes(x = week, y = n_calls)) + 
  geom_boxplot(fill = "blue", alpha = 0.5) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  ylab("N. Calls") 
p1

# Plot for species richness per week
p2 <- df_processed %>%
  ggplot(aes(x = week, y = n_species)) + 
  geom_boxplot(fill = "blue", alpha = 0.5) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  ylab("Sp. Rich") 
p2

# Plot for functional number of species per week
p4 <- df_processed %>%
  ggplot(aes(x = week, y = hill_number)) + 
  geom_boxplot(fill = "blue", alpha = 0.5) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  ylab("Hill Nb") 
p4

# Plot for indices
p3 <- df_processed %>% 
  pivot_longer(!c(n_calls, week, site), names_to = "indice", values_to = "value") %>% 
  ggplot(aes(y = value, x=week)) + 
  geom_boxplot(aes(color = indice, fill = indice), alpha = 0.7) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  xlab("Week") +
  ylab("Value of Indice") +
  facet_grid(rows = vars(indice), scales = "free_y")
p3
# Combine the overall plots
(p1 / p2 / p4 / p3) + plot_layout(heights = c(1, 1, 1, 6))

#################################################################################
# Make the df for number of calls, species richness and number of func. species #
#################################################################################

# Convert "None" to NA and prepare the data frame
df_model <- df_annotated %>% 
  mutate(tags = na_if(tags, "None")) %>%  # Convert "None" to NA
  filter(confidence > CONFIDENCE) %>% 
  select(doy, week, site, all_of(ACOUSTIC_VAR), tags) %>% 
  mutate(tags_number = ifelse(is.na(tags), 0, 1)) %>% 
  group_by(doy, week, site) %>% 
  summarise(across(all_of(ACOUSTIC_VAR), mean, na.rm=TRUE),
            n_calls = sum(tags_number),
            n_species = n_distinct(tags, na.rm = TRUE),
            hill_number = exp(diversity(table(tags), index = "shannon")))

# Standardize acoustic variables
preProcessParams <- preProcess(df_model[, ACOUSTIC_VAR], method = c("center", "scale"))
df_model[, ACOUSTIC_VAR] <- predict(preProcessParams, df_model[, ACOUSTIC_VAR])

# Add quadratic terms
df_model <- df_model %>%
  mutate(across(all_of(ACOUSTIC_VAR), ~ .^2, .names = "{.col}_sq"))

# Site with the highest number of calls
df_model %>% select(c(site, n_calls)) %>% arrange(desc(n_calls))
# Bergen sites seems to be the ones with the highest activity

# Filter sites with sufficient data
min_replications <- 2  # Minimum number of replications required per site
site_counts <- df_model %>% group_by(site) %>% summarise(count = n())
sufficient_sites <- site_counts %>% filter(count >= min_replications) %>% pull(site)
df_model <- df_model %>% filter(site %in% sufficient_sites)

# Combine data for modeling
df_model <- df_model %>%
  select(week, site, hill_number, n_calls, n_species, all_of(ACOUSTIC_VAR), ends_with("_sq"))

##########################################
# Build the Mixed-Effects Models         #
##########################################

# Model for hill_number
mixed_model_hill <- glmmTMB(hill_number ~ EVNsp + EVNsp_sq + MFC + MFC_sq + ROItotal + ROItotal_sq +
                              ACI + ACI_sq + BI + BI_sq + roughness + roughness_sq + (1 | site),
                            data = df_model)

# Model for n_calls
mixed_model_calls <- glmmTMB(n_calls ~ EVNsp + EVNsp_sq + MFC + MFC_sq + ROItotal + ROItotal_sq +
                               ACI + ACI_sq + BI + BI_sq + roughness + roughness_sq + (1 | site),
                             data = df_model)

# Model for n_species
mixed_model_species <- glmmTMB(n_species ~ EVNsp + EVNsp_sq + MFC + MFC_sq + ROItotal + ROItotal_sq +
                                 ACI + ACI_sq + BI + BI_sq + roughness + roughness_sq + (1 | site),
                               data = df_model)

# Summarize the mixed-effects models
summary(mixed_model_hill)
summary(mixed_model_calls)
summary(mixed_model_species)

# Extract model coefficients
coefficients_hill <- tidy(mixed_model_hill)
coefficients_calls <- tidy(mixed_model_calls)
coefficients_species <- tidy(mixed_model_species)
print(coefficients_hill)
print(coefficients_calls)
print(coefficients_species)

###################################
# Model validation and prediction #
###################################

# Predict and visualize for hill_number
df_model$predicted_hill_number <- predict(mixed_model_hill, newdata = df_model, allow.new.levels = TRUE)

# Plot observed vs predicted for hill_number
ggplot(df_model, aes(x = hill_number, y = predicted_hill_number)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Observed Hill Number", y = "Predicted Hill Number") +
  theme_minimal()

# Plot residuals for hill_number
ggplot(df_model, aes(x = week, y = residuals(mixed_model_hill))) +
  geom_point() +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Week", y = "Residuals") +
  theme_minimal()

# Predict and visualize for n_calls
df_model$predicted_n_calls <- predict(mixed_model_calls, newdata = df_model, allow.new.levels = TRUE)

# Plot observed vs predicted for n_calls
ggplot(df_model, aes(x = n_calls, y = predicted_n_calls)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Observed Number of Calls", y = "Predicted Number of Calls") +
  theme_minimal()

# Plot residuals for n_calls
ggplot(df_model, aes(x = week, y = residuals(mixed_model_calls))) +
  geom_point() +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Week", y = "Residuals") +
  theme_minimal()

# Predict and visualize for n_species
df_model$predicted_n_species <- predict(mixed_model_species, newdata = df_model, allow.new.levels = TRUE)

# Plot observed vs predicted for n_species
ggplot(df_model, aes(x = n_species, y = predicted_n_species)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Observed Number of Species", y = "Predicted Number of Species") +
  theme_minimal()

# Plot residuals for n_species
ggplot(df_model, aes(x = week, y = residuals(mixed_model_species))) +
  geom_point() +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Week", y = "Residuals") +
  theme_minimal()

##########################################
# Save Plots Per Site                    #
##########################################

# Create a directory to save the plots
dir.create("plots_per_site", showWarnings = FALSE)

# Loop through each site and create a plot
sites <- unique(df_model$site)
for (site in sites) {
  df_site_processed <- df_model %>% filter(site == !!site) 

  p1_site <- df_site_processed %>%
    ggplot(aes(x = doy)) + 
    geom_point(aes(y=n_calls), fill = "blue", alpha = 0.5) +
    geom_point(aes(y=predicted_n_calls), color="red")
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank()) +
    ylab("N. Calls") +
    ggtitle(paste("Site:", site))

  p2_site <- df_site_processed %>%
    ggplot(aes(x = doy)) + 
    geom_point(aes(y = n_species), fill = "blue", alpha = 0.5) +
    geom_point(aes(y=predicted_n_species), color="red")
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank()) +
    ylab("Sp. Rich") +
    ggtitle(paste("Site:", site))
  
  p4_site <- df_site_processed %>%
    ggplot(aes(x = doy)) + 
    geom_point(aes(y = hill_number), fill = "blue", alpha = 0.5) +
    geom_point(aes(y=predicted_hill_number), color="red")
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank()) +
    ylab("Hill Nb") +
    ggtitle(paste("Site:", site))

  p3_site <- df_site_processed %>%
    ungroup() %>% 
    select(!c(week, hill_number, n_species, predicted_hill_number, predicted_n_calls, predicted_n_species)) %>% 
    pivot_longer(!c(n_calls, doy, site), names_to = "indice", values_to = "value") %>% 
    ggplot(aes(y = value, x = doy)) + 
    geom_point(aes(color = indice, fill = indice), alpha = 0.7) +
    geom_line(aes(color = indice, fill = indice, group = indice), alpha = 0.7) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    xlab("Week") +
    ylab("Value of Indice") +
    facet_grid(rows = vars(indice), scales = "free_y") +
    ggtitle(paste("Site:", site    )) +
    ggtitle(paste("Site:", site))
  
  # Combine the plots for the site
  combined_site_plot <- (p1_site / p2_site / p4_site / p3_site) + plot_layout(heights = c(1, 1, 1, 6))
  
  # Save the plot
  ggsave(filename = paste0("plots_per_site/hill_number_site_", site, ".png"), plot = combined_site_plot, width = 20, height = 15)
}

