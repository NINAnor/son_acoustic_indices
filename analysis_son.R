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
library(usdm)

###################
# GET THE DATASET #
###################

# Run the data_wrangling.R

####################################################
# VARIABLE COMBINING ALL THE VARIABLES OF INTEREST #
####################################################

ACOUSTIC_VAR <- c("EVNspFract", "EVNspMean", "EVNspCount", 
                  "MFC", "ROItotal", "ROIcover", "ACI", "BI", 
                  "roughnessMean", "roughnessMedian", "roughnesStd",
                  "EVNtFraction", "EVNtMean",
                  "ACTtCount", "ACTtMean")
CONFIDENCE = 0.8

############################################################
# How do the acoustic indicators correlate with each other #
############################################################

df_annotated %>% 
  select(all_of(ACOUSTIC_VAR)) %>% 
  drop_na() %>% 
  cor(.)
df_annotated[sapply(df_annotated, is.infinite)] <- NA

#### Variance Inflation Factors ####

# From the manual for the usdm package "Collinearity causes instability in parameter estimation in regression-type 
#models. The VIF is based on the square of the multiple correlation coefficient resulting from regressing a 
#predictor variable against all other predictor variables. If a variable has a strong linear relationship 
#with at least one of the other variables, the correlation coefficient will be close to 1, and VIF for that 
#variable would be large."

# Two options for VIFs: vifcor and vifstep. The latter is the option used by Geue & Thomassen (2020) 
#Functions "exclude highly collinear variables in a stepwise procedure".

# "vifcor" first finds the pair of variables with the maximum linear correlation (greater than "th" threshold), 
#and excludes the one with the greater VIF. The procedure is repeated until no variable has a high 
#correlation coefficient (greater than threshold) with other variables.

# "vifstep" calculates the VIF for all variables, excludes the one with highest VIF (greater than threshold), 
#and repeats the procedure until no variables with VIF greater than "th" remains.

# In Species Distribution Modelling the norm is to retain predictors with VIFs of < 10. 
#In other modelling areas, then ViFs of 3 or 5 are more normal.

options(scipen = 999)
set.seed(123) # Not certain if there is random process? Results seem to vary?

PREDICTORS <- df_annotated %>% 
  select(all_of(ACOUSTIC_VAR), doy, latitude) %>% 
  drop_na()

# Predictors can only be numeric
VIF.STEP <- vifstep(PREDICTORS, th = 5)
VIF.COR <- vifcor(PREDICTORS, th = 0.7)

# Returns a VIF object, examine different outputs.
VIF.STEP.MATRIX <- data.frame(VIF.STEP@corMatrix)
VIF.COR.MATRIX <- data.frame(VIF.COR@corMatrix)

# Check which variables are still included.
colnames(PREDICTORS)
row.names(VIF.STEP.MATRIX)
row.names(VIF.COR.MATRIX)

PREDICTORS <- PREDICTORS[VIF.STEP@results$Variables]
ACOUSTIC_VAR_VIF <- PREDICTORS %>% 
  select(all_of(ACOUSTIC_VAR[ACOUSTIC_VAR %in% colnames(PREDICTORS)]))

ACOUSTIC_VAR_VIF <- colnames(ACOUSTIC_VAR_VIF)

# Selection 
ACOUSTIC_VAR_VIF <- ACOUSTIC_VAR


#################################################################################
# Make the df for number of calls, species richness and number of func. species #
#################################################################################

# Convert "None" to NA and prepare the data frame
df_model <- df_annotated %>% 
  mutate(tags = na_if(tags, "None")) %>%  # Convert "None" to NA
  filter(confidence > CONFIDENCE) %>% 
  select(doy, latitude, week, site, all_of(ACOUSTIC_VAR_VIF), tags) %>% 
  mutate(tags_number = ifelse(is.na(tags), 0, 1)) %>% 
  group_by(doy, week, site) %>% 
  summarise(across(all_of(ACOUSTIC_VAR_VIF), mean, na.rm=TRUE),
            n_calls = sum(tags_number),
            n_species = n_distinct(tags, na.rm = TRUE),
            latitude = mean(latitude), 
            hill_number = exp(diversity(table(tags), index = "shannon"))) %>% 
  drop_na()

# Standardize acoustic variables
preProcessParams <- preProcess(df_model[, ACOUSTIC_VAR_VIF], method = c("center", "scale"))
df_model[, ACOUSTIC_VAR_VIF] <- predict(preProcessParams, df_model[, ACOUSTIC_VAR_VIF])

# Add quadratic terms
df_model <- df_model %>%
  mutate(across(all_of(ACOUSTIC_VAR_VIF), ~ .^2, .names = "{.col}_sq"))

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
  select(doy, week, site, hill_number, n_calls, n_species, all_of(ACOUSTIC_VAR_VIF), ends_with("_sq"), latitude)

##########################################
# Build the Mixed-Effects Models         #
##########################################
HILL_FORMULA <- paste("hill_number~", paste(c(ACOUSTIC_VAR_VIF, df_model %>% ungroup() %>% select(ends_with("_sq")) %>% colnames(), "doy", "latitude", "(1|site)"), collapse="+"), sep="")
N_CALL_FORMULA <- paste("n_calls~", paste(c(ACOUSTIC_VAR_VIF, df_model %>% ungroup() %>% select(ends_with("_sq")) %>% colnames(), "doy", "latitude", "(1|site)"), collapse="+"), sep="")
N_SPECIES_FORMULA <- paste("n_species~", paste(c(ACOUSTIC_VAR_VIF, df_model %>% ungroup() %>% select(ends_with("_sq")) %>% colnames(), "doy", "latitude", "(1|site)"), collapse="+"), sep="")

mixed_model_hill <- glmmTMB(as.formula(HILL_FORMULA),
                            data = df_model,
                            control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000)),
                            family = Gamma(link = "log"))

# Model for n_calls
mixed_model_calls <- glmmTMB(as.formula(N_CALL_FORMULA),
                             data = df_model,
                             control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000)),
                             family=poisson)
# Model for n_species
mixed_model_species <- glmmTMB(as.formula(N_SPECIES_FORMULA),
                               data = df_model,
                               control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000)),
                               family = poisson)

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
hill_predictions <- predict(mixed_model_hill, newdata = df_model, allow.new.levels = TRUE, se.fit = TRUE)
df_model$predicted_hill_number <- hill_predictions$fit
df_model$hill_number_se <- hill_predictions$se.fit

# Plot observed vs predicted for hill_number
ggplot(df_model, aes(x = hill_number, y = predicted_hill_number)) +
  geom_point() +
  geom_errorbar(aes(ymin = predicted_hill_number - 1.96 * hill_number_se,
                    ymax = predicted_hill_number + 1.96 * hill_number_se), width = 0.2) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Observed Hill Number", y = "Predicted Hill Number") +
  theme_minimal()

# Plot residuals for hill_number
ggplot(df_model, aes(x = doy, y = residuals(mixed_model_hill))) +
  geom_point() +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Week", y = "Residuals") +
  theme_minimal()

# Predict and visualize for n_calls
calls_predictions <- predict(mixed_model_calls, newdata = df_model, allow.new.levels = TRUE, se.fit = TRUE)
df_model$predicted_n_calls <- calls_predictions$fit
df_model$calls_se <- calls_predictions$se.fit

# Plot observed vs predicted for n_calls
ggplot(df_model, aes(x = n_calls, y = predicted_n_calls)) +
  geom_point() +
  geom_errorbar(aes(ymin = predicted_n_calls - 1.96 * calls_se,
                    ymax = predicted_n_calls + 1.96 * calls_se), width = 0.2) +
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
species_predictions <- predict(mixed_model_species, newdata = df_model, allow.new.levels = TRUE, se.fit = TRUE)
df_model$predicted_n_species <- species_predictions$fit
df_model$species_se <- species_predictions$se.fit

# Plot observed vs predicted for n_species
ggplot(df_model, aes(x = n_species, y = predicted_n_species)) +
  geom_point() +
  geom_errorbar(aes(ymin = predicted_n_species - 1.96 * species_se,
                    ymax = predicted_n_species + 1.96 * species_se), width = 0.2) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Observed Number of Species", y = "Predicted Number of Species") +
  theme_minimal()

# Plot residuals for n_species
ggplot(df_model, aes(x = week, y = residuals(mixed_model_species))) +
  geom_point() +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Week", y = "Residuals") +
  theme_minimal()

