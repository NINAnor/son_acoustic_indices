library(tidyverse)
library(viridis)

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
  select(all_of(ACOUSTIC_VAR), doy) %>% 
  group_by(doy) %>% 
  summarise_all(.funs="mean") %>% 
  pivot_longer(!doy, names_to = "indice", values_to = "value") %>% 
  ggplot(., aes(x=doy, y=value)) +
  geom_point() +
  facet_wrap(~indice, scales="free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#################################
# PLOT PHENOLOGY OF HILL NUMBER #
#################################

df_model %>% 
  ggplot(aes(x=doy, y=site, fill=hill_number)) + geom_raster() + scale_fill_viridis() + theme_bw()
df_model %>% 
  ggplot(aes(x=doy, y=site, fill=n_calls)) + geom_tile() + scale_fill_viridis()
df_model %>% 
  ggplot(aes(x=doy, y=site, fill=n_species)) + geom_tile() + scale_fill_viridis()

############################################
######### MODEL PLOTS ######################
############################################

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
    geom_point(aes(y = n_calls), fill = "blue", alpha = 0.5) +
    geom_point(aes(y = predicted_n_calls), color = "red", alpha=.5) +
    geom_errorbar(aes(ymin = predicted_n_calls - 1.96 * calls_se,
                      ymax = predicted_n_calls + 1.96 * calls_se), 
                  width = 0.2, color="red", alpha=.5) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank()) +
    ylab("N. Calls") +
    ggtitle(paste("Site:", site))
  
  p2_site <- df_site_processed %>%
    ggplot(aes(x = doy)) + 
    geom_point(aes(y = n_species), fill = "blue", alpha = 0.5) +
    geom_point(aes(y = predicted_n_species), color = "red") +
    geom_errorbar(aes(ymin = predicted_n_species - 1.96 * species_se,
                      ymax = predicted_n_species + 1.96 * species_se), 
                  width = 0.2, color="red", alpha=.5) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank()) +
    ylab("Sp. Rich") +
    ggtitle(paste("Site:", site))
  
  p4_site <- df_site_processed %>%
    ggplot(aes(x = doy)) + 
    geom_point(aes(y = hill_number), fill = "blue", alpha = 0.5) +
    geom_point(aes(y = predicted_hill_number), color = "red") +
    geom_errorbar(aes(ymin = predicted_hill_number - 1.96 * hill_number_se,
                      ymax = predicted_hill_number + 1.96 * hill_number_se), 
                  width = 0.2, color="red", alpha=.5) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank()) +
    ylab("Hill Nb") +
    ggtitle(paste("Site:", site))
  
  p3_site <- df_site_processed %>%
    ungroup() %>% 
    select(!c(week, hill_number, n_species, starts_with("predicted"), ends_with("_sq"))) %>% 
    pivot_longer(!c(n_calls, doy, site), names_to = "indice", values_to = "value") %>% 
    ggplot(aes(y = value, x = doy)) + 
    geom_point(aes(color = indice, fill = indice), alpha = 0.7) +
    geom_line(aes(color = indice, fill = indice, group = indice), alpha = 0.7) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    xlab("Week") +
    ylab("Value of Indice") +
    facet_grid(rows = vars(indice), scales = "free_y") +
    ggtitle(paste("Site:", site))  # Combine the plots for the site
  combined_site_plot <- (p1_site / p2_site / p4_site / p3_site) + plot_layout(heights = c(1, 1, 1, 6))
  
  # Save the plot
  ggsave(filename = paste0("plots_per_site/hill_number_site_", site, ".png"), plot = combined_site_plot, width = 20, height = 15)
}

##########################################
# Global Plots                           #
##########################################

# Create a directory to save the global plots
dir.create("global_plots", showWarnings = FALSE)

df_global <- df_model %>% 
  group_by(doy) %>% 
  summarise(across(c(hill_number, predicted_hill_number, n_calls, predicted_n_calls, n_species, predicted_n_species),
                   list(mean = mean, sd = sd), .names = "{col}_{fn}"))

# Plot observed vs predicted for hill_number with mean and standard deviation
p1_global <- ggplot(df_global, aes(x = doy)) +
  geom_point(aes(y = hill_number_mean), fill = "blue", alpha = 0.5) +
  geom_errorbar(aes(ymin = hill_number_mean - hill_number_sd, ymax = hill_number_mean + hill_number_sd), width = 0.2, color="blue", alpha=.5) +
  geom_point(aes(y = predicted_hill_number_mean), color = "red") +
  geom_errorbar(aes(ymin = predicted_hill_number_mean - predicted_hill_number_sd, ymax = predicted_hill_number_mean + predicted_hill_number_sd), width = 0.2, color="red", alpha=.5) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  ylab("Hill Nb") +
  ggtitle("Global Hill Number")

# Plot observed vs predicted for n_calls with mean and standard deviation
p2_global <- ggplot(df_global, aes(x = doy)) +
  geom_point(aes(y = n_calls_mean), fill = "blue", alpha = 0.5) +
  geom_errorbar(aes(ymin = n_calls_mean - n_calls_sd, ymax = n_calls_mean + n_calls_sd), width = 0.2, color="blue", alpha=.5) +
  geom_point(aes(y = predicted_n_calls_mean), color = "red") +
  geom_errorbar(aes(ymin = predicted_n_calls_mean - predicted_n_calls_sd, ymax = predicted_n_calls_mean + predicted_n_calls_sd), width = 0.2, color="red", alpha=.5) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  ylab("N. Calls") +
  ggtitle("Global Number of Calls")

# Plot observed vs predicted for n_species with mean and standard deviation
p3_global <- ggplot(df_global, aes(x = doy)) +
  geom_point(aes(y = n_species_mean), fill = "blue", alpha = 0.5) +
  geom_errorbar(aes(ymin = n_species_mean - n_species_sd, ymax = n_species_mean + n_species_sd), width = 0.2, color="blue", alpha=.5) +
  geom_point(aes(y = predicted_n_species_mean), color = "red") +
  geom_errorbar(aes(ymin = predicted_n_species_mean - predicted_n_species_sd, ymax = predicted_n_species_mean + predicted_n_species_sd), width = 0.2, color="red", alpha=.5) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  ylab("Sp. Rich") +
  ggtitle("Global Species Richness")

# Plot acoustic indices
p4_global <- df_global %>%
  select(doy, ends_with("_mean")) %>%
  pivot_longer(-doy, names_to = "indice", values_to = "value") %>%
  ggplot(aes(y = value, x = doy)) +
  geom_point(aes(color = indice, fill = indice), alpha = 0.7) +
  geom_line(aes(color = indice, fill = indice, group = indice), alpha = 0.7) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  xlab("Day of Year") +
  ylab("Value of Indice") +
  facet_grid(rows = vars(indice), scales = "free_y") +
  ggtitle("Global Acoustic Indices")

# Combine the global plots
combined_global_plot <- (p1_global / p2_global / p3_global / p4_global) + plot_layout(heights = c(1, 1, 1, 6))

# Save the global plot
ggsave(filename = "global_plots/global_plot.png", plot = combined_global_plot, width = 20, height = 15)


################################################
# PLOTS PREDICTION PER SITE SORTED BY LATITUDE #
################################################

# Can we get phenology out of the indices?
df_model %>% 
  ggplot(aes(x = doy)) + 
  geom_point(aes(y = n_species), fill = "blue", alpha = 0.5) +
  geom_point(aes(y = predicted_n_calls), color = "red", alpha=.5) +
  geom_errorbar(aes(ymin = predicted_n_calls - 1.96 * calls_se,
                    ymax = predicted_n_calls + 1.96 * calls_se), 
                width = 0.2, color="red", alpha=.5) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) + 
  facet_wrap(~site, ncol=3, scales="free")


