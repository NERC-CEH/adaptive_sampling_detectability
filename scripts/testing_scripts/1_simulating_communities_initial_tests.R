rm(list = ls())

## testing simulate species
source("scripts/lotus/lotus_functions/slurm_simulate_species_function.R")


# ### testing
# env_data = "data/environmental_data/edited_data/envdata_1km_no_corr_noNA.grd"
# outPath = "outputs/communties/"
# seed = 1 # community number
# max_samp = 20000 # max number of observations per species
# n_env = 10 # number of environmental variables from which to sample for species generation
# n = 5 # number of species per community
# detect_prob = "uniform" # detection probability
# sample_across_species = TRUE # whether to sample the same locations across all species (i.e. list structure) or sample different locations for each species
# niche_breadth = "any"
# effort = "data/effort_layers/butterfly_1km_effort_layer.grd" # sampling effort layer
# background = "MeanDiRange" # a layer which contains a value for each cell in the region of interest
# community_version_name = "v1" # Which community version
# simulation_run_name = 'first_community'
# 
# extent = NULL
# weight_adj = 1
# beta = 0.5
# alpha = -0.05

# run script
simulated_community <- simulate_species(
  env_data = "data/environmental_data/edited_data/envdata_1km_no_corr_noNA.grd",
  outPath = "outputs/communities/", 
  seed = 1, # community number
  max_samp = 20000, # max number of observations per species
  n_env = 10, # number of environmental variables from which to sample for species generation
  n = 50, # number of species per community
  detect_prob = "uniform", # detection probability, "beta" or number between 0-1
  sample_across_species = TRUE, # whether to sample the same locations across all species (i.e. list structure) or sample different locations for each species
  niche_breadth = "any",
  effort = "data/effort_layers/butterfly_1km_effort_layer.grd", # sampling effort layer
  background = "MeanDiRange", # a layer which contains a value for each cell in the region of interest
  community_version_name = "v1", # Which community version
  simulation_run_name = 'equal_prevalence_community',
  write = TRUE
)


# simcom <- readRDS("outputs/communities/v1first_community/v1community_1_50_sim/v1community_1_50_sim_initial.rds")

### check to see distribution of detection probabilities
detdf <- do.call(rbind, lapply(1:length(simulated_community), function(x) {
  data.frame(species = paste0("sp", x),
             prevalence = simulated_community[[x]]$prevalence,
             det_prob = simulated_community[[x]]$detection_probability)
}))

library(tidyverse)

# detdf <- readRDS("outputs/communities/v1equal_prevalence_community/v1community_1_100_sim/v1community_1_100_sim_initial.rds")
# detdf

ggplot(detdf, aes(prevalence, det_prob)) +
  geom_hline(yintercept = median(detdf$det_prob)) +
  geom_vline(xintercept = median(detdf$prevalence)) +
  annotate("text", x = quantile(detdf$prevalence, probs = 0.25),
           y = quantile(detdf$det_prob, probs = 0.75), label = "rare / visible") +
  annotate("text",x = quantile(detdf$prevalence, probs = 0.25),
           y = quantile(detdf$det_prob, probs = 0.25), label = "rare / cryptic") +
  annotate("text", x = quantile(detdf$prevalence, probs = 0.75),
           y = quantile(detdf$det_prob, probs = 0.25), label = "common / cryptic") +
  annotate("text", x = quantile(detdf$prevalence, probs = 0.75),
           y = quantile(detdf$det_prob, probs = 0.75), label = "common / visible") +
  geom_point() +
  theme_bw()

