rm(list = ls())

## testing simulate species
source("scripts/lotus/lotus_functions/slurm_simulate_species_function.R")


# ### testing
# env_data = "data/environmental_data/edited_data/envdata_1km_no_corr_noNA.grd"
# outPath = "outputs/communties/" 
# seed = 1 # community number
# max_samp = 20000 # max number of observations per species
# n_env = 10 # number of environmental variables from which to sample for species generation
# n = 50 # number of species per community
# det_prob = 0.2 # detection probability
# sample_across_species = TRUE # whether to sample the same locations across all species (i.e. list structure) or sample different locations for each species
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
simulated_species <- simulate_species(
  env_data = "data/environmental_data/edited_data/envdata_1km_no_corr_noNA.grd",
  outPath = "outputs/communties/", 
  seed = 1, # community number
  max_samp = 20000, # max number of observations per species
  n_env = 10, # number of environmental variables from which to sample for species generation
  n = 10, # number of species per community
  det_prob = 0.2, # detection probability
  sample_across_species = TRUE, # whether to sample the same locations across all species (i.e. list structure) or sample different locations for each species
  effort = "data/effort_layers/butterfly_1km_effort_layer.grd", # sampling effort layer
  background = "MeanDiRange", # a layer which contains a value for each cell in the region of interest
  community_version_name = "v1", # Which community version
  simulation_run_name = 'first_community',
  write = TRUE
)


###