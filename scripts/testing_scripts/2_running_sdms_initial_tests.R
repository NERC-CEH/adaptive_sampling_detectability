
## running initial sdms
rm(list = ls())

dirs <- config::get("Cdrive_thomas")

source("scripts/lotus/lotus_functions/slurm_run_sim_sdm_function.R")

# ## testing for a single species/community
# 
# index = 1
# community_data = "outputs/communities/v1narrow_nichebreadth_community/v1community_1_100_sim/v1community_1_100_sim_initial.rds" # location of the community data
# model = "rf" # which models to use
# data_type = "initial"
# writeRas = FALSE
# GB = TRUE
# environmental_subset = (2/3) # what proportion of environmental layers should be used for modelling?
# community_version = "v1"
# simulation_run_name = "equal_prevalence_community"
# AS_version = "asv1"
# n_communities = 1
# n_species = 100
# function_path = dirs$functionpath
# outpath = dirs$commpath
# envpath = dirs$envpath

community_version_name = "v1" # Which community version
simulation_run_name = 'narrow_nichebreadth_community'
n_species = 1:100
n_communities = 1
as_version = "asv1"

## run the function for one community
sdm <- slurm_run_sim_sdm(
  index = 1,
  community_data = paste0(dirs$outpath, community_version_name, simulation_run_name, "/", 
                          community_version_name, 
                          sprintf(paste0("community_%i_%i_sim/", community_version_name, 
                                         "community_%i_%i_sim_initial.rds"), n_communities, max(n_species), 
                                  n_communities, max(n_species))), # location of the community data
  model = "rf", # which models to use
  data_type = "initial", 
  writeRas = FALSE,
  GB = TRUE,
  environmental_subset = (2/3), # what proportion of environmental layers should be used for modelling?
  community_version = community_version_name,
  simulation_run_name = simulation_run_name,
  AS_version = as_version,
  n_communities = n_communities,
  n_species = max(n_species),
  function_path = dirs$functionpath,
  outpath = dirs$commpath,
  envpath = dirs$envpath
)


## code to run for multiple communities -- use slurm script




