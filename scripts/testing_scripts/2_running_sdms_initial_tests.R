
## running initial sdms
rm(list = ls())

dirs <- config::get("Cdrive_thomas")

source("scripts/lotus/lotus_functions/slurm_run_sim_sdm_function.R")



community_version_name = "v1" # Which community version
simulation_run_name = 'narrow_nichebreadth_community'
n_species = 1:100
n_communities = 1
AS_version = "asv1"
data_type = "initial_AS_unc_plus_recs"
models = "rf"

# ## testing for a single species/community
# 
# index = 1
# community_data = sprintf(
#   paste0(dirs$commpath, community_version_name, simulation_run_name,"/", community_version_name,
#          "community_%i_%i_sim/", ifelse(data_type!='initial', paste0(AS_version, '_'), ''),
#          community_version_name, "community_%i_%i_sim_%s.rds"),
#   rep(rep(n_communities, each = length(models)), each = length(data_type)), max(n_species),
#   rep(rep(n_communities, each = length(models)), each = length(data_type)), max(n_species), data_type
# ) # location of the community data
# model = models# which models to use
# data_type = data_type
# writeRas = FALSE
# GB = TRUE
# environmental_subset = (2/3) # what proportion of environmental layers should be used for modelling?
# community_version = community_version_name
# simulation_run_name = simulation_run_name
# AS_version = AS_version
# n_communities = n_communities
# n_species = max(n_species)
# function_path = dirs$functionpath
# outpath = dirs$commpath
# envpath = dirs$envpat


## run the function for one community
sdm <- slurm_run_sim_sdm(
  index = 1,
  community_data = sprintf(
    paste0(dirs$commpath, community_version_name, simulation_run_name,"/", community_version_name,
           "community_%i_%i_sim/", ifelse(data_type!='initial', paste0(AS_version, '_'), ''), 
           community_version_name, "community_%i_%i_sim_%s.rds"), 
    rep(rep(n_communities, each = length(models)), each = length(data_type)), max(n_species), 
    rep(rep(n_communities, each = length(models)), each = length(data_type)), max(n_species), data_type
  ), # location of the community data
  model = models, # which models to use
  env_data = paste0(dirs$envpath, "envdata_1km_no_corr_noNA.grd"),
  data_type = data_type, 
  writeRas = FALSE,
  GB = TRUE,
  environmental_subset = (2/3), # what proportion of environmental layers should be used for modelling?
  community_version_name = community_version_name,
  simulation_run_name = simulation_run_name,
  AS_version = AS_version,
  n_communities = n_communities,
  n_species = max(n_species),
  function_path = dirs$functionpath,
  outpath = dirs$commpath,
  envpath = dirs$envpath
)


## code to run for multiple communities -- use slurm script




