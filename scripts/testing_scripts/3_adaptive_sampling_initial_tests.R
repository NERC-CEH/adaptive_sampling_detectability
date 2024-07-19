## script to submit adaptive sampling
rm(list = ls())

# source function to run
source("scripts/lotus/lotus_functions/slurm_adaptive_sample_function.R")

library(rslurm)

dirs <- config::get("Cdrive_thomas")


community_version_name = "v1" # Which community version
simulation_run_name = 'narrow_nichebreadth_community'
n_species = 1:100
n_communities = 1

asv_version = data.frame(as_version =  c("asv1", "asv2", "asv3", "asv4"),
                         uptake_value = c(0.1, 0.01, 0, 0.5))

asv = 1

# # ### testing
# community_file = paste0(dirs$commpath, community_version_name, simulation_run_name, "/",
#                         community_version_name,
#                         sprintf(paste0("community_%i_%i_sim/", community_version_name,
#                                        "community_%i_%i_sim_initial.rds"), n_communities, max(n_species),
#                                 n_communities, max(n_species)))
# sdm_path = paste0(dirs$commpath, community_version_name, simulation_run_name, "/",
#                   community_version_name,
#                   sprintf("community_%i_%i_sim/", n_communities, max(n_species)),
#                   community_version_name, "species_models/")
# effort = as.character(paste0(dirs$effortpath,"butterfly_1km_effort_layer.grd"))
# model = c("rf", "gam", "lr")
# background = "AnnualTemp"
# env_data = paste0(dirs$envpath,"envdata_1km_no_corr_noNA.grd")
# extent_crop = NULL
# probability_weight_adj = 1
# weight_adj = 1
# method = "unc_plus_recs"
# uptake = asv_version$uptake_value[asv]
# n = 2000
# community_version = community_version_name
# AS_version = asv_version$as_version[asv]
# outPath = paste0(dirs$commpath, community_version_name, simulation_run_name, "/",
#                  community_version_name, sprintf("community_%i_%i_sim/", n_communities, max(n_species)))



adptive_data <- slurm_adaptive_sample(
  community_file = paste0(dirs$commpath, community_version_name, simulation_run_name, "/", 
                          community_version_name, 
                          sprintf(paste0("community_%i_%i_sim/", community_version_name, 
                                         "community_%i_%i_sim_initial.rds"), n_communities, max(n_species), 
                                  n_communities, max(n_species))),
  sdm_path = paste0(dirs$commpath, community_version_name, simulation_run_name, "/", 
                        community_version_name, 
                        sprintf("community_%i_%i_sim/", n_communities, max(n_species)), 
                        community_version_name, "species_models/"), 
  effort = as.character(paste0(dirs$effortpath,"butterfly_1km_effort_layer.grd")), 
  background = "AnnualTemp", 
  env_data = paste0(dirs$envpath,"envdata_1km_no_corr_noNA.grd"),
  extent_crop = NULL,
  probability_weight_adj = 1,
  weight_adj = 1, 
  method = "unc_plus_detect",
  uptake = asv_version$uptake_value[asv],
  n = 2000,
  community_version = community_version_name,
  AS_version = asv_version$as_version[asv],
  outPath = paste0(dirs$commpath, community_version_name, simulation_run_name, "/", 
                       community_version_name, sprintf("community_%i_%i_sim/", n_communities, max(n_species)))
)

