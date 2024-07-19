## script to submit adaptive sampling
rm(list = ls())

# source function to run
source("scripts/lotus/lotus_functions/slurm_evaluate_function.R")

library(rslurm)

dirs <- config::get("Cdrive_thomas")


community_version = "v1" # Which community version
simulation_run_name = 'narrow_nichebreadth_community'
n_species = 1:100
n_communities = 1

as_versions = c('asv1', 'asv2', 'asv3', 'asv4')

# name of the adaptive sampling version we are looking to evaluate
AS_version = as_versions[1]

# the name of the simulation run - same as slurm_simulate species
simulation_run_name = 'narrow_nichebreadth_community'

n_communities = 1:50

n_species = 1:50

## testing
community_folder = paste0(dirs$commpath, community_version, 
                          simulation_run_name, "/", community_version, 
                          sprintf("community_%i_%i_sim/", n_communities, max(n_species)))[1]
model = "rf, gam, lr"
method =   "initial, none, uncertainty, coverage, detectability, prev_plus_detect, unc_plus_detect"
community_version = community_version
AS_version = AS_version

## index file
evaluation <- slurm_evaluate(community_folder = paste0(dirs$commpath, community_version, 
                                                       simulation_run_name, "/", community_version, 
                                                       sprintf("community_%i_%i_sim/", n_communities, max(n_species)))[1],
                             model = "rf, gam, lr", 
                             method =   "initial, none, uncertainty, coverage, detectability, prev_plus_detect, unc_plus_detect",
                             community_version = community_version,
                             AS_version = AS_version)
