
## jasmin debugging adaptive sample function

###### first stuff from R script
## script to submit adaptive sampling

# source function to run
source("scripts/slurm_adaptive_sample_function.R")

library(rslurm)

dirs <- config::get("LOTUSpaths_AS")

# name of the version we are running - so we're not overwriting things, keep same as for slurm_run_sim_sdm
version_name = 'v2'

# the name of the simulation run - same as slurm_simulate species
simulation_run_name = 'communities_1km'

# number ofc communities
n_communities = 1:10

# number of species in each community - used only in the parameter file to allow runs with different numbers of species
n_species = 50

# the adaptive sampling methods to use 
method = c("none", "uncertainty", "prevalence", "unc_plus_prev", "unc_plus_recs", "coverage")

# # set outpath and inputs for testing
# outpath = 'broom'
# inputs = 'handle'

pars <- data.frame(community_file = rep(paste0(dirs$outpath, version_name, simulation_run_name, "/", version_name, sprintf(paste0("community_%i_%i_sim/", version_name, "community_%i_%i_sim_initial.rds"), n_communities, max(n_species), n_communities, max(n_species))), each = length(method)), 
                   sdm_path = rep(paste0(dirs$outpath, version_name, simulation_run_name, "/", version_name, sprintf("community_%i_%i_sim/", n_communities, max(n_species)), version_name, "species_models/"), each = length(method)), 
                   effort = as.character(paste0(dirs$inputs,"butterfly_1km_effort_layer.grd")), 
                   background = "AnnualTemp", 
                   env_data = paste0(dirs$inputs,"envdata_1km_no_corr_noNA.grd"), 
                   weight_adj = 1, 
                   method = method, 
                   n = 2000,
                   version_name = version_name,
                   outPath = rep(paste0(dirs$outpath, version_name, simulation_run_name, "/", version_name, sprintf("community_%i_%i_sim/", n_communities, max(n_species))), each = length(method)))

pars$rownum <- 1:nrow(pars)

####### enter in parameters manually
x = 1

community_file = pars$community_file[x]
sdm_path = pars$sdm_path[x]
effort = pars$effort[x]
background = "AnnualTemp"
env_data = pars$env_data[x]
weight_adj = 1
method = pars$method[x]
n = 2000
probability_weight_adj = pars$probability_weight_adj[x]
community_version = pars$community_version[x]
AS_version = AS_version[x]
outPath = pars$outPath[x]
model = c("rf", "gam", "lr")
extent_crop = NULL
rownum = pars$rownum[x]

print(rownum)

#get rdata files with model outputs for each model/species (assuming communities are stored in separate folders) - only read initial models
models <- list.files(path = as.character(sdm_path), pattern = paste0("(",paste(model, sep = "", collapse = "|"),")*initial.rdata"))

#import simulated community data
community <- readRDS(as.character(community_file))








