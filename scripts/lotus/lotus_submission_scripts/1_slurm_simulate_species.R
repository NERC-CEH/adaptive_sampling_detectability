
library(rslurm)

source('../../scripts/functions/lotus_functions/slurm_simulate_species_function.R')

dirs <- config::get("LOTUSpaths")

# a version name that follows all the way through the community
community_version_name = 'v1'

# number of communities
n_communities = 1:10

# detection probability
detect_prob = 0.2

pars <- data.frame(env_data = paste0(dirs$envpath, "envdata_1km_no_corr_noNA.grd"),
                   sample_across_species = TRUE, # whether to sample the same locations across all species (i.e. list structure) or sample different locations for each species
                   
                   n = 50, # number of species per community
                   
                   outPath = dirs$commpath, 
                   seed = n_communities, # community number
                   n_env = 10, # number of environmental variables from which to sample for species generation
                   max_samp = 20000, # max number of observations per species
                   detect_prob = detect_prob, #"uniform", # detection probability, "beta", "uniform" or number between 0-1
                   niche_breadth = "narrow",
                   effort = paste0(dirs$effortpath, "butterfly_1km_effort_layer.grd"), # sampling effort layer
                   background = "MeanDiRange", # a layer which contains a value for each cell in the region of interest
                   community_version_name = community_version_name, # Which community version
                   simulation_run_name = paste0('narrow_breadth_', detect_prob, '_detect_community'),
                   write = TRUE) # the name of the run name - don't change unless changing the resolution of the area of interest.



sjob <- slurm_apply(simulate_species, pars, 
                    jobname = paste0(community_version_name, 'sim_spp'),
                    nodes = nrow(pars), 
                    cpus_per_node = 1, 
                    submit = TRUE,
                    slurm_options = list(partition = "short-serial-4hr",
                                         time = "3:59:59",
                                         mem = "10000",
                                         output = "sim_spp_%a.out",
                                         error = "sim_spp_%a.err",
                                         account = "short4hr"),
                    sh_template = "jasmin_submit_sh.txt")

