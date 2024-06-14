
## running initial sdms
rm(list = ls())

dirs <- config::get("Cdrive_thomas")

source("scripts/lotus/lotus_functions/slurm_run_sim_sdm_function.R")

# ## testing for a single species/community
# index = 1
# community_data = "outputs/communities/v1narrow_nichebreadth_community/v1community_1_100_sim/v1community_1_100_sim_initial.rds" # location of the community data
# model = "rf" # which models to use
# data_type = "initial"
# writeRas = FALSE
# GB = TRUE
# environmental_subset = (2/3) # what proportion of environmental layers should be used for modelling? If NULL, use all
# community_version = "v1"
# simulation_run_name = "equal_prevalance_community"
# AS_version = "asv1"
# n_communities = 1
# n_species = 100
# 
# function_path = dirs$functionpath
# 
# outpath = dirs$outputpath
# 
# envdata = dirs$envpath



## run the function for one community
sdm <- slurm_run_sim_sdm(
  index = 1,
  community_data = "outputs/communities/v1narrow_nichebreadth_community/v1community_1_100_sim/v1community_1_100_sim_initial.rds", # location of the community data
  model = "rf", # which models to use
  data_type = "initial", 
  writeRas = FALSE,
  GB = TRUE,
  environmental_subset = (2/3), # what proportion of environmental layers should be used for modelling?
  community_version = "v1",
  simulation_run_name = "equal_prevalance_community",
  AS_version = "asv1",
  n_communities = 1,
  n_species = 100,
  function_path = dirs$functionpath,
  outpath = dirs$outputpath,
  envpath = dirs$envpath
)


## code to run for multiple communities

# name of the versions we are running - so we're not overwriting things
# one for community-level which includes the community folders and species models folders
community_version = 'v1'

## run multiple adaptive sampling versions at once
asv_version = data.frame(as_version =  c("asv1"),
                         uptake_value = c(0.1))

asv = 1

# One name for the adaptive sampling round to allow us to create different sampling methods of the same initial community
# This doesn't get used if running only the initial model, but does get used when running the adaptive sampling methods.
AS_version = asv_version$as_version[asv]
print(AS_version)

# the name of the simulation run - same as slurm_simulate species
simulation_run_name = 'equal_prevalance_community'


n_species = 1:100 # vector of number of species in each community
n_communities = 1 # number of communities to go through # can submit 11 at a time
models = c('lr', 'gam', 'rf') # models to run
data_type = "initial" # one of "initial" or c("initial_AS_none", "initial_AS_uncertainty", "initial_AS_prevalence", "initial_AS_unc_plus_prev", "initial_AS_unc_plus_recs", "initial_AS_coverage") # 'initial'

# name of the versions we are running - so we're not overwriting things
# one for community-level which includes the community folders and species models folders
community_version = community_version

# One name for the adaptive sampling round to allow us to create different sampling methods of the same initial community
# This doesn't get used if running only the initial model, but does get used when running the adaptive sampling methods.
AS_version = AS_version

# the name of the simulation run - same as slurm_simulate species
simulation_run_name = simulation_run_name

# # set commpath to something for testing (then delete 'dirs$')
# dirs <- data.frame(commpath = 'blob')




#### testing  index = rep(n_species, length(n_communities)*length(models)*length(data_type)),
spdata = rep(sprintf(
  paste0(dirs$commpath, community_version, simulation_run_name,"/", community_version,
         "community_%i_%i_sim/", ifelse(data_type!='initial', paste0(AS_version, '_'), ''), community_version, "community_%i_%i_sim_%s.rds"), 
  rep(rep(n_communities, each = length(models)), each = length(data_type)), max(n_species), 
  rep(rep(n_communities, each = length(models)), each = length(data_type)), max(n_species), data_type
), each = length(n_species))[1] # location of all the community data

model = rep(rep(rep(models, length(n_communities)), each = length(data_type)), each = length(n_species))[1] # which models to use

data_type = rep(rep(data_type, length(n_communities)*length(models)), each = length(n_species))[1] 

writeRas = FALSE

GB = TRUE

community_version = community_version

simulation_run_name = simulation_run_name

AS_version = AS_version

n_communities = rep(rep(rep(n_communities, each = length(models)), each = length(data_type)), each = length(n_species))[1]

n_species = max(n_species)

function_path = dirs$functionpath

outpath = dirs$outputpath

envdata = dirs$envpath


