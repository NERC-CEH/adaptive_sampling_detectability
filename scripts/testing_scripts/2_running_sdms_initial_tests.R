
## running initial sdms

slurm_run_sim_sdm(index = rep(n_species, length(n_communities)*length(models)*length(data_type)),
           spdata = rep(sprintf(
             paste0(dirs$commpath, community_version, simulation_run_name,"/", community_version,
                    "community_%i_%i_sim/", ifelse(data_type!='initial', paste0(AS_version, '_'), ''), community_version, "community_%i_%i_sim_%s.rds"), 
             rep(rep(n_communities, each = length(models)), each = length(data_type)), max(n_species), 
             rep(rep(n_communities, each = length(models)), each = length(data_type)), max(n_species), data_type
           ), each = length(n_species)), # location of all the community data
           model = rep(rep(rep(models, length(n_communities)), each = length(data_type)), each = length(n_species)), # which models to use
           data_type = rep(rep(data_type, length(n_communities)*length(models)), each = length(n_species)), 
           writeRas = FALSE,
           GB = TRUE,
           community_version = community_version,
           simulation_run_name = simulation_run_name,
           AS_version = AS_version,
           n_communities = rep(rep(rep(n_communities, each = length(models)), each = length(data_type)), each = length(n_species)),
           n_species = max(n_species)
)

