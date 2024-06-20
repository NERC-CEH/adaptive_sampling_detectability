#' # Run all simulated species

#' @param index Index of the species the models are being run for (used to pick out the species from the species list that's created from the community data)
#' @param community_data Location on disc of where the community data are stored
#' @param model Models to run, one or multiple of lr, gam, rf
#' @param env_data Path to the environmental data
#' @param data_type The type of data being used to distinguish between initial models and the adaptive sampling methods. One or multiple of c("initial_AS_none", "initial_AS_uncertainty", "initial_AS_prevalence", "initial_AS_unc_plus_prev", "initial_AS_unc_plus_recs", "initial_AS_coverage") # 'initial'
#' @param writeRas Logical. Write the model predictions and standard deviations as a spatRast
#' @param extent_crop extent to crop raster - must be able to be converted to a SpatialPolygons object
#' @param extent_crs coordinates for the extent
#' @param environmental_subset what proportion of environmental layers should be used for modelling. If NULL, use all
#' @param community_version_name Name of the community version we are running - so no overwriting occurs
#' @param AS_version  One name for the adaptive sampling round to allow us to create different sampling methods of the same initial community
#' @param simulation_run_name The name of the simulation run - same as slurm_simulate species
#' @param n_communities Vector of the number of communities to go through
#' @param n_species vector of the number of species in each community
#' @param function_path location of the functions to source for modelling - can remove if turn into a package
#' @param outpath Where to store outputs
#' @param envpath where the environmental data are stored
#' @export

slurm_run_sim_sdm <- function(index, 
                              community_data, 
                              model,
                              env_data,
                              data_type, 
                              writeRas, 
                              extent_crop,
                              extent_crs,
                              environmental_subset = (2/3), # what proportion of environmental layers should be used for modelling? If NULL, use all
                              community_version_name, 
                              AS_version, 
                              simulation_run_name, 
                              n_communities, 
                              n_species,
                              function_path,
                              outpath){
  
  library(terra)
  library(virtualspecies)
  library(dismo)
  library(tidyverse)
  library(Rfast)
  library(mgcv)
  library(randomForest)
  library(ENMTools)
  
  source(paste0(function_path,"reformat_simulated_data.R"))
  source(paste0(function_path, "Edited_Rob_Functions.R"))
  source(paste0(function_path,"getpredictions_dfsd.R"))
  
  
  # remove this
  # dirs <- config::get("LOTUSpaths")
  
  ## 1. Simulate distributions (or read in simulated spp)
  
  # load community data
  print("! reading community file")
  community <- readRDS(as.character(community_data))
  
  #' ## 2. Create input data for models
  #' 
  #' We need to extract the virtual species data we simulated and combine into a community dataset. 
  #' For each species we can create a new dataset with pseudoabsences generated from the 
  #' other species in the community
  #' 
  #create a pseudo-absence dataset
  
  #read in raster data for env data
  #read in env data frame
  
  full_env_data <- terra::rast(as.character(env_data))
  full_env_data_df <- as.data.frame(full_env_data, xy = TRUE)
  
  #crop to extent if specified
  if(!is.null(extent_crop)){
    
    e <- as(extent_crop, "SpatialPolygons")
    sp::proj4string(e) <- extent_crs
    
    e.geo <- sp::spTransform(e, terra::crs(env)) # CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs"))

        full_env_data <- terra::crop(env, e.geo)
    full_env_data_df <- as.data.frame(full_env_data, xy = TRUE)
    
  }
  
  presences_df <- reformat_simulated_data(sim_species_out = community, 
                                          year = 2015, 
                                          species_name = 'Sp')
  #head(presences_df)
  
  # get pseudoabsences ofr each species
  print("! Generating pseudoabsences")
  
  species_list <- unique(presences_df$species)
  
  pres_abs <- vector('list', length = length(species_list))
  
  
  for(s in 1:length(species_list)){
    
    pres_abs[[s]] <- cpa(spdat = presences_df, species = species_list[s], 
                         matchPres = FALSE, nAbs = 10000,
                         minYear = 2000, maxYear = 2017, recThresh = 1)
    
  }
  
  names(pres_abs) <- species_list
  
  
  #' Pseudoabsence weighting dependent on model, Thomas' code has this accounted for
  
  #' 3. Run DECIDE models 
  #' 
  #' Initially we can run just the logistic regression models as a test
  
  print("! starting modelling")
  
  #loop over all 10 species - set up for LOTUS
  
  #number of bootstraps
  k = 10
  
  #species index
  sp_list <- names(pres_abs)
  species <- sp_list[index]
  
  #subset envdata for species of interest  
  env_data <- subset(full_env_data, subset = community[[index]]$variables)
  
  # use subset of data for modelling to reduce model performance?
  if(!is.null(environmental_subset)) {
    
    print(paste("! using only a random subset of",  round(dim(env_data)[3]*environmental_subset), 
                "of the",  dim(env_data)[3], "environmental layers"))
    
    #use only a proportion of the environmental data to run the models to reduce model fit
    env_index <- sample(1:dim(env_data)[3], 
                        size = round(dim(env_data)[3]*environmental_subset), 
                        replace = FALSE)
    
    # subset environmental data for the model run
    env_data <- env_data[[env_index]]
  }
  
  # run species distribution model
  sdm <- fsdm(species = species, model = model,
              climDat = env_data, spData = pres_abs, knots_gam = 4,
              k = k, 
              write =  FALSE, outPath = outpath)
  
  #predictions
  print("! predicting to full dataset")
  
  # when model is rf - it's normal to have two columns - accounted for later!!!
  preds1 <- get_predictions_dfsd(model_outs = sdm, 
                                 model = model, # model that was run to create the models in model_outs
                                 env_data = full_env_data_df)
  
  ## save files ##
  species_name <- gsub(pattern = ' ', replacement = '_', species) # get species name without space
  
  # create the output path to be in the community that the species belongs to
  outPath <- paste0(outpath, community_version_name, simulation_run_name,'/', community_version_name, "community_",n_communities,"_", n_species, "_sim/")
  
  #' Calculate very simple DECIDE score - prediction * standard deviation
  
  DECIDE_score <- preds1$mean_predictions*preds1$sd_predictions
  
  
  if(writeRas == TRUE){
    
    # create a new location to save all rasters
    dir.create(paste0(outPath, community_version_name, 'sdm_plots/'))
    
    # save prediction raster
    terra::writeRaster(x = terra::rast(cbind(full_env_data_df$x,full_env_data_df$y,preds1$mean_predictions), type = "xyz"), 
                       filename = paste0(outPath, community_version_name, 'sdm_plots/', model, "_SDMs_", species_name, "_meanpred.grd"),
                       format = 'raster', overwrite = T)
    
    # save sd raster
    terra::writeRaster(x = terra::rast(cbind(full_env_data_df$x, full_env_data_df$y, preds1$sd_predictions), type = "xyz"), 
                       filename = paste0(outPath, community_version_name, 'sdm_plots/', model, "_SDMs_", species_name, "_sdpred.grd"),
                       format = 'raster', overwrite = T)
    
  }
  
  # write AUC to file for easy-access
  #write.csv(x = data.frame(raw_AUC = sdm$AUC,
  #                         meanAUC = sdm$meanAUC),
  #          file = paste0(outPath, model, "_SDMs_", species_name, "_AUC_values.csv"))
  
  # write data to file too
  #write.csv(x = sdm$Data,
  #         file = paste0(outPath, model, "_SDMs_", species_name, "_Data.csv"))
  
  # save subset model output
  # remove data from model output
  sdm$Data <- NULL
  
  community_name <- basename(community_data) #strsplit(as.character(community_data),"\\/")[[1]][10]
  
  # output of model to store
  model_output <- list(community_version_name,
                       AS_version,
                       community = community_name, 
                       species = species_name,
                       model = model,
                       sdm_output = lapply(sdm$Bootstrapped_models, function(x) summary(x)),
                       number_validations = k,
                       meanAUC = sdm$meanAUC,
                       predictions = data.frame(x = full_env_data_df$x, y = full_env_data_df$y, mean = preds1$mean_predictions, sd = preds1$sd_predictions, DECIDE_score = DECIDE_score))
  
  # create a new directory to store species SDMs
  dir.create(paste0(outPath, community_version_name, "species_models/"))
  
  print("#####     Saving output     #####") ## to check if the process is hanging on lotus
  
  save(model_output, file = paste0(outPath, community_version_name, "species_models/", ifelse(data_type!='initial', paste0(AS_version, '_'), ''), community_version_name, model, "_SDMs_GBnew_", species_name, "_", data_type, ".rdata"))
  
  print("#####     Output saved     #####")
  
  
  print(paste("Read files using:", 
              paste0("'readRDS(", 
                     paste0(outPath, community_version_name, "species_models/", ifelse(data_type!='initial', paste0(AS_version, '_'), ''), community_version_name, model, "_SDMs_GBnew_", species_name, "_", data_type, ".rdata'"))))
  
  
}
