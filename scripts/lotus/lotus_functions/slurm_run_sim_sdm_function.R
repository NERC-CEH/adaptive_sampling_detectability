#' # Run all simulated species

#' @param index
#' @param community_data
#' @param model
#' @param data_type
#' @param writeRas
#' @param GB
#' @param environmental_subset = (2/3) # what proportion of environmental layers should be used for modelling? If NULL, use all
#' @param community_version
#' @param AS_version
#' @param simulation_run_name
#' @param n_communities
#' @param n_species
#' @param function_path
#' @param outpath
#' @param envpath
#' @export


slurm_run_sim_sdm <- function(index, 
                              community_data, 
                              model, 
                              data_type, 
                              writeRas, 
                              GB, 
                              environmental_subset = (2/3), # what proportion of environmental layers should be used for modelling? If NULL, use all
                              community_version, 
                              AS_version, 
                              simulation_run_name, 
                              n_communities, 
                              n_species,
                              function_path,
                              outpath,
                              envpath){
  
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
  
  if(GB) {
    hbv_y <- terra::rast(paste0(envpath,"envdata_1km_no_corr_noNA.grd"))
    hbv_df <- as.data.frame(hbv_y, xy = TRUE) # read.csv(paste0(envdata, "hbv_df_1km.csv"))
  } else { ### sort this out to crop extent for testing - not toooooo sure what this does?
    hbv_y <- terra::rast(paste0(dirs$inpath,"hbv_y.grd")) 
    hbv_df <- readRDS(paste0(dirs$inpath, "hbv_df.rds"))
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
  env_data <- subset(hbv_y, subset = community[[index]]$variables)
  
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
                                 env_data = hbv_df)
  
  ## save files ##
  species_name <- gsub(pattern = ' ', replacement = '_', species) # get species name without space
  
  # create the output path to be in the community that the species belongs to
  outPath <- paste0(outpath, community_version, simulation_run_name,'/', community_version, "community_",n_communities,"_", n_species, "_sim/")
  
  #' Calculate very simple DECIDE score - prediction * standard deviation
  
  DECIDE_score <- preds1$mean_predictions*preds1$sd_predictions
  
  
  if(writeRas == TRUE){
    
    # create a new location to save all rasters
    dir.create(paste0(outPath, community_version, 'sdm_plots/'))
    
    # save prediction raster
    writeRaster(x = rast(cbind(hbv_df$x,hbv_df$y,preds1$mean_predictions), type = "xyz"), 
                filename = paste0(outPath, community_version, 'sdm_plots/', model, "_SDMs_", species_name, "_meanpred.grd"),
                format = 'raster', overwrite = T)
    
    # save sd raster
    writeRaster(x = rasterFromXYZ(cbind(hbv_df$x, hbv_df$y, preds1$sd_predictions)), 
                filename = paste0(outPath, community_version, 'sdm_plots/', model, "_SDMs_", species_name, "_sdpred.grd"),
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
  model_output <- list(community_version,
                       AS_version,
                       community = community_name, 
                       species = species_name,
                       model = model,
                       sdm_output = lapply(sdm$Bootstrapped_models, function(x) summary(x)),
                       number_validations = k,
                       meanAUC = sdm$meanAUC,
                       predictions = data.frame(x = hbv_df$x, y = hbv_df$y, mean = preds1$mean_predictions, sd = preds1$sd_predictions, DECIDE_score = DECIDE_score))
  
  # create a new directory to store species SDMs
  dir.create(paste0(outPath, community_version, "species_models/"))
  
  print("#####     Saving output     #####") ## to check if the process is hanging on lotus
  
  save(model_output, file = paste0(outPath, community_version, "species_models/", ifelse(data_type!='initial', paste0(AS_version, '_'), ''), community_version, model, "_SDMs_GBnew_", species_name, "_", data_type, ".rdata"))
  
  print("#####     Output saved     #####")
  
  
  print(paste("Read files using:", 
              paste0("'readRDS(", 
                     paste0(outPath, community_version, "species_models/", ifelse(data_type!='initial', paste0(AS_version, '_'), ''), community_version, model, "_SDMs_GBnew_", species_name, "_", data_type, ".rdata'"))))
  
  
}
