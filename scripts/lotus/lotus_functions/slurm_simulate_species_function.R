
#' Simulate species


#' @param env_data 
#' @param sample_across_species
#' @param extent = NULL
#' @param n = 10
#' @param outPath
#' @param seed = NULL
#' @param n_env = NULL
#' @param beta = 0.5
#' @param alpha = -0.05
#' @param max_samp = 1000
#' @param detect_prob = 0.5
#' @param effort = NULL 
#' @param weight_adj = 1 
#' @param background = NULL
#' @param niche_breadth = "narrow"
#' @param community_version_name
#' @param simulation_run_name
#' @param write = TRUE

simulate_species <- function(env_data, 
                             sample_across_species, 
                             extent = NULL, 
                             n = 10, 
                             outPath, 
                             seed = NULL, 
                             n_env = NULL, 
                             beta = 0.5, 
                             alpha = -0.05, 
                             max_samp = 1000, 
                             detect_prob = 0.5,
                             effort = NULL, 
                             weight_adj = 1, 
                             background = NULL,
                             niche_breadth = "narrow",
                             community_version_name, 
                             simulation_run_name,
                             write = TRUE){
  
  library(terra)
  library(virtualspecies)
  
  #read in data
  env <- terra::rast(as.character(env_data))
  
  #set seed if specified
  if(!is.null(seed)){set.seed(seed)}
  # 
  # for(i in 1:10) if(TRUE) print(rbeta(1, 2, 5))
  
  #crop to extent if specified
  if(!is.null(extent)){
    e <- as(extent, "SpatialPolygons")
    sp::proj4string(e) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    
    e.geo <- sp::spTransform(e, CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs"))
    
    env_extent <- terra::crop(env, e.geo)
    
  } else {env_extent <- env}
  
  # subset env data layers
  
  # decide on number of layers to use or range of layers from which to select, 
  # if no limit given use all layers in raster
  if(is.null(n_env)){ 
    env_min = nlayers(env_extent); env_max = nlayers(env_extent)
  } else if (length(n_env) == 1) {
    env_min = n_env; env_max = n_env
  } else if(length(n_env) == 2) {
    env_min = n_env[1]; env_max = n_env[2]
  } else if (length(n_env) > 2) {
    stop("n_env greater than 2: Input either a single number or a range from 
         which to select the number of environmental layers to use in species generation")
  }
  
  # set background if given, can indicate a layer in env_data or be a filepath to a raster
  if(is.numeric(background)){
    bg_layer <- env_extent[[background]]
  } else if ((is.character(background)|is.factor(background)) & !grepl("\\.", background)) {
    bg_layer <- terra::subset(env_extent, as.character(background))
  } else if ((is.character(background)|is.factor(background)) & grepl("\\.", background)) {
    bg_layer <- terra::rast(as.character(background))
  } else { bg_layer <- NULL }
  
  # extract effort layer from raster if provided (note currently uses layers in 
  # existing raster stack, could read in other layers)
  if(is.numeric(effort)){
    eff_layer <- env_extent[[effort]]
  } else if((is.character(effort)|is.factor(effort)) & !grepl("\\.", effort)) {
    eff_layer <- terra::subset(env_extent,as.character(effort))
  } else if ((is.character(effort)|is.factor(effort)) & grepl("\\.", effort)) {
    eff_layer <- terra::rast(as.character(effort))
  } else  {eff_layer <- NULL}
  
  if(is.null(eff_layer)){
    eff_weights <- (env_extent[[1]]*0)+1
  } else if (is.null(bg_layer)){
    eff_weights <- eff_layer/weight_adj
  } else {
    
    # ensure extents match
    if(ext(bg_layer)!=ext(eff_layer)) {
      
      # extend both to match each other
      bg_layer <- extend(bg_layer, eff_layer)
      eff_layer <- extend(eff_layer, bg_layer)
      
    }
    
    eff_weights <- (bg_layer/bg_layer) + (eff_layer/weight_adj)
  }
  
  
  community <- list()
  
  #for each species generate observations
  for (i in 1:n){
    
    # simulate random detection probability from different distributions if specified
    if(detect_prob == "beta") {
      det_prob <- rbeta(1, 2,5)
    } else if(detect_prob == "uniform") {
      det_prob <- runif(1, 0, 1)
    } else if(is.numeric(detect_prob)) {
      det_prob <- detect_prob
    } else {
      stop("!! 'detect_prob' must be a number between 0-1, 'beta' or 'uniform'") # could add option to add custom function?
    }
    
    #subset env raster
    my.stack <- env_extent[[sample(1:nlyr(env_extent),
                                   size = runif(1,env_min,env_max), replace = FALSE)]]
    
    # make edits to specify number of rare vs common species
    #generate a suitability raster
    my.pca.species <- generateSpFromPCA(raster.stack = my.stack, 
                                        sample.points = TRUE, 
                                        nb.points = 10000, 
                                        plot = FALSE, 
                                        niche.breadth = niche_breadth)
    
    #convert to presence-absence
    pa <- convertToPA(my.pca.species, beta = beta, alpha = alpha, plot = FALSE)
    
    #extract prevalence
    prevalence <- as.numeric(pa$PA.conversion[5])
    
    
    #determine maximum number of observations based on prevalence
    #max_obs <- round(prevalence*max_samp)
    max_obs <- max_samp #set max no of observations - could use data?
    
    #sample observations based on bias and detection prob
    occs <- sampleOccurrences(pa, 
                              n = max_obs, 
                              type = "presence-absence", 
                              detection.probability = det_prob, 
                              bias = "manual", 
                              weights = eff_weights, 
                              plot = FALSE)
    
    #rename columns of occurrences data
    if(nrow(occs$sample.points) > 0){ 
      names(occs$sample.points) <- c("lon", "lat", "Real", "Observed")
    }
    
    #subset to PO data
    occs$sample.points <- occs$sample.points[occs$sample.points$Real == 1,]
    occs$sample.points$Observed[occs$sample.points$Observed == 0] <- NA
    
    #store required outputs to list
    community[[i]] <- list(true_prob_occ = wrap(pa$probability.of.occurrence), # need to wrap rasters for saving, then use terra::unwrap
                           pres_abs = wrap(pa$pa.raster), 
                           observations = occs$sample.points, 
                           variables = pa$details$variables, 
                           prevalence = prevalence,
                           detection_probability = det_prob,
                           niche_breadth = niche_breadth)
  }
  
  #return(community)
  
  
  ####   Optional code to sample the same locations across all species in a community
  
  ## sample the same locations across all species in a community
  if(sample_across_species){
    
    ## first, get the sampling locations to be used across all species - biased by butterfly recording effort
    ## use the first community as a template to get cell numbers and coordinates
    
    ## start by getting sampling locations biased by effort
    # convert effort raster to data frame
    eff_df <- as.data.frame(eff_weights, xy = TRUE, na.rm = TRUE)
    
    # get random sample of locations in effort layer - sampling bias only related to effort layer
    # so it's okay to get coordinates from this layer only
    row_ind <- sample(1:nrow(eff_df), size = max_samp, prob = eff_df$layer)
    
    # get the coordinates of sampled cells
    sampled_locs <- eff_df[row_ind,]
    
    # rename columns
    colnames(sampled_locs) <- c("lon", "lat", "layer")
    
    # use the sampled locations to get the corresponding cell numbers in the community file
    # - easier to sample across communities using cell numbers rather than coordinates
    # can't use cell numbers directly from the effort layer, just in case effort layer and community layer 
    # are slightly different extents
    cell_nums <- cellFromXY(terra::unwrap(community[[1]]$pres_abs), xy = sampled_locs[,1:2])
    
    # get presence absence at chosen locations for all species
    comms_sampled <- lapply(community, FUN = function(x) {
      data.frame(sampled_locs[,1:2], 
                 Real = terra::extract(x=terra::unwrap(x$pres_abs), 
                                       y=cell_nums)[,1])
    })
    
    # determine if species is detected - if a species is present at a sampling site,
    # sample between a 0 and 1 according to the detection probability
    comms_observed <- lapply(1:length(comms_sampled), function(com) {
      
      data.frame(comms_sampled[[com]], 
                 Observed = sapply(comms_sampled[[com]]$Real, FUN = function(x) {
                   ifelse(x == 1, sample(c(NA,1), size = 1, 
                                         prob = c(1-community[[com]]$detection_probability, 
                                                  community[[com]]$detection_probability)), NA)
                 }))
      
    })
    
    community2 <- list()
    
    ## bind it all back to the original data format in 'community'
    for(com_n in 1:n){
      
      # get the species list from both sampling methods
      comm_cross_spp <- comms_observed[[com_n]]
      comm_each_spp <- community[[com_n]]
      
      community2[[com_n]] <- list(true_prob_occ = comm_each_spp$true_prob_occ, 
                                  pres_abs = comm_each_spp$pres_abs, 
                                  observations = comm_cross_spp[comm_cross_spp$Real == 1,], 
                                  variables = comm_each_spp$variables,
                                  model_variables = comm_each_spp$model_variables,
                                  prevalence = comm_each_spp$prevalence,
                                  detection_probability = comm_each_spp$detection_probability,
                                  niche_breadth = niche_breadth)
      
      
      
    }
    
    ## replace original community file
    community <- community2
    
  }
  
  # name the community from the parameters
  community_name <- paste0("community_",seed,"_", n, "_sim")
  
  if(!dir.exists(paste0(outPath, community_version_name, simulation_run_name,"/", community_version_name, community_name,"/"))){
    dir.create(paste0(outPath, community_version_name, simulation_run_name,"/", community_version_name, community_name,"/"), recursive = T)
  }
  
  if(write){
    saveRDS(community, file = paste0(outPath, community_version_name, simulation_run_name,"/", community_version_name, community_name,"/", community_version_name, community_name, "_initial.rds"))
  
    print(paste("! Files have been saved, read them using:", 
                paste0("'readRDS(", 
                       paste0(outPath, community_version_name, simulation_run_name,"/", community_version_name, community_name,"/", community_version_name, community_name, "_initial.rds'"))))
    
    }
  
  return(community)
  
  
}
