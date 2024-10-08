
rm(list = ls())

# deleting species models to free up space on lotus

# list all files in the v4species_models folders
list.files(Sys.glob("outputs/communities/v1narrow_breadth_uniform_detect_community/*/v1species_models/"), 
           pattern = "_AS_", # find the _AS_ pattern - only files I want deleted 
           full.names = TRUE)

# test
unlink(list.files(Sys.glob("outputs/communities/v1narrow_breadth_uniform_detect_community/*/v1species_models/"), 
                  pattern = "_AS_", # find the _AS_ pattern - only files I want deleted 
                  full.names = TRUE)[1:2])

# do it
unlink(list.files(Sys.glob("outputs/communities/v1narrow_breadth_uniform_detect_community/*/v1species_models/"), 
                  pattern = "_AS_", # find the _AS_ pattern - only files I want deleted 
                  full.names = TRUE))

