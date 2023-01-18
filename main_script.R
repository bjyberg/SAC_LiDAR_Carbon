library(terra)
library(lidR)
library(sf)
source('helper_functions.R') #script with the custom functions used
#remotes::install_github("Jean-Romain/lidRplugins") --adds extra seg. algorithms 


Lidar_folder <- '~/lidar' #directory holding the lidar data if processing all
Lidar_file <- '~/lidar/treed.las' #for processing individual point clouds
output_path <- '~/lidar' #directory for outputs
tree_type <- 'Angiosperm'

# lidR uses half available cores by default, can be changed below
# N_cores can be any value (up to the number of cores in available)
N_cores <- get_lidr_threads() *1.5 #changes to .75 available cores

####  the code for processing individual point clouds ####
start_time <- Sys.time()
set_lidr_threads(N_cores)
cat(paste('Using', get_lidr_threads(), 'cores'), sep = '\n')

las <- readLAS(Lidar_file)

dtm <- create_dtm(las, classify = T, plot = F)

norm_las <- normalize_height(las, knnidw()) #other algorithms available

chm <- create_chm(norm_las, smooth = TRUE)

variable <- function(x) {x * 0.1 + 3} #function for variable window size
#Seg_algo <- li2012() #there are many options for this in lidR 
tree_segmentation <- detect_trees(norm_las, chm, window_size = variable,
                                  plot = TRUE)

agb <- calculate_biomass(tree_segmentation, tree_type = tree_type)

end_time <- Sys.time()
elapsed_time <- end_time - start_time
cat("Total time:", paste(round(elapsed_time/60, 2), 'minutes'), sep = '\n')


#### the code for processing a full folder of point clouds ####
start_time <- Sys.time()
set_lidr_threads(N_cores)
cat(paste('Using', get_lidr_threads(), 'cores'), sep = '\n')

las_files <- list.files(Lidar_folder, '*.las', #places paths to data in a list
                        recursive = F, #can be True if in sub-folders
                        full.names = T)

all_AGB <- lapply(las_files, function(x) { #function args can still changed
  las <- readLAS(x)
  dtm <- create_dtm(las, classify = F)
  norm_las <- normalize_height(las, knnidw())
  chm <- create_chm(norm_las, smooth = TRUE)
  variable <- function(x) {x * 0.1 + 3}
  tree_segmentation <- detect_trees(norm_las, chm, window_size = variable)
  agb <- calculate_biomass(tree_segmentation, tree_type = tree_type)
  name <- tools::file_path_sans_ext(basename(x))
  assign(paste0("agb_site_", name), agb)
  writeVector(agb, paste0(output_path, name))
})

end_time <- Sys.time()
elapsed_time <- end_time - start_time
cat("Total time:", paste(round(elapsed_time/60, 2), 'minutes'), sep = '\n')

assign(paste0('z', pi), N_cores)
