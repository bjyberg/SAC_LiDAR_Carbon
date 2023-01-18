library(terra)
library(lidR)
library(sf)
source('helper_functions.R') #script with the custom functions used
#remotes::install_github("Jean-Romain/lidRplugins") --adds extra seg. algorithms 



#make sure parallel is set/specify cores
#lidR uses half available cores by default, can be changed below
#N_cores can be any value (up to the number of cores in avaliable)
N_cores <- get_lidr_threads() *1.5 #changes to .75 available cores


Lidar_folder <- '' #directory holding the lidar data if processing all
Lidar_file <- '~/lidar/treed.las' #for processing individual point clouds
output_path <- '' #directory for outputs






####  the code for processing individual point clouds ####
start_time <- Sys.time()
set_lidr_threads(N_cores)
cat(paste('Using', get_lidr_threads(), 'cores'), sep = '/n')

las <- readLAS(Lidar_file)

dtm <- create_dtm(las, classify = T, plot = F)

norm_las <- normalize_height(las, knnidw()) #other algorithms available

chm <- create_chm(norm_las, smooth = TRUE)

variable <- function(x) {x * 0.1 + 3} #function for variable window size
#Seg_algo <- li2012() #there are many options for this in lidR 
tree_segmentation <- detect_trees(norm_las, chm, window_size = variable,
                                  plot = TRUE)

plot(tree_segmentation['convhull_area'])

agb <- calculate_biomass(


#### the code for processing a full folder of point clouds ####
las_files <- list.files(Lidar_folder, '*.las', #places paths to data in a list
                        recursive = T,
                        full.names = T)

#add crop function?
norm_las <- normalize_height(las)

chm <- rasterize_canopy(las, res = 0.5, pitfree(subcircle = 0.15))
