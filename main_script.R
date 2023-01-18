## ######################################################################### ##
# This script can be used to run above ground biomass estimates on a LiDAR 
# point cloud. The script can be run solely by changing the options within
# the set-up section to the desired files, the segmentation algorithm, etc.
# Changing the algorithms used for DTM and CHM creation, tree top detection,
# etc. will require changes in the 2 processing sections. Similarly, it is
# possible to save all of these individual outputs to a folder through basic
# modifications to the script. It should be relatively easy to implement these 
# by reading and changing this script and the 'helper_functions.R' script. I
# would suggest testing a few different algorithms for segmentation and tree
# detections, as well as different internal parameters (i.e., window size,
# smoothing). The renv files are to ensure reproducibility by managing
# packages and dependencies. 
#
#
# This script was developed as part of a larger soil carbon project at
# SAC consulting which was funded by the Scottish Government.
# Author: Brayden Youngberg 
# Contact: brayden.youngberg@sruc.ac.uk
## ######################################################################### ##

#### Set-up ####
#the full script can be run just by changing these to the required data
library(terra)
library(lidR)
library(sf)
source('helper_functions.R') #script with the custom functions used
#remotes::install_github("Jean-Romain/lidRplugins") --adds extra seg. algorithms

Lidar_folder <- '~/lidar' #directory holding the lidar data if processing all
#or 
Lidar_file <- '~/lidar/treed.las' #for processing individual point clouds

output_folder <- '~/lidar/' #Required -- directory for outputs
tree <- 'Angiosperm' #or Gymnosperm
#segmentation_algorithm <- li2012() #Not required -- defaults to dalponte2016()

# lidR uses half available cores by default, can be changed below
# N_cores can be any value (up to the number of cores in available)
# careful to only run this once, as it will use all cores if run again
N_cores <- get_lidr_threads() *1.5 #changes to .75 available cores

#### Code for processing individual point clouds ####
if (exists('Lidar_file')) {
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
  
  agb <- calculate_biomass(tree_segmentation, tree_type = tree)
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  cat("Total time:", paste(round(elapsed_time, 2), 'minutes'), sep = '\n')
}

#### Code for processing a full folder of point clouds ####
if (exists('Lidar_folder')) {
  start_time <- Sys.time()
  set_lidr_threads(N_cores)
  cat(paste('Using', get_lidr_threads(), 'cores'), sep = '\n')
  
  las_files <- list.files(Lidar_folder, '*.las', #places paths to data in a list
                          recursive = F, #can be True if in sub-folders
                          full.names = T)
  
  all_AGB <- lapply(las_files, function(i) { #function args can still be changed
    las <- readLAS(i)
    dtm <- create_dtm(las, classify = F)
    norm_las <- normalize_height(las, knnidw())
    chm <- create_chm(norm_las, smooth = TRUE)
    variable <- function(x) {x * 0.1 + 3}
    tree_segmentation <- detect_trees(norm_las, chm, window_size = variable)
    name <- tools::file_path_sans_ext(basename(i)) #get site from file name
    agb <- calculate_biomass(tree_segmentation, tree_type = tree,
                             output_path = output_folder, site = name)
  })
  
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  cat("Total time:", paste(round(elapsed_time, 2), 'minutes'), sep = '\n')
}
