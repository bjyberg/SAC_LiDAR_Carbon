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
# smoothing). If working with large areas & different tree types,
#I would crop out the species/hedgerows and saving them to different folders
# with the site as the file name, but other options could easily be implemented.
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
library(spatstat)
source('helper_functions.R') #script with the custom functions used
#remotes::install_github("Jean-Romain/lidRplugins") --adds extra seg. algorithms

#Lidar_folder <- '~/lidar' #directory holding the lidar data if processing all
#or 
Lidar_file <- '~/lidar/treed.las' #for processing individual point clouds
output_folder <- '~/lidar/' #Required -- directory for outputs
tree <- 'Angiosperm' #Gymnosperm, angiosperm, or hedgerow_random 
#(hedgerow also available as tree type to test watershed segmentation)
site_name <- tools::file_path_sans_ext(basename(Lidar_file)) #name for outputs
#segmentation_algorithm <- li2012() #Not required -- defaults to dalponte2016()

# lidR uses half available cores by default, can be changed below
# N_cores can be any value (up to the number of cores in available)
# careful to only run this once, as it will use all cores if run again
N_cores <- 18 #get_lidr_threads() *1.5 #changes to .75 available cores (18/24)

#### Code for processing individual point clouds ####
if (exists('Lidar_file')) {
  start_time <- Sys.time()
  set_lidr_threads(N_cores)
  cat(paste('Using', get_lidr_threads(), 'cores'), sep = '\n')
  
  las <- readLAS(Lidar_file)
  
  dtm <- create_dtm(las, classify = F, plot = F, output_path = output_folder,
                    site = site_name)
  
  norm_las <- normalize_height(las, knnidw()) #other algorithms available
  
  chm <- create_chm(norm_las, smooth = T, output_path = output_folder,
                    site = site_name)
  
  if (tree == 'hedgerow_random' | tree == 'Hedgerow_random') {
    hedgerow_agb <- random_hedgerows(chm, point_distance = 2, iterations = 30)
    #point distance and number of iterations should be tuned, and
    #the study using this method performed a sensitivity analysis for distance
  } else { 
    variable <- function(x) {x * 0.1 + 5} #function for variable window size
    Seg_algo <- watershed(chm, th_tree = 5, tol = .2, ext = 1)
    tree_segmentation <- detect_trees(norm_las, chm, window_size = variable,
                                      plot = F,
                                      segmentation_algorithm = Seg_algo)
    agb <- calculate_biomass(tree_segmentation, tree_type = tree, DBH = T,
                             output_path = output_folder, site = site_name)
  }
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
    name <- tools::file_path_sans_ext(basename(i)) #get site from file name
    dtm <- create_dtm(las, classify = F)
    norm_las <- normalize_height(las, knnidw())
    chm <- create_chm(norm_las, smooth = TRUE)
    variable <- function(x) {x * 0.1 + 3}
    tree_segmentation <- detect_trees(norm_las, chm, window_size = variable)
    agb <- calculate_biomass(tree_segmentation, tree_type = tree,
                             output_path = output_folder, site = name)
  })
  
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  cat("Total time:", paste(round(elapsed_time, 2), 'minutes'), sep = '\n')
}
