## ######################################################################### ##
# Below is a set of functions which are used when running main_script.R. The 
# parameters for these functions can be changed relatively easy, and I would
# suggest changing the algorithms used by the lidR package functions and their
# related parameters (specifically for the tree segmentation, which can be done
# in main_script.R, rather than here), to test which is best for different sites
# and datasets. Other possible changes would be the algorithm for tree
# detection, the resolution and smoothing of the chm, and the window sizes.
#
#
# This script was developed as part of a larger soil carbon project at
# SAC consulting which was funded by the Scottish Government.
# Author: Brayden Youngberg 
# Contact: brayden.youngberg@sruc.ac.uk
## ######################################################################### ##
library(terra)
library(lidR)
library(sf)


create_dtm <- function(point_cloud, output_path, site,
                       classify = FALSE, plot = FALSE) {
  if (isTRUE(classify)) {
    #if true, this adds a ground class into the existing point cloud,
    #or overwrites the current ground class if it was already classified
    ground_alg <- csf() #can be pmf() or mcc(), see lidR package for details
    point_cloud <- classify_ground(point_cloud, algorithm = ground_alg)
  }
  dtm_alg <- knnidw(k = 10L, p = 2) #can be kreiging or tin(), see lidR package
  dtm <- rasterize_terrain(point_cloud, algorithm = dtm_alg)
  if (!missing(output_path)) {
    writeRaster(dtm, paste0(output_path, site, '_dtm.tif'))
  }
  if (isTRUE(plot)) {
    plot_dtm3d(dtm)
  }
  return(dtm)
}

create_chm <- function(point_cloud, output_path, site,
                       smooth = FALSE) {
  chm <- rasterize_canopy(point_cloud, res = .5, pitfree(subcircle = 0.15))
  if (isTRUE(smooth)) {
    kernel <- matrix(1,3,3) #can easily be changed to different size
    chm <- terra::focal(chm, w = kernel, median, na.rm = TRUE)
  }
  if (!missing(output_path)) {
    writeRaster(dtm, paste0(output_path, site, '_dtm.tif'))
  }
  return(chm)
}

detect_trees <- function(point_cloud, chm, output_path, site,
                         segmentation_algorithm,
                         window_size = 3,
                         plot = FALSE) {
  window <- window_size #this can be a fixed size or variable
  t_tops <- locate_trees(chm, lmf(ws = window)) 
  if (missing(segmentation_algorithm)) {
    t_tops <- locate_trees(point_cloud, lmf(ws = window)) 
    algo <- dalponte2016(chm, t_tops) #default, if not provided
  } else {
    algo <- segmentation_algorithm
  }
  point_cloud <- segment_trees(point_cloud, algo, attribute = 'tree_id')
  crowns <- crown_metrics(point_cloud, func = .stdtreemetrics,
                          geom = 'convex', 
                          attribute = 'tree_id')
  if (!missing(output_path)) {
    tryCatch( #added function so the script doesn't stop with a write error
      {
        write_sf(crowns, paste0(output_path, site, '_crowns.gpkg'))
      } error = function(e) {
        message('An error occured while saving the segmented tree polygons:')
        print(e)
        print('outputs not saved to file, check object instead')
      } warning = function(w) {
        message('An warning occured while saving the segmented tree polygons:')
        print(w)
      }
    )
  }
  if (isTRUE(plot)) {
    plot(chm)
    plot(crowns['convhull_area'], alpha = 0.5, add = TRUE)
  }
  
  return(crowns)
}
  
calculate_biomass <- function(crown_polygons, tree_type, output_path, site){
  tree_z <- crown_polygons[['Z']]
  c_area <- crown_polygons[['convhull_area']]
  c_diameter <- sqrt(c_area / pi)
  crown_polygons$Crown_diameter <- c_diameter
  if (tree_type == 'angiosperm' | tree_type == 'Angiosperm') {
    a <- 0
    b <- 0
  } else if (tree_type == 'gymnosperm' | tree_type == 'Gymnosperm') {
    a <- 0.093
    b <- -0.223
  } else {
    warning("tree_type must be either 'angiosperm' or 'gymnosperm'")
  }
  
  crown_polygons$AGB_pred <- ((0.016 + a) 
                                 * (tree_z * c_diameter)^(2.013 + b)
                                 * exp(0.204^2 / 2))
  #equation from Jucker, T., et al. (2016) Allometric equations for integrating 
  #remote sensing imagery into forest monitoring programmes. 
  #Global Change Biology. 23(1). 177-190
  #https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.13388
  
  total_AGB <- sum(crown_polygons$AGB_pred)
  AGB_summary <- summary(crown_polygons$AGB_pred)
  z_summary <- summary(crown_polygons$Z)
  cd_summary <- summary(crown_polygons$Crown_diameter)
  
  if (!missing(output_path)) {
    tryCatch( #added function so the script doesn't stop with a write error
      {
        write_sf(crown_polygons, paste0(output_path, site, '_AGBcrowns.gpkg'))
      } error = function(e) {
        message('An error occured while saving the segmented tree polygons:')
        print(e)
        print('outputs not saved to file, check object instead')
      } warning = function(w) {
        message('An warning occured while saving the segmented tree polygons:')
        print(w)
        print('outputs not saved to file, check object instead')
      }
    )
  }
  cat(sep = '\n')
  cat('Summary Statistics for height:', sep = '\n')
  print(z_summary)
  cat('Summary Statistics for crown diameter:', sep = '\n')
  print(cd_summary)
  cat('Summary Statistics for AGB:', sep = '\n')
  print(AGB_summary)
  cat(paste('The total AGB for the site is:', total_AGB, 'kg.'), sep = '\n')
  return(crown_polygons)
}


