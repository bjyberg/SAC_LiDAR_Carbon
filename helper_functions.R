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
library(spatstat)


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
                       smooth = FALSE, plot = FALSE) {
  chm <- rasterize_canopy(point_cloud, res = .5, pitfree(subcircle = 0.15))
  if (isTRUE(smooth)) {
    kernel <- matrix(1,5,5) #can easily be changed to different size
    chm <- terra::focal(chm, w = kernel, median, na.rm = TRUE)
  }
  if (!missing(output_path)) {
    tryCatch(
      {
        writeRaster(dtm, paste0(output_path, site, '_chm.tif'))
      },
      error = function(e) {
        message('An error occured while saving the chm:')
        print(e)
        print('outputs not saved to file, check object instead')
      },
      warning = function(w) {
        message('An warning occured while saving the chm:')
        print(w)
        print('outputs not saved to file, check object instead')
      }
    )
  }
  if (isTRUE(plot)) {
    plot(chm)
  }
  return(chm)
}

detect_trees <- function(point_cloud, chm, output_path, site,
                         segmentation_algorithm, min_tree = 5,
                         window_size = 3, save_las = FALSE,
                         plot = FALSE) {
  window <- window_size #this can be a fixed size or variable
  t_tops <- locate_trees(chm, lmf(ws = window)) 
  if (missing(segmentation_algorithm)) {
    t_tops <- locate_trees(point_cloud, lmf(ws = window)) 
    algo <- dalponte2016(chm, t_tops, th_tree = min_tree) #default
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
        if (isTRUE(save_las)) {
          writeLAS(point_cloud, paste0(output_path, site, 'tree_segmented.las'))
        }
        write_sf(crowns, paste0(output_path, site, '_crowns.gpkg'))
      },
      error = function(e) {
          message('An error occured while saving the segmented trees:')
          print(e)
          print('outputs not saved to file, check object instead')
      },
      warning = function(w) {
          message('An warning occured while saving the segmented trees:')
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
  
calculate_biomass <- function(crown_polygons, tree_type, output_path, site,
                              DBH = TRUE){
  tree_z <- crown_polygons[['Z']]
  c_area <- crown_polygons[['convhull_area']]
  c_diameter <- sqrt(c_area / pi)
  crown_polygons$Crown_diameter <- c_diameter
  if (isTRUE(DBH)) {
    crown_polygons$Tree_DBH <- (0.557 * (tree_z * c_diameter)^0.809
                          * exp((0.056^2) / 2))
    Tree_DBH <- crown_polygons[['Tree_DBH']]
    DBH_summary <- summary(crown_polygons$Tree_DBH)
    # equation from Jucker, T., et al. (2016) Allometric equations for
    # integrating remote sensing imagery into forest monitoring programmes.
    # Global Change Biology. 23(1). 177-190
    # https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.13388
    p <- 0.525 # Wood Density
    crown_polygons$DBH_Biomass <- (0.0673
    * (p * (Tree_DBH^2) * tree_z)^0.976
      * exp((0.367^2) / 2))
    # NEEDS CITED
    DBH_Biomass_summary <- summary(crown_polygons$DBH_Biomass)
    DBH_AGB_total <- sum(crown_polygons$DBH_Biomass)
    DBH_agb_txt <- c(paste('The total DBH calculated AGB for',
                  tree_type, 'on the site is:'),
                  paste(DBH_AGB_total, 'kg.')) # Make sure unit is correct... 
  }
  if (tree_type == 'angiosperm' | tree_type == 'Angiosperm') {
    a <- 0
    b <- 0
    crown_polygons$AGB_pred <- ((0.016 + a) 
                                * (tree_z * c_diameter)^(2.013 + b)
                                * exp(0.204^2 / 2))
  } else if (tree_type == 'gymnosperm' | tree_type == 'Gymnosperm') {
    a <- 0.093
    b <- -0.223
    crown_polygons$AGB_pred <- ((0.016 + a) 
                                * (tree_z * c_diameter)^(2.013 + b)
                                * exp(0.204^2 / 2))
    #equation from Jucker, T., et al. (2016) Allometric equations for 
    #integrating remote sensing imagery into forest monitoring programmes. 
    #Global Change Biology. 23(1). 177-190
    #https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.13388
    
  } else if (tree_type == 'hedgerow' | tree_type == 'Hedgerow') {
    crown_polygons$AGB_pred <- (0.179 * tree_z^3.3)
    #equation from Black, K., et al. (2014) Carbon Sequestration by Hedgerows in
    #the Irish Landscape.
    #Ireland Environmental Protection Agency, CCRP Report
  } else {
    warning("tree_type must be 'angiosperm', 'gymnosperm', 'hedgerow'")
  }
  
  total_AGB <- sum(crown_polygons$AGB_pred)
  AGB_summary <- summary(crown_polygons$AGB_pred)
  z_summary <- summary(crown_polygons$Z)
  cd_summary <- summary(crown_polygons$Crown_diameter)
  t_agb_txt <- c(paste('The total AGB for', tree_type, 'on the site is:'),
                paste(total_AGB, 'kg.'))

  txt <- capture.output(
    cat(sep = '\n'),
    if (exists("DBH_summary")) {
      cat("Summary Statistics for DBH:", sep = "\n")
      print(DBH_summary)
      cat("Summary Statistics for DBH-based AGB:", sep = "\n")
      print(DBH_Biomass_summary)
      cat(DBH_agb_txt, sep = '\n')
    },
    cat(sep = '\n'),
    cat('Summary Statistics for height:', sep = '\n'),
    print(z_summary),
    cat('Summary Statistics for crown diameter:', sep = '\n'),
    print(cd_summary),
    cat('Summary Statistics for AGB:', sep = '\n'),
    print(AGB_summary),
    cat(t_agb_txt, sep = '\n'),
    cat(paste(nrow(crown_polygons), tree_type, 'trees identified on site')))
  
  if (!missing(output_path)) {
    tryCatch( #added function so the script doesn't stop with a write error
      {
        write_sf(crown_polygons, paste0(output_path, site, '_', tree_type,
                                        '_AGBcrowns.gpkg'))
        writeLines(txt, paste0(output_path, site, '_', tree_type,
                                        '_AGBcrowns.txt'))
      },
      error = function(e) {
        message('An error occured while saving the AGB tree vector:')
        print(e)
        print('outputs not saved to file, check object instead')
      },
      warning = function(w) {
        message('An warning occured while saving the AGB tree vector:')
        print(w)
        print('outputs not saved to file, check object instead')
      }
    )
  }
  
  cat(txt, sep = '\n')
  
  return(crown_polygons)
}

random_hedgerows <- function(chm, point_distance, iterations,
                             output_path, site) {
  clean_chm <- clamp(chm, lower = 1.5, values = F)
  chm_extent <- as.polygons(clean_chm > -Inf)
  hedge_agb_iterations <- data.frame()
  for (i in 1:iterations) {
    points <- st_sample(st_as_sf(chm_extent), type = 'SSI', r = point_distance)
    random_heights <- extract(chm, points, ID = F)
    agb <- (0.179 * random_heights^3.3)
    site_agb <- sum(agb)
    hedge_agb_iterations <- rbind(hedge_agb_iterations, site_agb)
    cat(paste("Iteration", i, 'complete -', iterations - i, 'remain'), sep = '\n')
  }
  names(hedge_agb_iterations) <- 'site_agb'
  hedge_agb_stats <- summary(hedge_agb_iterations)
  hedge_txt <- c('Random hedgerow method', '\n', 'Summary for hedgerow agb:',
                 '\n', hedge_agb_stats)
  cat(hedge_txt, sep = '\n')
  if (!missing(output_path)) {
    tryCatch( #added function so the script doesn't stop with a write error
      {
        writeLines(hedge_agb_stats, paste0(output_path, site, '_', tree_type,
                                        'Hedge_AGB_Summary.txt'))
        write.csv(hedge_agb_iterations, paste0(output_path, site,
                                         'Hedge_AGB_iterations.csv'))
      },
      error = function(e) {
        message('An error occured while saving the random hedge output:')
        print(e)
        print('not saved to file, check object instead')
      },
      warning = function(w) {
        message('An warning occured while saving therandom hedge output:')
        print(w)
        print('not saved to file, check object instead')
      }
    )
  return(hedge_agb_iterations)
  }
}

