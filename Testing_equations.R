DBH_bio <- function(p, Tree_DBH, tree_z) {
  bio <- (0.0673
          * (p * (Tree_DBH^2) * tree_z)^0.976
          * exp(0.367^2) / 2)
  return(bio)
}

bio <- function(tree_z, c_diameter, a = 0.093, b = -0.223) {
  bio <- ((0.016 + a)
          * (tree_z * c_diameter)^(2.013 + b)
          * exp(0.204^2 / 2))
  return(bio)
}

spruce <- function(tree_z) {
  bio <- 1.32 * (tree_z^1.7) * 1.38
  return(bio)
}
spruce2 <- function(tree_z, dbh) {
  bio <- .23 * (dbh^2.12) + (5*10^-7) * (tree_z^4.99)
  return(bio)
}

hedge <- function(tree_z) {
  bio <- (0.179 * tree_z^3.3)
  return(bio)
}

conif <- function(tree_z) {
  bio_z <- (0.005 * tree_z^1.58 * 1.12)
  return(bio_z)
}

conif_DBH <- function(tree_z, dbh) {
  bio_dbh <- (0.022 * dbh^2.73 + 0.19 * tree_z^2.06)
  return(bio_dbh)
}

conif_DBH(10,40)
broad <- function(tree_z) {
  bio_z <- 0.031 * tree_z^1.72
  return(bio_z)
}
broad_dbh <- function(dbh) {
  bio_dbh <- 0.08 + ((25000 * dbh^2.5)/(dbh^2.5 + 246872))
  return(bio_dbh)
}


calc_dbh <- function(tree_z, c_diameter) {
  (0.557 * (tree_z * c_diameter)^0.809
  * exp((0.056^2) / 2))
}

calc_dbh(6.8, 5)
