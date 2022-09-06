
#' @import ape
#' @import raster
#' @import stats
NULL



tip_indices <- function(tree) which(tree$edge[,2] %in% setdiff(tree$edge[,2], tree$edge[,1]))

parentProb <- function(x) 1 - prod(1 - x)

# for a given tree edge index, calculate occurrence probability for each grid cell
build_clade_range <- function(e, phylo, sxt){
      node <- phylo$edge[e,2]
      if(node <= length(phylo$tip.label)){
            otu <- phylo$tip.label[node]
            prob <- sxt[,otu]
      } else{
            clade <- extract.clade(phylo, node)
            otu <- clade$tip.label
            prob <-  apply(sxt[,otu], 1, parentProb)
      }
      return(prob)
}

#' Derive ranges for the internal nodes of a tree
#'
#' @param tree Phylogeny (object of class "phylo").
#' @param tip_occs Site-by-taxon matrix of occurrences for terminal taxa.
#'
#' @return A site-by-taxon matrix with a column for every branch, including terminals and internal edges.
#' @export
build_tree_ranges <- function(tree, tip_occs){
      sapply(1:nrow(tree$edge), build_clade_range, phylo=tree, sxt=tip_occs)
}




#' Convert a site-by-variable matrix into a raster brick
#'
#' @param m Matrix.
#' @param template Raster layer, with number of cells equal to the number of rows in m.
#'
#' @return A raster brick.
#' @export
to_raster <- function(m, template){
      a <- array(m,
                 c(ncol(template), nrow(template), ncol(m)),
                 list(NULL, NULL, colnames(m)))
      brick(aperm(a, c(2, 1, 3)))
}




#' Simulate a toy spatial phylogenetic dataset
#'
#' @param n_tips Number of terminals on phylogeny.
#' @param n_x Number of raster cells in x dimension of landscape.
#' @param n_y Number of raster cells in y dimension of landscape.
#' @param boolean_ranges Should simulated ranges be boolean (T) or probabilities (F, default)?
#'
#' @return A list containing a phylogenetic tree and a raster stack of terminal taxon ranges.
#' @export
sphy_simulate <- function(n_tips = 10,
                          n_x = 20,
                          n_y = 20,
                          boolean_ranges = F
){
      tree <- ape::rtree(n_tips)
      tip_ranges <- stack(lapply(1:n_tips, function(i){
            # distance decay from range center
            x <- sqrt(matrix(rep(1:n_x, n_y) - sample(n_x, 1), n_x, n_y) ^ 2 +
                            matrix(rep(1:n_y, each = n_x) - sample(n_y, 1), n_x, n_y) ^ 2)
            r <- max(c(n_x, n_y)) / runif(1, 1, 5) # range radius
            x[x > r] <- r
            p <- runif(1, .1, 1) # prevalence at range center
            raster(t((max(x) - x) / max(x))) * p
      }))
      names(tip_ranges) <- tree$tip.label
      if(boolean_ranges) tip_ranges <- tip_ranges > 0

      return(sphy(tree, tip_ranges))

      # return(list(tree = tree,
      #             tip_ranges = tip_ranges))
}



#' Create a spatial phylogeny object
#'
#' @param tree Phylogeny of class "phylo".
#' @param occ Occurrence data across taxa and sites. Can be a matrix with a column per taxon and a row per site, or a RasterBrick/RasterStack with a layer per taxon. There must either be one taxon per edge in the tree, or one taxon per tip.
#' @param spatial An optional raster layer indicating site locations. Ignored if \code{occ} is a raster object.
#'
#' @return A "spatialphy" object
#' @export
sphy <- function(tree, occ, spatial = NULL){

      if(!inherits(tree, "phylo")) stop("tree must be an object of class 'phylo'")
      if(!inherits(occ, c("matrix", "RasterBrick", "RasterStack"))) stop("occurrence data must be a matrix or a raster stack/brick")
      if(!is.null(spatial) & !inherits(spatial, "RasterLayer")) stop("spatial reference must be a RasterLayer")
      if(inherits(spatial, "RasterLayer")) if(ncell(spatial) != nrow(occ)) stop("raster has different number of sites from occurrence data")

      if(inherits(occ, c("RasterBrick", "RasterStack"))){
            spatial <- occ[[1]]
            spatial[] <- NA
            occ <- values(occ)
      }

      if(! ncol(occ) %in% c(length(tree$tip.label), length(tree$edge.length))){
            stop("occurrence dataset does not match tree: must have exactly one variable per terminal, or per edge")
      }

      if(ncol(occ) == length(tree$tip.label)){
            occ <- build_tree_ranges(tree, occ)
      }

      structure(list(tree = tree,
                     occ = occ,
                     spatial = spatial,
                     dist = NULL),
                class = "spatialphy")
}
