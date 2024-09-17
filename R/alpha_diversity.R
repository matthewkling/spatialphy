
#' Calculate spatial phylogenetic diversity and endemism metrics
#'
#' @param sp spatialphy object (created by \code{sphy()} or \code{simulate_sphy()}).
#' @param spatial Boolean: should the function return a spatial object (TRUE, default) or a vector (FALSE).
#'
#' @details The function calculates the following metrics:
#' * TR: Terminal richness, i.e. richness of terminal taxa (in many cases these are species)
#' * CR: Clade richness, i.e. richness of taxa at all levels (equivalent to PD on a cladogram)
#' * PD: Phylogenetic diversity
#' * TE: Terminal endemism, i.e. total endemism-weighted diversity of terminal taxa (a.k.a. "weighted endemism")
#' * CE: Clade endemism, i.e. total endemism-weighted diversity of taxa at all levels (equivalent to PE on a cladrogram)
#' * PE: Phylogenetic endemism, i.e. endemism-weighted PD
#' * Em: Mean endemism (equivalent to CE / CR)
#' * RPD: Relative phylogenetic diversity, i.e. branch length of mean resident (equivalent to PD / CR)
#' * PEm: Mean phylogenetic endemism, i.e. branch length / range size of mean resident (equivalent to PE / CR)
#' * RPE: Relative phylogenetic endemism, i.e. mean endemism-weighted branch length (equivalent to PE / CE)
#'
#' @return A matrix or raster stack with a column or layer (respectively) for each diversity metric.
#' @export
sphy_diversity <- function(sp, spatial = T){

      ## taxon variables ##
      V <- sp$tree$edge.length / sum(sp$tree$edge.length)
      V[V == Inf] <- max(V[V != Inf])
      R <- apply(sp$occ, 2, sum, na.rm = T) # range sizes
      tips <- tip_indices(sp$tree)

      ## site variables ##
      div <- cbind(TR =  apply(sp$occ[, tips], 1, sum, na.rm=T),
                   CR =  apply(sp$occ, 1, sum, na.rm=T),
                   PD =  apply(sp$occ, 1, function(p) sum(p * V, na.rm = T)),
                   TE =  apply(sp$occ[, tips], 1, function(p) sum(p / R[tips], na.rm = T)),
                   CE =  apply(sp$occ, 1, function(p) sum(p / R, na.rm = T)),
                   PE =  apply(sp$occ, 1, function(p) sum(p * V / R, na.rm = T)),
                   Em =  apply(sp$occ, 1, function(p) weighted.mean(1 / R, w = p, na.rm = T)),
                   RPD = apply(sp$occ, 1, function(p) weighted.mean(V, w = p, na.rm = T)),
                   PEm = apply(sp$occ, 1, function(p) weighted.mean(V / R, w = p, na.rm = T)),
                   RPE = apply(sp$occ, 1, function(p) weighted.mean(V, w = p / R, na.rm = T)) )

      if(spatial & !is.null(sp$spatial)) div <- to_raster(div, sp$spatial)
      return(div)
}





