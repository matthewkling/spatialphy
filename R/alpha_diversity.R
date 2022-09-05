
#' Calculate spatial phyologenetic diversity metrics
#'
#' @param sp spatialphy object (created by \code{sphy()} or \code{simulate_sphy()}).
#' @param spatial Boolean: should the function return a spatial object (TRUE, default) or a vector (FALSE).
#'
#' @return A matrix or raster stack with a column or layer (respectively) for each diversity metric.
#' @export
diversity <- function(sp, spatial = T){

      ## taxon variables ##

      V <- sp$tree$edge.length / sum(sp$tree$edge.length)
      V[V==Inf] <- max(V[V!=Inf])

      R <- apply(sp$occ, 2, sum, na.rm=T) # range sizes


      ## site variables ##

      # CR: Clade richness, i.e. richness of taxa at all levels
      # PD: Phylogenetic diversity
      # E: Endemism, i.e. total endemic diversity, aka WE
      # PE: Phylogenetic endemism
      # Em: Mean endemism (derivation is equivalent to E / CR)
      # PDm: Mean phylogenetic diversity, i.e. branch length of mean resident (derivation is equivalent to PD / CR)
      # PEm: Mean phylogenetic endemism, i.e. branch length / range size of mean resident (derivation is equivalent to PE / CR
      # BEm: Mean branch length of the endemics

      div <- cbind(CR =  apply(sp$occ, 1, sum, na.rm=T),
                   PD =  apply(sp$occ, 1, function(p) sum(p * V, na.rm = T)),
                   E =   apply(sp$occ, 1, function(p) sum(p / R, na.rm = T)),
                   PE =  apply(sp$occ, 1, function(p) sum(p * V / R, na.rm = T)),
                   Em =  apply(sp$occ, 1, function(p) weighted.mean(1 / R, w = p, na.rm = T)),
                   PDm = apply(sp$occ, 1, function(p) weighted.mean(V, w = p, na.rm = T)),
                   PEm = apply(sp$occ, 1, function(p) weighted.mean(V / R, w = p, na.rm = T)),
                   BEm = apply(sp$occ, 1, function(p) weighted.mean(V, w = p / R, na.rm = T)) )

      if(spatial & !is.null(sp$spatial)) div <- to_raster(div, sp$spatial)
      return(div)
}
