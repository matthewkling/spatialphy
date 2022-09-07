
#' Calculate spatial phyologenetic diversity metrics
#'
#' @param sp spatialphy object (created by \code{sphy()} or \code{simulate_sphy()}).
#' @param spatial Boolean: should the function return a spatial object (TRUE, default) or a vector (FALSE).
#'
#' @details The function calculates the following metrics:
#' * TR: Terminal richness, i.e. richness of terminal taxa (in many cases these are species)
#' * CR: Clade richness, i.e. richness of taxa at all levels
#' * PD: Phylogenetic diversity
#' * TE: Terminal endemism, i.e. total endemic diversity of terminal taxa, aka WE
#' * CE: Clade endemism, i.e. total endemic diversity of taxa at all levels
#' * PE: Phylogenetic endemism
#' * Em: Mean endemism (derivation is equivalent to E / CR)
#' * PDm: Mean phylogenetic diversity, i.e. branch length of mean resident (derivation is equivalent to PD / CR)
#' * PEm: Mean phylogenetic endemism, i.e. branch length / range size of mean resident (derivation is equivalent to PE / CR
#' * BEm: Mean branch length of the endemics
#'
#' @return A matrix or raster stack with a column or layer (respectively) for each diversity metric.
#' @export
sphy_diversity <- function(sp, spatial = T){

      ## taxon variables ##

      V <- sp$tree$edge.length / sum(sp$tree$edge.length)
      V[V==Inf] <- max(V[V!=Inf])

      R <- apply(sp$occ, 2, sum, na.rm=T) # range sizes

      tips <- tip_indices(sp$tree)


      ## site variables ##

      div <- cbind(TR =  apply(sp$occ[, tips], 1, sum, na.rm=T),
                   CR =  apply(sp$occ, 1, sum, na.rm=T),
                   PD =  apply(sp$occ, 1, function(p) sum(p * V, na.rm = T)),
                   TE =  apply(sp$occ[, tips], 1, function(p) sum(p / R[tips], na.rm = T)),
                   CE =  apply(sp$occ, 1, function(p) sum(p / R, na.rm = T)),
                   PE =  apply(sp$occ, 1, function(p) sum(p * V / R, na.rm = T)),
                   Em =  apply(sp$occ, 1, function(p) weighted.mean(1 / R, w = p, na.rm = T)),
                   PDm = apply(sp$occ, 1, function(p) weighted.mean(V, w = p, na.rm = T)),
                   PEm = apply(sp$occ, 1, function(p) weighted.mean(V / R, w = p, na.rm = T)),
                   BEm = apply(sp$occ, 1, function(p) weighted.mean(V, w = p / R, na.rm = T)) )

      if(spatial & !is.null(sp$spatial)) div <- to_raster(div, sp$spatial)
      return(div)
}


#' Randomization tests including CANAPE
#'
#' This function is a wrapper around the \code{cpr_rand_test} function from Joel Nitta's \code{canaper} package.
#'
#' @param sp spatialphy object
#' @param null_model see \code{canaper::cpr_rand_test}
#' @param spatial Boolean: should the function return a spatial object (TRUE, default) or a vector (FALSE).
#' @param ... further arguments passed to \code{canaper::cpr_rand_test}
#'
#' @details This function runs \code{canaper::cpr_rand_test}; see the help for that function for details.
#'
#' It also runs \code{canaper::cpr_classify_endem} on the result, and includes the resulting classification as an additional variable, 'endem_type', in the output. 'endem_type' values 0-4 correspond to not-significant, neo, paleo, mixed, and super endemesim, respectively.
#'
#' @return A matrix or raster stack with a column or layer (respectively) for each metric.
#' @export
sphy_rand <- function(sp, null_model = "swap", spatial = T, ...){

      cpr <- require("canaper")
      if(!cpr) stop("The sphy_canape() function requires the canaper library, but couldn't find it; please see https://github.com/joelnitta/canaper for info and installation.")

      phy <- sp$tree
      comm <- sp$occ[, tip_indices(phy)]
      colnames(comm) <- phy$tip.label
      rownames(comm) <- paste0("s", 1:nrow(comm))
      cm <- comm[rowSums(comm) > 0, ]

      r <- as.matrix(canaper::cpr_rand_test(cm, phy, null_model = null_model, ...))

      ro <- matrix(NA, nrow(comm), ncol(r))
      ro[rowSums(comm) > 0, ] <- r
      rownames(ro) <- rownames(comm)
      colnames(ro) <- colnames(r)

      ro <- canaper::cpr_classify_endem(as.data.frame(ro))
      ro$endem_type <- as.integer(factor(ro$endem_type,
                                         levels = c("not significant", "neo", "paleo", "mixed", "super"))) - 1
      ro <- as.matrix(ro)

      if(spatial & !is.null(sp$spatial)) ro <- to_raster(ro, sp$spatial)
      return(ro)
}


