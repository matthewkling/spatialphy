
phylo_sorensen <- function(b1, b2, b12) 2 * b12 / (b1 + b2)

#' Compute a pairwise community phylogenetic similarity matrix
#'
#' @param M Site-by-taxon matrix.
#' @param branch_lengths Vector of branch lenghts for each taxon.
#'
#' @return Pairwise matrix of sorensen similarity between sites
#' @export
beta_diversity <- function(M, branch_lengths){

      M_bl <- t(apply(M, 1, function(x) x * branch_lengths))

      n <- nrow(M)
      beta <- matrix(NA, n, n)

      for(i in 1:n){
            for(j in 1:n){
                  if(j > i) next()
                  beta[i, j] <- beta[j, i] <- phylo_sorensen(b1 = sum(pmax(M_bl[i,] - M_bl[j,], 0)),
                                                             b2 = sum(pmax(M_bl[j,] - M_bl[i,], 0)),
                                                             b12 = sum(pmin(M_bl[i,], M_bl[j,])))
            }
      }

      return(beta)

      # image(matrix(1/beta[sample(n, 1),], n^.5, n^.5))
}




#' Compute a phyologenetic turnover matrix
#'
#' @param sp spatialphy object.
#' @param endemism Logical indicating whether occurrence values should be divided by column (species) totals.
#' @param normalize Logical indicating whether occurrence values should be divided by row (community) totals.
#'    If so, this happens after endemism division.
#' @param add Logical indicating whether results should be added to the spatialphy object and returned (TRUE) or returned alone (FALSE)?
#'
#' @return A pairwise phylogenetic difference matrix, with values of 1/PhyloSor, either on its own or as an element of \code{sp}.
#' @export
sphy_dist <- function(sp, endemism = FALSE, normalize = TRUE, add = TRUE){

      if(!is.null(sp$dist)){
            message("distance already included in dataset; skipping calculation")
            if(add) return(sp)
            if(!add) return(sp$dist)
      }

      occ <- sp$occ
      if(endemism) occ <- apply(occ, 2, function(x) x / sum(x))
      if(normalize) occ <- t(apply(occ, 1, function(x) x / sum(x)))
      occ[!is.finite(occ)] <- 0
      dist <- 1 / beta_diversity(occ, sp$tree$edge.length / sum(sp$tree$edge.length))
      dist = as.dist(dist)

      if(add){
            sp$dist <- dist
            return(sp)
      }else{
            return(dist)
      }
}


#' Phyologenetic regionalization
#'
#' @param sp spatialphy object.
#' @param k Number of spatial clusters to divide the region into (Positive integer).
#' @param method Clustering method. Options include "kmeans", and the methods listed under \link[stats]{hclust}.
#' @param endemism Logical indicating whether occurrence values should be divided by column (species) totals.
#' @param normalize Logical indicating whether occurrence values should be divided by row (community) totals.
#'    If so, this happens after endemism division.
#'
#' @return A raster or matrix with an integer indicating which of the \code{k} regions each site belongs to.
#' @export
sphy_regions <- function(sp, k = 5, method = "kmeans", endemism = FALSE, normalize = TRUE){

      # sites with taxa
      a <- rowSums(sp$occ) > 0

      if(method == "kmeans"){
            occ <- sp$occ[a,]
            if(endemism) occ <- apply(occ, 2, function(x) x / sum(x))
            if(normalize) occ <- t(apply(occ, 1, function(x) x / sum(x)))
            occ[!is.finite(occ)] <- 0
            regions <- kmeans(occ, k)$cluster
      }else{
            if(is.null(sp$dist)) sp <- sphy_dist(sp, endemism, normalize, add = T)

            d <- as.matrix(sp$dist)
            rownames(d) <- colnames(d) <- paste("cell", 1:ncol(d))
            da <- d[a, a]

            # sites fully segregated by the 2 basal clades have Inf distance;
            # set distance to value greater than max observed distance
            da[is.infinite(da)] <- max(da[!is.infinite(da)]) + 1000
            da <- as.dist(da)

            clust <- hclust(da, method = method)
            regions <- cutree(clust, k)
      }


      a[a] <- regions
      a[!a] <- NA

      if(!is.null(sp$spatial)){
            r <- sp$spatial
            r[] <- a
            return(setNames(r, "phyloregion"))
      }else{
            return(a)
      }

}



#' Color space ordination of spatial phylogenetic composition
#'
#' @param sp spatialphy object.
#' @param method Ordination method, either "pca" (principal component analysis), "cmds" (classical MDS),
#'    or "nmds" (nonmetric MDS, the default, which is slower but often preferred).
#' @param trans A function giving the transformation to apply to each dimension of the ordinated data.
#'    The default is the identity funciton. Specifying \code{rank} generates a more uniform color distribution.
#'
#' @return A raster or matrix with layers or columns (respectively) containing RGB color values in the range 0-1.
#' @export
sphy_rgb <- function(sp, method = "nmds", trans = identity){

      d <- as.matrix(sp$dist)
      rownames(d) <- colnames(d) <- paste("cell", 1:ncol(d))
      a <- !is.na(rowSums(d)) # sites with taxa
      da <- d[a, a]

      # sites fully segregated by the 2 basal clades have Inf distance;
      # set distance to value greater than max observed distance
      da[is.infinite(da)] <- max(da[!is.infinite(da)]) + 1000

      # ordinate
      if(method == "cmds") ord <- cmdscale(da, k = 3)
      if(method == "nmds") ord <- vegan::metaMDS(da, k = 3, trace = 0)$points
      if(method == "pca") ord <- prcomp(sp$occ[a,], scale. = T)$x[,1:3]

      # scale
      rgb <- apply(ord, 2, function(x){
            x <- trans(x)
            (x - min(x)) / (max(x) - min(x))
      })

      # reinsert NA values
      b <- cbind(a, a, a)
      b[a,] <- rgb
      b[!a,] <- NA

      if(!is.null(sp$spatial)){
            # convert to raster
            r <- stack(sp$spatial, sp$spatial, sp$spatial)
            r[] <- b
            r <- setNames(r, c("r", "g", "b"))
            return(r)
      }else{
            return(b)
      }
}


