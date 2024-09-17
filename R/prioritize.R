
#' Calculate taxon conservation benefit
#'
#' Nonlinear functions to convert proportion of range conserved into conservation "benefit."
#'
#' @param x Fraction of taxon range protected (value between 0 and 1).
#' @param lambda Shape parameter.
#'
#' @return Value between 0 and 1.
#' @export
benefit <- function(x, lambda = 1){
      lambda <- 2^lambda
      (1-(1-x)^lambda)^(1/lambda)
}

#' Phylogenetic conservation prioritization
#'
#' Create a ranking of conservation priorities using greedy forward stepwise optimization. This
#'
#' @param sp spatialphy object.
#' @param protection Starting protection status, either a raster or vector, with values between 0 and 1 representing their conservation effectiveness. If this argument is not supplied, it is assumed no existing reserves are present.
#' @param lambda Shape parameter for taxon conservation benefit function.
#' @param level Effectiveness level of proposed new reserves (number between 0 and 1, with same meaning as starting \code{protection}).
#'
#' @return A ranking of conservation priorities, with low values representing higher priorities.
#' @export
sphy_prioritize <- function(sp,
                       protection = NULL,
                       lambda = 1,
                       level = 1){

      e <- sp$tree$edge.length / sum(sp$tree$edge.length) # edges evolutionary value

      if(is.null(protection)){
            p <- rep(0, nrow(sp$occ))
      }else{
            p <- protection[] # protected
      }

      r <- rep(NA, length(p)) # prioritization rankings
      m <- apply(sp$occ, 2, function(x) x / sum(x)) # normalize to fraction of range

      pb <- txtProgressBar(min = 0, max = sum(rowSums(m) > 0), initial = 0, style = 3)
      for(i in 1:length(p)){
            setTxtProgressBar(pb, i)

            # value of current reserve network iteration
            b <- apply(m, 2, function(x) sum(x * p)) # range protection
            v <- sum(e * benefit(b, lambda))

            # marginal value of protecting each cell
            mp <- pmax(0, (level - p)) # marginal protection boost
            u <- apply(m, 2, function(x) x * mp) # unprotected value per taxon*cell
            mv <- apply(u, 1, function(x) sum(e * benefit(x + b, lambda))) - v

            if(min(mv) < 0) stop("bork")

            # protect optimal site
            o <- which(mv == max(mv)) # identify optimal site(s)
            if(length(o) > 1) o <- sample(o, 1) # random tiebreaker
            if(mv[o] == 0) break()
            p[o] <- level # protect site
            r[o] <- i # record ranking
      }
      close(pb)

      # return prioritization
      if(!is.null(sp$spatial)) priorities <- sp$spatial
      priorities[] <- r
      setNames(priorities, "priority")
}
