


#' Stratified randomization of community matrix
#'
#' This is a community null model method for quantitative data (e.g. abundance or occurrence probability)
#' that preserves row and column totals, and also approximately preserves the marginal distributions of
#' rows and columns. For each randomization, the data set is split into strata representing numerical ranges
#' of the input quantities, a separate binary randomization is done for each stratum, and the results are
#' combined to produce a randomized, quantitative community matrix. See `vegan::commsim()` for details about
#' other binary and quantitative null models.
#'
#' @param x Community matrix with species in rows, sites in columns, and nonnegative occurrence quantities in cells.
#' @param n_strata Integer giving the number of strata to split the data into. Must be 2 or greater. Larger values
#'    will result in randomizations with less mixing but higher fidelity to marginal distributions.
#' @param transform A function used to transform the values in \code{x} before assigning them to \code{n_strata}
#'    equal intervals. Examples include \code{sqrt}, \code{log}, \code{rank}, etc.; the default is \code{identity}.
#' @param algorithm Null model algorithm, passed to `vegan::nullmodel`. Testing has only been done with the
#'    "curveball" algorithm, so other options should be use with caution. Only binary methods should be used.
#' @param n_iter Positive integer giving the number of iterations to use for stepwise algorithms.
#' @param jitter Number between 0 and 1, indicating how much to randomly jitter the location of stratum boundaries.
#' @param priority Either "rows", "cols", or "neither", indicating whether randomization within strata should
#'    prioritize maintaining the marginal distributions of the rows or columns of the input matrix. The default,
#'    "neither", doesn't give precedence to either dimension.
#' @return A randomized version of \code{x}.
#' @export
strand <- function(x, n_strata = 5, transform = identity, algorithm = "curveball",
                   jitter = .99, priority = "neither", n_iter = 10000, ...){

      nc <- ncol(x)
      nr <- nrow(x)

      # convert to stratified binary community array
      s <- transform(x)
      bw <- diff(range(s)) / n_strata
      breaks <- c(-Inf, seq(min(s)+bw, max(s)-bw, bw), Inf)
      if(jitter > 0){
            offset <- seq(-jitter, jitter, length.out = 1000)
            offset <- sample(offset, 1, prob = jitter - abs(offset))
            breaks <- breaks + offset * bw
      }
      s[] <- as.integer(cut(s, breaks))
      b <- apply(s, 1:2, function(x) replace(rep(0, n_strata), x, 1))

      # quantities shuffled within strata, and within rows or columns or neither
      r <- s
      resample <- function(x, ...) x[sample.int(length(x), ...)]
      if(priority == "rows") for(i in 1:n_strata) for(j in 1:nr) r[j, s[j,] == i] <- resample(x[j, s[j,] == i])
      if(priority == "cols") for(i in 1:n_strata) for(j in 1:nc) r[s[,j] == i, j] <- resample(x[s[,j] == i, j])
      if(priority == "neither") for(i in 1:n_strata) r[s == i] <- resample(x[s == i])
      tf <- function(x) x
      if(priority == "rows") tf <- t

      # randomize community
      for(i in 1:n_strata){
            # bb <- tf(cpr_rand_comm(b[i,,], algorithm, n_iterations = n_iter, ...))
            null <- vegan::nullmodel(b[i,,], method = algorithm)
            bb <- tf(stats::simulate(null, nsim = 1, thin = 1, burnin = n_iter - 1, seed = NULL))
            bb[bb == 1] <- tf(r)[tf(s) == i]
            b[i,,] <- tf(bb)
      }
      b <- apply(b, 2:3, sum) # sum across strata
      b
}


#' Null model randomization analysis of alpha diversity metrics
#'
#' This function compares to diversity metrics calculated in \link{sphy_diversity} to their null distributions
#' computed by randomizing the community matrix. Randomization is done using the \link{strand} method for
#' community matrices containing continuous occurrence quantities such as occurrence probabilities or abundances.
#'
#' @param sp `spatialphy` object.
#' @param n_rand Integer giving the number of random communities to generate.
#' @param spatial Boolean: should the function return a spatial object (TRUE, default) or a matix (FALSE).
#' @param n_cores Integer giving the number of compute cores to use for parallel processing.
#' @param ... Additional arguments passed to \link{strand}, such as \code{n_strata}, \code{jitter}, \code{transform},
#'    \code{priority}, etc.
#' @return A matrix with a row for every row of \code{x}, a column for every metric in \link{sphy_diversity}, and
#'    values indicating the proportion of randomizations in which the observed diversity metric was greater than
#'    the randomized metric.
#' @export
sphy_rand <- function(sp, n_rand = 100, spatial = T, n_cores = 1, ...){
      phy <- sp$tree
      tip_occs <- get_tip_occs(sp)

      div <- sphy_diversity(sp, spatial = F)
      rand <- array(NA, c(dim(div), n_rand + 1))
      rand[,,1] <- div

      perm <- function(comm, tree, ...){
            rcomm <- strand(comm, ...)
            rsp <- sphy(tree, rcomm)
            sphy_diversity(rsp)
      }

      if(n_cores == 1){
            pb <- txtProgressBar(min = 0, max = n_rand, initial = 0, style = 3)
            for(i in 1:n_rand){
                  rand[,,i+1] <- perm(tip_occs, phy, ...)
                  setTxtProgressBar(pb, i)
            }
            close(pb)
      }else{
            require(furrr)
            plan(multisession, workers = n_cores)
            rnd <- future_map(1:n_rand,
                              function(i) perm(tip_occs, phy, ...),
                              .progress = TRUE,
                              .options = furrr_options(seed = TRUE))
            for(i in 1:n_rand) rand[,,i+1] <- rnd[[i]]
      }

      q <- apply(rand, 1:2, function(x) mean(x[1] > x[2:(n_rand+1)], na.rm = T) )
      colnames(q) <- paste0("q", colnames(div))

      if(spatial & !is.null(sp$spatial)) q <- to_raster(q, sp$spatial)
      return(q)
}



#' Binary randomization tests including CANAPE
#'
#' This function is a wrapper around the \code{cpr_rand_test} function from Joel Nitta's \code{canaper} package.
#' In keeping with CANAPE, it requires a binary community matrix and uses a hard significance threshold (for an
#' alternative that utilizes continuous occurrences values and retains a significance gradient, see
#' \link{snape} and \link{sphy_rand}).
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
sphy_canape <- function(sp, null_model = "curveball", spatial = T, ...){

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
