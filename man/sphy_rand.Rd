% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rand.R
\name{sphy_rand}
\alias{sphy_rand}
\title{Null model randomization analysis of alpha diversity metrics}
\usage{
sphy_rand(sp, n_rand = 100, spatial = T, n_cores = 1, ...)
}
\arguments{
\item{sp}{\code{spatialphy} object.}

\item{n_rand}{Integer giving the number of random communities to generate.}

\item{spatial}{Boolean: should the function return a spatial object (TRUE, default) or a matix (FALSE).}

\item{n_cores}{Integer giving the number of compute cores to use for parallel processing.}

\item{...}{Additional arguments passed to \link{strand}, such as \code{n_strata}, \code{jitter}, \code{transform},
\code{priority}, etc.}
}
\value{
A matrix with a row for every row of \code{x}, a column for every metric in \link{sphy_diversity}, and
values indicating the proportion of randomizations in which the observed diversity metric was greater than
the randomized metric.
}
\description{
This function compares to diversity metrics calculated in \link{sphy_diversity} to their null distributions
computed by randomizing the community matrix. Randomization is done using the \link{strand} method for
community matrices containing continuous occurrence quantities such as occurrence probabilities or abundances.
}
