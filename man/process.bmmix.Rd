% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/processing.R
\name{process.bmmix}
\alias{process.bmmix}
\title{Extract results and make graphics from bmmix outputs}
\usage{
process.bmmix(x, what = "post", burnin = 10000, ggplot = TRUE)
}
\arguments{
\item{x}{a data.frame output by \code{bmmix}.}

\item{what}{a character string indicating which result is seeked, matched
against the names of the columns of \code{x}.}

\item{burnin}{an integer indicating the burnin, i.e. the number of MCMC
iterations to be discarded.}

\item{ggplot}{a logical indicating whether graphics using \code{ggplot2}
should be returned.}
}
\value{
A list containing processed results, and optionally ggplot graphics.
}
\description{
This function and model are under development. Do not use them, contact the
author if interested.
}
\author{
Thibaut Jombart \email{t.jombart@imperial.ac.uk}
}

