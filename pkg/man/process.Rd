\name{Process bmmix outputs}
\alias{process.bmmix}
\title{Extract results and make graphics from bmmix outputs}
\description{
  This function and model are under development. Do not use them,
  contact the author if interested.
}
\usage{
process.bmmix(x, what="post", burnin=1e4, ggplot=TRUE)
}
\arguments{
  \item{x}{a data.frame output by \code{bmmix}.}

  \item{what}{a character string indicating which result is seeked,
    matched against the names of the columns of \code{x}.}

  \item{burnin}{an integer indicating the burnin, i.e. the number of
    MCMC iterations to be discarded.}

  \item{ggplot}{a logical indicating whether graphics using
    \code{ggplot2} should be returned.}
}
\value{
  A list containing processed results, and optionally ggplot graphics.
}
\author{
  Thibaut Jombart \email{t.jombart@imperial.ac.uk}
}

\examples{
## GENERATE TOY DATA ##
## ST frequencies in 3 origins:
## dogs, cows, asymptotic cases in human
f.dogs <- c(.5, .3, .1, .1, 0)
f.cows <- c(.6, .1, .1, .1, .1)
f.asymp <- c(0, .1, .2, 0, .7)

## mixture (y would be symptomatic cases)
f.y <- .1*f.dogs + .1*f.cows + .8*f.asymp

set.seed(1)
dogs <- rmultinom(1, 30, f.dogs)
cows <- rmultinom(1, 50, f.cows)
asymp <- rmultinom(1, 35, f.asymp)
X <- data.frame(dogs, cows, asymp,
   row.names=paste("ST", letters[1:5]))
X
y <- rmultinom(1, 40, f.y)
y

cbind(X,y)


## RUN BMMIX ##

## BASIC MODEL
## note: small n for this example only!
set.seed(1)
res <- bmmix(X,y, n=3e4)
head(res)


## VISUALIZE RESULTS ##
if(require("ggplot2") && require("reshape2")){
## trace of log posterior
process.bmmix(res)

## mixture coefficients
temp <- process.bmmix(res, "alpha")
names(temp)
temp$alpha # values 
temp$trace # graphics: trace
temp$hist # graphics: histograms
temp$dens # graphics: densities
temp$violin # graphics: violinplot

}
