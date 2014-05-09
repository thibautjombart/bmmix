\name{Bayesian multinomial mixture model}
\alias{bmmix}
\title{Bayesian multinomial mixture model estimation using MCMC}
\description{
  This function and model are under development. Do not use them,
  contact the author if interested.
}
\usage{
bmmix(x, y, n=5e4, sample.every=200,
      move.alpha=TRUE, move.phi=FALSE,
      sd.alpha=0.1, sd.phi=0.05, move.phi.every=10,
      model.unsampled=FALSE, prior.unsampled.contrib=0.1,
      file.out="mcmc.txt", quiet=FALSE)
}
\arguments{
  \item{x}{a matrix containing multinomial data in columns used to
    compose the mixture (i.e., each column is an 'origin').}
  \item{y}{a vector of the same length as the number of rows of \code{x}
    containing the response variable.}
  \item{sample.every}{an integer indicating the frequency at which to
    save MCMC samples.}
  \item{move.alpha}{a logical indicating whether the mixture
    coefficients (alpha) should be estimated.}
  \item{move.phi}{a logical indicating whether the frequencies in \code{x}
    (phi) should be estimated; see details.}
  \item{sd.alpha}{the standard deviation of the normal distribution used
    as proposal for alpha.}
  \item{sd.phi}{the standard deviation of the normal distribution used
    as proposal for phi.}
  \item{move.phi.every}{the frequency at which values of phi should be moved.}
  \item{model.unsampeld}{a logical indicating whether an 'unsampled
    origin' should be allowed; if TRUE, then \code{move.phi} should be
    TRUE as well, to allow for frequencies in this group to be estimated.}
  \item{prior.unsampled.contrib}{the mean of the exponential
    distribution used as prior for the contribution of the unsampled
    origin in the mixture; all other origins have flat priors.}
  \item{file.out}{the name of the file used to store the outputs.}
  \item{quiet}{a logical indicating whether output messages should be hidden.}
}
\value{
  A data.frame containing the MCMC outputs.
}
\author{
  Thibaut Jombart \email{t.jombart@imperial.ac.uk}
}
\examples{
## GENERATE TOY DATA ##
set.seed(1)
f.dogs <- c(.5, .3, .1, .1, 0)
f.cows <- c(.6, .1, .1, .1, .1)
f.asymp <- c(0, .1, .2, 0, .7)
f.y <- .1*f.dogs + .1*f.cows + .8*f.asymp

dogs <- rmultinom(1, 30, f.dogs)
cows <- rmultinom(1, 50, f.cows)
asymp <- rmultinom(1, 35, f.asymp)
X <- data.frame(dogs, cows, asymp, row.names=paste("ST", letters[1:5]))
y <- rmultinom(1, 40, f.y)
cbind(X,y)


## RUN BMMIX ##
set.seed(1)
res <- bmmix(X,y)


## VISUALIZE RESULTS ##
if(require("ggplot2") && require("reshape2")){
ggplot(dat=res) + geom_line(aes(x=step, y=post))

fig.dat <- melt(res, id=1:4)
ggplot(dat=fig.dat, aes(x=step)) + geom_line(aes(y=value, colour=variable))
ggplot(dat=fig.dat, aes(x=step)) + geom_density(aes(x=value, fill=variable),alpha=.2)
}

}
