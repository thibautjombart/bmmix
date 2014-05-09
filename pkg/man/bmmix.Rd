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
  
  \item{n}{the length of the MCMC.}
  
  \item{sample.every}{an integer indicating the frequency at which to
    save MCMC samples.}
  
  \item{move.alpha}{a logical indicating whether the mixture
    coefficients (alpha) should be estimated.}
  
  \item{move.phi}{a logical indicating whether the frequencies in
    \code{x} (phi) should be estimated; see details.}
  
  \item{sd.alpha}{the standard deviation of the normal distribution used
    as proposal for alpha.}
  
  \item{sd.phi}{the standard deviation of the normal distribution used
    as proposal for phi.}
  
  \item{move.phi.every}{the frequency at which values of phi should be
    moved.}
  
  \item{model.unsampled}{a logical indicating whether an 'unsampled
    origin' should be allowed; if TRUE, then \code{move.phi} should be
    TRUE as well, to allow for frequencies in this group to be
    estimated.}
  
  \item{prior.unsampled.contrib}{the mean of the exponential
    distribution used as prior for the contribution of the unsampled
    origin in the mixture; all other origins have flat priors.}
  
  \item{file.out}{the name of the file used to store the outputs.}
  
  \item{quiet}{a logical indicating whether output messages should be
    hidden.}
}
\value{
  A data.frame containing the MCMC outputs: step, log-posterior, log-
  likelihood, log-prior, alpha values (mixture coefficients), and
  optionally frequencies for each group and origin (phi).
}
\author{
  Thibaut Jombart \email{t.jombart@imperial.ac.uk}
}
\details{
  There are essentially 3 variants of the model implemented by
  \code{bmmix}:
  \itemize{
    \item{estimate mixture only (default)}{only the mixture coefficients
      are estimated; frequencies (phi) are fixed to their maximum
      likelihood estimate from the data; this model has 'K' parameters
      (where 'K' is the number of putative origins, i.e. the number of
      columns in 'x').}
    
    \item{estimate mixture and frequencies}{both mixture coefficients
      and frequencies for each group and origin are estimated; this
      model has (N+1)K parameters (N being the number of rows in 'x');
      to use this model, specify \code{move.phi=TRUE}.}
    
    \item{estimate mixture, frequencies, and allow unsampled
      origin}{this is the most complex model; in addition to the
      previous one, an unsampled origin is allowed, and its frequencies
      are estimated; this model therefore has (N+1)(K+1) parameters; to
      use this model, specify \code{move.phi=TRUE} and
      \code{model.unsampled=TRUE}; note that if frequencies are not
      estimated (\code{move.phi=FALSE}), the frequencies in the
      unsampled origin will be fixed to their initial value in which all
      groups have the same frequency.}
  }
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

## chech log-posterior
ggplot(dat=res) + geom_line(aes(x=step, y=post)) +
   labs(title="Trace of log-posterior values")

## check mixture coefficients
fig.dat <- melt(res, id=1:4)

ggplot(dat=fig.dat, aes(x=step)) +
   geom_line(aes(y=value, colour=variable)) +
   labs(title="Trace of mixture coefficients")

ggplot(dat=fig.dat, aes(x=step)) +
   geom_density(aes(x=value, fill=variable),alpha=.2) +
   labs(title="Densities of mixture coefficients")
}


\dontrun{
## MODEL WITH ESTIMATED FREQUENCIES
set.seed(1)
res <- bmmix(X,y, move.phi=TRUE)
head(res)

## VISUALIZE RESULTS
if(require("ggplot2") && require("reshape2")){

## chech log-posterior
ggplot(dat=res) + geom_line(aes(x=step, y=post)) +
   labs(title="Trace of log-posterior values")

fig.dat <- melt(res[,1:7], id=1:4)

## check mixture coefficients
ggplot(dat=fig.dat, aes(x=step)) +
   geom_line(aes(y=value, colour=variable)) +
   labs(title="Trace of mixture coefficients")

## check ST frequencies, i.e. in dogs:
fig.dat <- melt(res[,c(1,grep("dogs", names(res))[-1])], id=1)

ggplot(dat=fig.dat) +
   geom_line(aes(x=step, y=value, colour=variable)) +
   labs(title="Estimate of ST frequencies in dogs")

ggplot(dat=fig.dat) +
   geom_density(aes(x=value, fill=variable),alpha=.2) +
   labs(title="Estimate of ST frequencies in dogs")
}
}

}
