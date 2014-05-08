\name{Bayesian multinomial mixture model}
\alias{bmmix}
\title{Bayesian multinomial mixture model estimation using MCMC}
\description{
  This function and model are under development. Do not use them,
  contact the author if interested.
}
\usage{
bmmix()
}
\arguments{
  \item{x}{}
}
\value{
  A data.frame containing main outputs of the MCMC.
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
res <- bmmix(X,y, n=4e4, sample.every=1000, move.sd=.3)


## VISUALIZE RESULTS ##
if(require("ggplot2") && require("reshape2")){
ggplot(dat=res) + geom_line(aes(x=step, y=post))

fig.dat <- melt(res, id=1:4)
ggplot(dat=fig.dat, aes(x=step)) + geom_line(aes(y=value, colour=variable))
ggplot(dat=fig.dat, aes(x=step)) + geom_density(aes(x=value, fill=variable),alpha=.2)
}
}
