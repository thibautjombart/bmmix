#' Bayesian multinomial mixture model estimation using MCMC
#'
#' This function and model are under development. Do not use them, contact the
#' author if interested.
#'
#' There are essentially 4 variants of the model implemented by \code{bmmix}:
#' \itemize{  \item estimate mixture only (default)only the mixture coefficients
#' are estimated; frequencies (phi) are fixed to their maximum likelihood
#' estimate from the data; this model has 'K' parameters (where 'K' is the
#' number of putative origins, i.e. the number of columns in 'x').
#'
#'  \item estimate mixture and frequenciesboth mixture coefficients and
#' frequencies for each group and origin are estimated; this model has (N+1)K
#' parameters (N being the number of rows in 'x'); to use this model, specify
#' \code{move.phi=TRUE}.
#'
#'  \item estimate mixture, allowing unsampled originmixture coefficients are
#' estimated with an additional 'unsampled' origin whose frequencies are
#' estimated; this model has K+N+1 + parameters (N being the number of rows in
#' 'x'); to use this model, specify \code{move.phi=FALSE,model.unsampled=TRUE};
#' this is the only practical model allowing unsampled origin for medium-sized
#' or large datasets..
#'
#'  \item estimate mixture, frequencies, and allow unsampled originthis is the
#' most complex model; in addition to the previous one, an unsampled origin is
#' allowed, and its frequencies are estimated; this model therefore has
#' (N+1)(K+1) parameters; to use this model, specify \code{move.phi=TRUE} and
#' \code{model.unsampled=TRUE}; note that if frequencies are not estimated
#' (\code{move.phi=FALSE}), the frequencies in the unsampled origin will be
#' fixed to their initial value in which all groups have the same frequency;
#' this model quickly becomes hard to fit for medium-sized to large datasets.
#'
#' }
#'
#' @aliases bmmix
#'
#' @param x a matrix containing multinomial data in columns used to compose the
#' mixture (i.e., each column is an 'origin').
#' @param y a vector of the same length as the number of rows of \code{x}
#' containing the response variable.
#' @param n the length of the MCMC.
#' @param sample.every an integer indicating the frequency at which to save
#' MCMC samples.
#' @param move.alpha a logical indicating whether the mixture coefficients
#' (alpha) should be estimated.
#' @param move.phi a logical indicating whether the frequencies in \code{x}
#' (phi) should be estimated; see details.
#' @param sd.alpha the standard deviation of the normal distribution used as
#' proposal for alpha.
#' @param sd.phi the standard deviation of the normal distribution used as
#' proposal for phi.
#' @param move.phi.every the frequency at which values of phi should be moved.
#' @param model.unsampled a logical indicating whether an 'unsampled origin'
#' should be allowed; if TRUE, then \code{move.phi} should be TRUE as well, to
#' allow for frequencies in this group to be estimated.
#' @param prior.unsampled.contrib the mean of the exponential distribution used
#' as prior for the contribution of the unsampled origin in the mixture; all
#' other origins have flat priors.
#' @param min.ini.freq the default minimum frequency of unobserved items in
#' \code{x} used for the initial frequency estimate.
#' @param file.out the name of the file used to store the outputs.
#' @param quiet a logical indicating whether output messages should be hidden.

#' @return A data.frame with class 'bmmix', containing the MCMC outputs: step,
#' log-posterior, log- likelihood, log-prior, alpha values (mixture
#' coefficients), and optionally frequencies for each group and origin (phi).
#'
#' @importFrom gtools rdirichlet ddirichlet
#'
#' @export bmmix
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @examples
#'
#' ## GENERATE TOY DATA ##
#' ## ST frequencies in 3 origins:
#' ## dogs, cows, asymptotic cases in human
#' f.dogs <- c(.5, .3, .1, .1, 0)
#' f.cows <- c(.6, .1, .1, .1, .1)
#' f.asymp <- c(0, .1, .2, 0, .7)
#'
#' ## mixture (y would be symptomatic cases)
#' f.y <- .1*f.dogs + .1*f.cows + .8*f.asymp
#'
#' set.seed(1)
#' dogs <- rmultinom(1, 30, f.dogs)
#' cows <- rmultinom(1, 50, f.cows)
#' asymp <- rmultinom(1, 35, f.asymp)
#' X <- data.frame(dogs, cows, asymp,
#'    row.names=paste("ST", letters[1:5]))
#' X
#' y <- rmultinom(1, 40, f.y)
#' y
#'
#' cbind(X,y)
#'
#'
#' ## RUN BMMIX ##
#'
#' ## BASIC MODEL
#' ## note: small n for this example only!
#' set.seed(1)
#' res <- bmmix(X,y, n=3e4)
#' head(res)
#'
#'
#' ## VISUALIZE RESULTS ##
#' if(require("ggplot2") && require("reshape2")){
#'
#' ## manually ##
#' ## chech log-posterior
#' ggplot(dat=res) + geom_line(aes(x=step, y=post)) +
#'    labs(title="Trace of log-posterior values")
#'
#' ## check mixture coefficients
#' fig.dat <- melt(res, id=1:4)
#'
#' ggplot(dat=fig.dat, aes(x=step)) +
#'    geom_line(aes(y=value, colour=variable)) +
#'    labs(title="Trace of mixture coefficients")
#'
#'
#' ## with process.bmmix ##
#' ## mixture coefficients
#' temp <- process.bmmix(res, "alpha")
#' names(temp)
#' temp$alpha # values
#' temp$trace # graphics: trace
#' temp$hist # graphics: histograms
#' temp$dens # graphics: densities
#' temp$violin # graphics: violinplot
#'
#' }
#'
#'
#' \dontrun{
#' ## MODEL WITH ESTIMATED FREQUENCIES
#' set.seed(1)
#' res <- bmmix(X,y, move.phi=TRUE)
#' head(res)
#'
#' ## VISUALIZE RESULTS
#' if(require("ggplot2") && require("reshape2")){
#'
#' ## chech log-posterior
#' ggplot(dat=res) + geom_line(aes(x=step, y=post)) +
#'    labs(title="Trace of log-posterior values")
#'
#' fig.dat <- melt(res[,1:7], id=1:4)
#'
#' ## check mixture coefficients
#' ggplot(dat=fig.dat, aes(x=step)) +
#'    geom_line(aes(y=value, colour=variable)) +
#'    labs(title="Trace of mixture coefficients")
#'
#' ## check ST frequencies, i.e. in dogs:
#' fig.dat <- melt(res[,c(1,grep("dogs", names(res))[-1])], id=1)
#'
#' ggplot(dat=fig.dat) +
#'    geom_line(aes(x=step, y=value, colour=variable)) +
#'    labs(title="Estimate of ST frequencies in dogs")
#'
#' ggplot(dat=fig.dat) +
#'    geom_density(aes(x=value, fill=variable),alpha=.2) +
#'    labs(title="Estimate of ST frequencies in dogs")
#' }
#' }


##########
## bmmix
##########
bmmix <- function(x, y, n=5e4, sample.every=200,
                  move.alpha=TRUE, move.phi=FALSE,
                  sd.alpha=0.1, sd.phi=0.05, move.phi.every=10,
                  model.unsampled=FALSE, prior.unsampled.contrib=0.1,
                  min.ini.freq=0.01,
                  file.out="mcmc.txt", quiet=FALSE){
  ## CHECKS ##
  NEARZERO <- 1e-20
  if(n/sample.every < 10) warning("less than 10 samples are going to be produced")
  x <- as.matrix(x)
  if(model.unsampled){
    x <- cbind(x, unsampled=rep(0, nrow(x)))
    rate.alpha.prior <- 1/prior.unsampled.contrib
  }
  K <- ncol(x)
  N <- nrow(x)
  if(N != length(y)) stop("The number of rows in x does not match the length of y")


  ## MAKE SURE THERE'S NO EMPTY SOURCE
  if (any(!is.finite(x))) {
    stop("Non-finite values detected in 'x'")
  }
  if (any(!is.finite(y))) {
    stop("Non-finite values detected in 'y'")
  }
  sources_sums <- colSums(x)
  if (any(sources_sums < 1)) {
    stop("Some of the sources have no observation")
  }

  
  ## LIKELIHOOD FUNCTIONS ##
  ## LIKELIHOOD OF CASE DATA 'Y' ##
  ## y is a vector of N numbers
  ## phi is a NxK matrix of numbers
  ## alpha is a vector of K mixture coefficients
  LL.y <- function(y, phi, alpha){
    phi.y <- phi %*% (alpha/sum(alpha))
    return(dmultinom(y, prob=phi.y, log=TRUE))
  }


  ## LIKELIHOOD OF PUTATIVE ORIGIN DATA 'x'
  ## y is a vector of N numbers
  ## phi is a NxK matrix of numbers
  ## alpha is a vector of K mixture coefficients
  LL.x <- function(x, phi){
    return(sum(sapply(1:ncol(x), function(i) dmultinom(x[,i], prob=phi[,i], log=TRUE))))
  }


  ## LIKELIHOOD OF ALL DATA ##
  LL.all <- function(y, x, phi, alpha){
    return(LL.y(y, phi, alpha) + LL.x(x, phi))
  }


  ## PRIOR FUNCTIONS ##
  if(model.unsampled){
    LPrior.alpha <- function(alpha){
      return(dexp((alpha/sum(alpha))[K], rate=rate.alpha.prior, log=TRUE))
    }
  } else {
    LPrior.alpha <- function(alpha){
      return(0)
    }
  }




  ## POSTERIOR FUNCTIONS ##
  LPost.all <- function(y, x, phi, alpha){
    return(LL.all(y,x, phi, alpha) + LPrior.alpha(alpha))
  }




  ## MOVEMENT FUNCTIONS ##
  ## MOVE ALPHA
  ALPHA.ACC <- 0
  ALPHA.REJ <- 0
  SCALE.ALPHA <- 10
  alpha.move <- function(alpha, sigma=sd.alpha){
    ## propose new vector  ##
    newval <- as.vector(rdirichlet(1, alpha*SCALE.ALPHA))

    ## correction factor for MH ##
    ## LL.corr = log(p(new->old)) - log(p(old->new))
    LL.correc <- log(ddirichlet(alpha, newval*SCALE.ALPHA)) - log(ddirichlet(newval, alpha*SCALE.ALPHA))

    if(all(newval>=0 & newval<=1)){
      metro.ratio <- LL.y(y, phi, newval) - LL.y(y, phi, alpha) +
        LPrior.alpha(newval) - LPrior.alpha(alpha) +
        LL.correc
      if((r <- log(runif(1))) <=  metro.ratio){
        alpha <- newval # accept
        ALPHA.ACC <<- ALPHA.ACC+1
      } else {
        ALPHA.REJ <<- ALPHA.REJ+1
      }
    } else {
      ALPHA.REJ <<- ALPHA.REJ+1
    }

    ## return moved vector
    return(alpha)
  }


  ## MOVE PHI
  PHI.ACC <- 0
  PHI.REJ <- 0
  SCALE.PHI <- 10
  phi.move <- function(phi, sigma=sd.phi){
    ## check which one must move
    if(move.phi) {
      idx.toMove <- 1:K
    } else if(model.unsampled){
      idx.toMove <- K
    } else return(phi)

    ## for all frequencies to move...
    for(tomove in idx.toMove){

      ## generate all proposals ##
      ## newval <- rnorm(n=nrow(phi), mean=phi[,tomove], sd=sigma)
      ## newval <- newval/sum(newval)
      newval <- as.vector(rdirichlet(1, phi[,tomove]*SCALE.PHI))
      newphi <- phi
      newphi[,tomove] <- newval

      ## correction factor for MH ##
      ## LL.corr = log(p(new->old)) - log(p(old->new))
      LL.correc <- log(ddirichlet(phi[,tomove], newval*SCALE.PHI)) - log(ddirichlet(newval, phi[,tomove]*SCALE.PHI))

      if(all(newval>=0 & newval<=1)){
        if((r <- log(runif(1))) <=  (LL.all(y, x, newphi, alpha) - LL.all(y, x, phi, alpha))){
          phi <- newphi # accept
          PHI.ACC <<- PHI.ACC+1
        } else {
          PHI.REJ <<- PHI.REJ+1
        }
      } else {
        PHI.REJ <<- PHI.REJ+1
      }
    }

    ## return moved vector
    return(phi)
  }




  ## MAIN MCMC FUNCTION ##
  ## INITIALIZE MCMC
  ## initial alpha
  alpha <- rep(1/K, K)

  ## initial phi
  if(model.unsampled){
    phi <- cbind(prop.table(x[,-K],2), "unsampled"=rep(1/nrow(x), nrow(x)))
  } else {
    phi <- prop.table(x,2)
  }

  ## handle 'zero' replacement
  nb.toreplace <- apply(phi,2, function(e) sum(e<NEARZERO))
  replace.freq <- min.ini.freq/nb.toreplace
  freq.tosubstract <- 0.01/(nrow(x)-nb.toreplace)

  for(j in 1:ncol(phi)){
    ## replace zeros with small freq
    phi.arezero <- phi[,j] < NEARZERO
    phi[phi.arezero,j] <- replace.freq[j]
    phi[!phi.arezero,j] <- phi[!phi.arezero,j] - freq.tosubstract[j]
  }


  ## ADD HEADER TO THE OUTPUT FILE
  ## basic header
  header <- "step\tpost\tlikelihood\tprior"

  ## header for alpha
  if(move.alpha) header <- c(header, paste("alpha", colnames(x), sep=".", collapse="\t"))

  ## header for phi
  if(move.phi) {
    annot.phi <- paste("phi",
                       rownames(phi)[as.vector(row(phi))],
                       colnames(phi)[as.vector(col(phi))],
                       sep=".", collapse="\t")
    header <- c(header, annot.phi)
  } else if(model.unsampled){
    annot.phi <- paste("phi",
                       rownames(phi),
                       "unsampled",
                       sep=".", collapse="\t")
    header <- c(header, annot.phi)
  }

  ## collapse everything and write to file
  header <- paste(header, collapse="\t")
  cat(header, file=file.out)


  ## add first line
  ## temp: c(loglike, logprior)
  temp <- c(LL.all(y, x, phi, alpha), LPrior.alpha(alpha))

  ## check that initial LL is not -Inf
  if(!is.finite(temp[1])) warning("Initial likelihood is zero")

  ## write to file
  cat("\n", file=file.out, append=TRUE)
  cat(c(1, sum(temp), temp), sep="\t", append=TRUE, file=file.out)
  if(move.alpha) cat("", alpha/sum(alpha), sep="\t", append=TRUE, file=file.out)
  if(move.phi){
    cat("", as.vector(phi), sep="\t", append=TRUE, file=file.out)
  } else if(model.unsampled){
    cat("", as.vector(phi[,K]), sep="\t", append=TRUE, file=file.out)
  }

  if(!quiet) cat("\nStarting MCMC: 1")

  ## mcmc ##
  for(i in 1:n){
    ## move stuff ##
    ## move alpha if needed
    if(move.alpha) alpha <- alpha.move(alpha, sd.alpha)

    ## move phi if needed (phi.move makes the necessary moves)
    if((move.phi|model.unsampled) && (i %% move.phi.every == 0)) phi <- phi.move(phi, sd.phi)

    ## if retain this sample ##
    if(i %% sample.every ==0){
      temp <- c(LL.all(y, x, phi, alpha), LPrior.alpha(alpha))
      cat("\n", file=file.out, append=TRUE)
      cat(c(i, sum(temp), temp), sep="\t", append=TRUE, file=file.out)
      if(move.alpha) cat("", alpha/sum(alpha), sep="\t", append=TRUE, file=file.out)
      if(move.phi) {
        cat("", as.vector(phi), sep="\t", append=TRUE, file=file.out)
      } else if(model.unsampled){
        cat("", as.vector(phi[,K]), sep="\t", append=TRUE, file=file.out)
      }
      if(!quiet) cat("..",i)
    }
  }

  if(!quiet) cat("..done!\nResults were saved in file:",file.out,"\n")


  ## re-read output file ##
  out <- read.table(file.out, header=TRUE, colClasses="numeric", sep="\t")
  out$step <- as.integer(out$step)

  ## format using coda ##
  if(!quiet){
    ## acceptance rates for alpha
    if(move.alpha){
      cat("\nacceptance rate for alpha: ", ALPHA.ACC/(ALPHA.ACC+ALPHA.REJ))
      cat("\naccepted: ", ALPHA.ACC)
      cat("\nreject: ", ALPHA.REJ)
      cat("\n")
    }

    ## acceptance rates for phi
    if(move.phi || model.unsampled){
      cat("\nacceptance rate for phi: ", PHI.ACC/(PHI.ACC+PHI.REJ))
      cat("\naccepted: ", PHI.ACC)
      cat("\nreject: ", PHI.REJ)
      cat("\n")
    }

  }


  class(out) <- c("data.frame", "bmmix")
  return(out)

} # end bmmix
