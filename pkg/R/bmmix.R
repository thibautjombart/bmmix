
##########
## bmmix
##########
bmmix <- function(x, y, n=5e4, sample.every=200,
                  move.alpha=TRUE, move.phi=FALSE,
                  sd.alpha=0.1, sd.phi=0.05, move.phi.every=10,
                  model.unsampled=FALSE, prior.unsampled.contrib=0.1,
                  file.out="mcmc.txt", quiet=FALSE){
    ## CHECKS ##
    if(n/sample.every < 10) warning("less than 10 samples are going to be produced")
    x <- as.matrix(x)
    if(model.unsampled){
        x <- cbind(x, unsampled=rep(0, nrow(x)))
        rate.alpha.prior <- 1/prior.unsampled.contrib
    }
    K <- ncol(x)
    N <- nrow(x)
    if(N != length(y)) stop("The number of rows in x does not match the length of y")
    if(model.unsampled && !move.phi) warning("It is strongly recommended to move phi when allowing for unsampled origins; \notherwise all ST frequencies will be equal in the unsampled origin")



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
    alpha.move <- function(alpha, sigma=sd.alpha){
        ## generate all proposals ##
        newval <- rnorm(n=length(alpha), mean=alpha, sd=sigma)
        newval <- newval/sum(newval)

        if(all(newval>=0 & newval<=1)){
            metro.ratio <- LL.y(y, phi, newval) - LL.y(y, phi, alpha) + LPrior.alpha(newval) - LPrior.alpha(alpha)
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
    phi.move <- function(phi, sigma=sd.phi){
        ## ## select phi moved ##
        ## tomove <- sample(1:ncol(phi), 1)
        for(tomove in 1:ncol(phi)){

            ## generate all proposals ##
            newval <- rnorm(n=nrow(phi), mean=phi[,tomove], sd=sigma)
            newval <- newval/sum(newval)
            temp <- phi
            temp[,tomove] <- newval

            if(all(newval>=0 & newval<=1)){
                if((r <- log(runif(1))) <=  (LL.all(y, x, temp, alpha) - LL.all(y, x, phi, alpha))){
                    phi <- temp # accept
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
    alpha <- rep(1,K)

    ## initial phi
    phi <- prop.table(x,2)
    if(model.unsampled){
        phi[,K] <- rep(1/nrow(x), nrow(x))
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
    }

    ## collapse everything and write to file
    header <- paste(header, collapse="\t")
    cat(header, file=file.out)


    ## add first line
    ## temp: c(loglike, logprior)
    temp <- c(LL.all(y, x, phi, alpha), LPrior.alpha(alpha))
    cat("\n", file=file.out, append=TRUE)
    cat(c(1, sum(temp), temp), sep="\t", append=TRUE, file=file.out)
    if(move.alpha) cat("", alpha/sum(alpha), sep="\t", append=TRUE, file=file.out)
    if(move.phi) cat("", as.vector(phi), sep="\t", append=TRUE, file=file.out)

    if(!quiet) cat("\nStarting MCMC: 1")

    ## mcmc ##
    for(i in 1:n){
        ## move stuff ##
        ## move alpha if needed
        if(move.alpha) alpha <- alpha.move(alpha, sd.alpha)

        ## move phi if needed
        if(move.phi && (i %% move.phi.every == 0)) phi <- phi.move(phi, sd.phi)

        ## if retain this sample ##
        if(i %% sample.every ==0){
            temp <- c(LL.all(y, x, phi, alpha), LPrior.alpha(alpha))
            cat("\n", file=file.out, append=TRUE)
            cat(c(i, sum(temp), temp), sep="\t", append=TRUE, file=file.out)
            if(move.alpha) cat("", alpha/sum(alpha), sep="\t", append=TRUE, file=file.out)
            if(move.phi) cat("", as.vector(phi), sep="\t", append=TRUE, file=file.out)
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
        if(move.phi){
            cat("\nacceptance rate for phi: ", PHI.ACC/(PHI.ACC+PHI.REJ))
            cat("\naccepted: ", PHI.ACC)
            cat("\nreject: ", PHI.REJ)
            cat("\n")
        }

    }

    return(out)

} # end bmmix
