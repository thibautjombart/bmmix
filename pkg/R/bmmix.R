
##########
## bmmix
##########
bmmix <- function(X, y, n=1e5, sample.every=500,
                          phi.move=FALSE, alpha.prior=NULL,
                          move.sd=0.5,
                          file.out="mcmc.txt", quiet=FALSE){
    ## CHECKS ##
    if(n/sample.every < 10) warning("less than 10 samples are going to be produced")
    X <- as.matrix(X)
    K <- ncol(X)
    N <- nrow(X)
    if(N != length(y)) stop("The number of rows in X does not match the length of y")
    if(phi.move) stop("phi cannot be moved in this version")
    ## if(is.null(alpha.prior)){
    ##     alpha.prior <- rep(1/ncol(X), ncol(X))
    ## }


    ## LIKELIHOOD FUNCTIONS ##
    ## LIKELIHOOD OF CASE DATA 'Y' ##
    ## y is a vector of N numbers
    ## phi is a NxK matrix of numbers
    ## alpha is a vector of K mixture coefficients
    LL.y <- function(y, phi, alpha){
        phi.y <- phi %*% (alpha/sum(alpha))
        return(dmultinom(y, prob=phi.y, log=TRUE))
    }


    ## LIKELIHOOD OF PUTATIVE ORIGIN DATA 'X'
    ## y is a vector of N numbers
    ## phi is a NxK matrix of numbers
    ## alpha is a vector of K mixture coefficients
    LL.X <- function(X, phi){
        return(sum(sapply(1:ncol(X), function(i) dmultinom(X[,i], prob=phi[,i], log=TRUE))))
    }


    ## LIKELIHOOD OF ALL DATA ##
    LL.all <- function(y, X, phi, alpha){
        return(LL.y(y, phi, alpha) + LL.X(X, phi))
    }


    ## PRIOR FUNCTIONS ##
    LPrior.alpha <- function(alpha){
        return(0)
    }



    ## POSTERIOR FUNCTIONS ##
    LPost.all <- function(y, X, phi, alpha){
        return(LL.all(y,X, phi, alpha) + LPrior.alpha(alpha))
    }


    ## MOVEMENT FUNCTIONS ##
    ## MOVE ALPHA
    COUNT.ACC <- 0
    COUNT.REJ <- 0
    move.alpha <- function(alpha, sigma=move.sd){
        ## generate all proposals ##
        newval <- rnorm(n=length(alpha), mean=alpha, sd=sigma)
        newval <- newval/sum(newval)

        if(all(newval>=0 & newval<=1)){
            if((r <- log(runif(1))) <=  (LL.y(y, phi, newval) - LL.y(y, phi, alpha))){
                alpha <- newval # accept
                COUNT.ACC <<- COUNT.ACC+1
            } else {
                COUNT.REJ <<- COUNT.REJ+1
            }
        } else {
            COUNT.REJ <<- COUNT.REJ+1
        }


        ## ## for all alpha values (origins)
        ## for(i in 1:length(alpha)){
        ##     if(newval[i]>=0){
        ##         temp[i] <- newval[i]

        ##         ## Metropolis acceptance rule
        ##         ## ! assumes flat prior on alpha for now
        ##         if((r <- log(runif(1))) <=  (LL.y(y, phi, temp) - LL.y(y, phi, alpha))){
        ##             ## debugging ##
        ##             ## if(temp[i]<0) {
        ##             ##     cat("!!accepted a negative alpha!!")
        ##             ##     cat("\n\nalpha:")
        ##             ##     print(alpha)
        ##             ##     cat("\n\ntemp:")
        ##             ##     print(temp)
        ##             ##     cat("\n\nr:")
        ##             ##     print(r)
        ##             ##     cat("\n\nLL temp:")
        ##             ##     print(LL.y(y, phi, temp))
        ##             ##     cat("\n\nLL old:")
        ##             ##     print(LL.y(y, phi, alpha))
        ##             ##     return(list(LL.y,y,phi,temp,alpha))
        ##             ## }
        ##             alpha[i] <- temp[i] # accept
        ##             COUNT.ACC <<- COUNT.ACC+1
        ##         } else {
        ##             temp[i] <- alpha[i] # reject
        ##             COUNT.REJ <<- COUNT.REJ+1
        ##         }
        ##     }
        ## }

        ## return moved vector
        return(alpha)
    }


    ## MAIN MCMC FUNCTION ##
    ## INITIALIZE MCMC
    ## initial alpha
    alpha <- rep(1,K)

    ## initial phi
    phi <- prop.table(X,2)

    ## add header to the output file
    header <- "step\tpost\tlikelihood\tprior"
    header <- c(header, paste("alpha", 1:K, sep=".", collapse="\t"))
    header <- paste(header, collapse="\t")
    cat(header, file=file.out)


    ## add first line
    ## temp: c(loglike, logprior)
    temp <- c(LL.all(y, X, phi, alpha), LPrior.alpha(alpha))
    cat("\n", file=file.out, append=TRUE)
    cat(c(1, sum(temp), temp, alpha/sum(alpha)), sep="\t", append=TRUE, file=file.out)
    if(!quiet) cat("\nStarting MCMC: 1")

    ## mcmc ##
    for(i in 1:n){
        ## move stuff ##
        alpha <- move.alpha(alpha)

        ## if retain this sample ##
        if(i %% sample.every ==0){
            temp <- c(LL.all(y, X, phi, alpha), LPrior.alpha(alpha))
            cat("\n", file=file.out, append=TRUE)
            cat(c(i, sum(temp), temp, alpha/sum(alpha)), sep="\t", append=TRUE, file=file.out)
            if(!quiet) cat("..",i)
        }
    }

    if(!quiet) cat("..done!\nResults were saved in file:",file.out,"\n")


    ## re-read output file ##
    out <- read.table(file.out, header=TRUE, colClasses="numeric")
    out$step <- as.integer(out$step)

    ## format using coda ##
    if(!quiet){
        cat("\nacceptance rate: ", COUNT.ACC/(COUNT.ACC+COUNT.REJ))
        cat("\naccepted: ", COUNT.ACC)
        cat("\nreject: ", COUNT.REJ)
        cat("\n")
    }

    return(out)

} # end bmmix
