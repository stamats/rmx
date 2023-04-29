###############################################################################
## RMX estimator for Poisson mean
###############################################################################
rmx.pois <- function(x, eps.lower=0, eps.upper=NULL, eps=NULL, k = 3L, 
                     initial.est=NULL, fsCor = FALSE, na.rm = TRUE, 
                     M = 10000, parallel = FALSE, ncores = NULL, 
                     aUp = 100*max(x), cUp = 1e4, delta = 1e-9){
    if(length(x) <= 2){
        stop("A sample size of at least 3 is required!")
    }
    if(is.null(initial.est)){
        lambda <- cvm.pois(x, mu = "model")
    }else{
        stopifnot(is.numeric(initial.est))
        if(length(initial.est) != 1)
            stop("'initial.est' needs to be a numeric vector of length 1 or NULL")
        if(initial.est < 0)
            stop("initial estimate for lambda < 0, which is not valid")
        lambda <- initial.est
    }

    lcr <- .lcr.pois(lambda = lambda)
    if(!is.null(eps)){
        r <- sqrt(length(x))*eps
        if(fsCor){ 
            r.as <- r
            if(r < lcr){
                r <- fsRadius.pois(r = r, n = length(x), lambda = lambda,
                                   M = M, parallel = parallel, ncores = ncores)
            }
        }
    }else{
        sqrtn <- sqrt(length(x))
        rlo <- sqrtn*eps.lower
        rup <- sqrtn*eps.upper
        if(rlo >= lcr){
            r <- (rlo + rup)/2
            r.as <- r
        }else{
            r <- uniroot(.getInterval.pois, lower = rlo+1e-8, upper = rup, 
                         tol = .Machine$double.eps^0.25, rlo = rlo, 
                         rup = rup, lambda = lambda, aUp = aUp, 
                         cUp = cUp, delta = delta)$root
            r.as <- r
        }
        if(fsCor){
            r.as <- r
            if(r < lcr){
                r <- fsRadius.pois(r = r, n = length(x), lambda = lambda,
                                   M = M, parallel = parallel, ncores = ncores)
            }
        }
    }
    if(!is.null(eps)){
        if(eps == 0){
            rmxEst <- mean(x)
            names(rmxEst) <- "lambda"
            IF <- .getMLIF.pois(lambda = rmxEst)
            IF$radius <- 0
            Info.matrix <- matrix(c("rmx", "ML estimate"), ncol = 2, 
                                  dimnames = list(NULL, c("method", "message")))

            RMX <- list(rmxEst = rmxEst, rmxIF = IF, initial.est = NULL, 
                        Infos = Info.matrix, n = length(x))
            class(RMX) <- "rmx"
            return(RMX)
        }
    }
    if(!is.null(eps)){
        if(r >= lcr){
            IF <- .getMBIF.pois(lambda = lambda)
            IF$radius <- r
            IF$asBias <- r^2*IF$b^2
            IF$asMSE <- IF$asVar + IF$asBias
        }else{
            LM <- .getLM.pois(r0 = r, lambda = lambda, 
                              aUp = aUp, cUp = cUp, delta = delta)
            IF <- .getOptIF.pois(lambda = lambda, radius = r, 
                                 A = LM$A, a = LM$A*LM$z, b = LM$A*LM$c0)
            IF$radius <- r
        }
        rmxEst.all <- .kstep.pois(x = x, IF, k = k, aUp = aUp, cUp = cUp,
                                  delta = delta)
        rmxEst <- rmxEst.all$lambda
        names(rmxEst) <- "lambda"
        if(fsCor){
            Info.matrix <- matrix(c("rmx", 
                                    paste("fs-corrected estimate for 'eps' =", 
                                          round(eps, 3))),
                                  ncol = 2, dimnames = list(NULL, c("method", "message")))
        }else{
            Info.matrix <- matrix(c("rmx", 
                                    paste("asymptotic estimate for 'eps' =", 
                                          round(eps, 3))),
                                  ncol = 2, dimnames = list(NULL, c("method", "message")))
        }
        IF <- rmxEst.all$IF
        IF$radius <- r
    }else{
        if(r >= lcr){
            IF <- .getMBIF.pois(lambda = lambda)
            IF$radius <- r
            IF$asBias <- r^2*IF$b^2
            IF$asMSE <- IF$asVar + IF$asBias
        }else{
            LM <- .getLM.pois(r0 = r, lambda = lambda, 
                              aUp = aUp, cUp = cUp, delta = delta)
            IF <- .getOptIF.pois(lambda = lambda, radius = r, 
                                 A = LM$A, a = LM$A*LM$z, b = LM$A*LM$c0)
            IF$radius <- r
        }
        if(rlo == 0){
            b <- LM$A*LM$c0
            ineff <- (LM$A - b^2*r.as^2)/lambda
        }else{
            if(rlo >= lcr){
                ineff <- 1
            }else{
                b <- LM$A*LM$c0
                ineff <- (LM$A - b^2*(r.as^2 - rlo^2))/LM$A
            }
        }
        rmxEst.all <- .kstep.pois(x = x, IF = IF, k = k, aUp = aUp, cUp = cUp,
                                  delta = delta)
        rmxEst <- rmxEst.all$lambda
        names(rmxEst) <- "lambda"
        if(fsCor){
            Info.matrix <- matrix(c(rep("rmx", 3), 
                                  paste("fs-corrected rmx estimate for 'eps' in [", 
                                    round(eps.lower, 3), ", ", round(eps.upper, 3), "]", sep = ""),
                                  paste("least favorable (uncorrected) contamination: ", 
                                        100*signif(r.as/sqrtn, 3), " %", sep = ""),
                                  paste("maximum asymptotic MSE-inefficiency: ", signif(ineff, 3), sep = "")), 
                                  ncol = 2, dimnames = list(NULL, c("method", "message")))
        }else{
            Info.matrix <- matrix(c(rep("rmx", 3), 
                                  paste("rmx estimate for 'eps' in [", 
                                    round(eps.lower, 3), ", ", round(eps.upper, 3), "]", sep = ""),
                                  paste("least favorable contamination: ", 
                                        100*signif(r/sqrtn, 3), " %", sep = ""),
                                  paste("maximum asymptotic MSE-inefficiency: ", 
                                        signif(ineff, 3), sep = "")), 
                                  ncol = 2, dimnames = list(NULL, c("method", "message")))
        }
        IF <- rmxEst.all$IF
        IF$radius <- r
    }

    RMX <- list(rmxEst = rmxEst, rmxIF = IF, initial.est = lambda, 
                Infos = Info.matrix, n = length(x))
    class(RMX) <- "rmx"
    RMX
}

###############################################################################
## computation of radius-minimax IC
###############################################################################
.getInterval.pois <- function(r, rlo, rup, lambda, aUp, cUp, delta){
    lcr <- .lcr.pois(lambda = lambda)

    LM <- .getLM.pois(r0 = r, lambda = lambda, aUp = aUp, 
                       cUp = cUp, delta = delta)
    A <- LM$A
    c0 <- LM$c0
    b <- A*c0
    
    if(r >= lcr){
        effre <- 1
    }else{
        if(rup >= lcr){
            m0 <- qpois(.5,lambda)
            k0 <- m0/lambda
            m1 <- ceiling(m0)
            m2 <- floor(m0)
            p1 <- 1-ppois(m2-1,lambda)
            p2 <- 1-ppois(m2,lambda)
            p3 <- ppois(m1-1,lambda)
            p4 <- ppois(m1-2,lambda)
            
            inte <- p1 - p4 + k0*(p3-p2)
            bmin <- 1/inte
            names(bmin) <- NULL

            effre <- b^2/bmin^2
        }
        else{
            LMre <- .getLM.pois(r0 = rup, lambda = lambda, aUp = aUp, 
                                cUp = cUp, delta = delta)
            effre <- (A - b^2*(r^2 - rup^2))/LMre$A
        }
    }
    
    if(identical(all.equal(rlo, 0), TRUE))
        effli <- (A - b^2*r^2)/lambda
    else{
        LMli <- .getLM.pois(r0 = rlo, lambda = lambda, aUp = aUp, 
                            cUp = cUp, delta = delta)
        effli <- (A - b^2*(r^2 - rlo^2))/LMli$A
    }
    
    effre-effli
}
###############################################################################
## computation of k-step construction
###############################################################################
.onestep.pois <- function(x, IF){
    lambda <- IF$parameter["lambda"]
    max(.Machine$double.eps, lambda + mean(IF$IFun(x)))
}
.updateIF.pois <- function(lambda, IF, aUp, cUp, delta){
    r <- IF$radius
    if(r == 0){
        IF <- .getMLIF.pois(lambda = lambda)
        IF$radius <- r
        return(IF)
    }
    lcr <- .lcr.pois(lambda = lambda)
    if(r >= lcr){
        IF <- .getMBIF.pois(lambda = lambda)
        IF$radius <- r
        if(r != Inf){
            IF$asMSE <- IF$asVar + r^2*IF$b^2
            IF$asBias <- r^2*IF$b^2
        }
        return(IF)
    }
    
    LM <- .getLM.pois(r0 = r, lambda = lambda, 
                      aUp = aUp, cUp = cUp, delta = delta)
    IF <- .getOptIF.pois(lambda = lambda, radius = r, 
                         A = LM$A, a = LM$A*LM$z, b = LM$A*LM$c0)
    IF$radius <- r
    IF
}
.kstep.pois <- function(x, IF, k, aUp, cUp, delta){
    lambda <- .onestep.pois(x = x, IF = IF)
    if(k > 1){
        for(i in 2:k){
            IF <- .updateIF.pois(lambda = lambda, IF = IF, aUp = aUp, cUp = cUp,
                                 delta = delta)
            lambda <- .onestep.pois(x = x, IF = IF)
        }
    }
    IF <- .updateIF.pois(lambda = lambda, IF = IF, aUp = aUp, cUp = cUp,
                         delta = delta)

    list(lambda = lambda, IF = IF)
}
