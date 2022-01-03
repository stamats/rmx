###############################################################################
## RMX estimator for probability of success of a binomial model
###############################################################################
rmx.binom <- function(x, eps.lower, eps.upper, eps, initial.est = NULL, k = 3L, 
                      fsCor = FALSE, size, M = 10000, parallel = FALSE, ncores = NULL, 
                      aUp = 100*size, cUp = 1e4, delta = 1e-9){
    if(length(x) <= 2){
        stop("A sample size of at least 3 is required!")
    }
    if(is.null(initial.est)){
        prob <- cvm.binom(x, size = size, mu = "model")
    }else{
        stopifnot(is.numeric(initial.est))
        if(length(initial.est) != 1)
            stop("'initial.est' needs to be a numeric vector of length 1 or NULL")
        if(0 < initial.est || initial.est < 1)
            stop("initial estimate for prob <= 0 or >=1 which is not valid")
        prob <- initial.est
    }

    lcr <- .lcr.binom(prob = prob, size = size)
    if(!is.null(eps)){
        r <- sqrt(length(x))*eps
        if(fsCor){ 
            r.as <- r
            if(r < lcr){
                r <- fsRadius.binom(r = r, n = length(x), prob = prob,
                                    size = size, M = M, parallel = parallel, 
                                    ncores = ncores)
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
            r <- uniroot(.getInterval.binom, lower = rlo+1e-8, upper = rup, 
                         tol = .Machine$double.eps^0.25, rlo = rlo, 
                         rup = rup, prob = prob, size = size, aUp = aUp, 
                         cUp = cUp, delta = delta)$root
            r.as <- r
        }
        if(fsCor){
            r.as <- r
            if(r < lcr){
                r <- fsRadius.binom(r = r, n = length(x), prob = prob,
                                    size = size, M = M, parallel = parallel, 
                                    ncores = ncores)
            }
        }
    }
    if(!is.null(eps)){
        if(eps == 0){
            rmxEst <- mean(x)/size
            names(rmxEst) <- "prob"
            IF <- .getMLIF.binom(size = size, prob = rmxEst)
            IF$radius <- 0
            Info.matrix <- matrix(c("rmx", "ML estimate"), ncol = 2, 
                                  dimnames = list(NULL, c("method", "message")))

            RMX <- list(rmxEst = rmxEst, rmxIF = IF, initial.est = NULL, 
                        Infos = Info.matrix)
            class(RMX) <- "rmx"
            return(RMX)
        }
    }
    if(!is.null(eps)){
        if(r >= lcr){
            IF <- .getMBIF.binom(size = size, prob = prob)
            IF$radius <- r
        }else{
            LM <- .getLM.binom(r0 = r, prob = prob, size = size, 
                               aUp = aUp, cUp = cUp, delta = delta)
            IF <- .getOptIF.binom(size = size, prob = prob, radius = r, 
                                  A = LM$A, a = LM$A*LM$z, b = LM$A*LM$c0)
            IF$radius <- r
        }
        rmxEst.all <- .kstep.binom(x = x, IF, k = k, aUp = aUp, cUp = cUp,
                                   delta = delta)
        rmxEst <- rmxEst.all$prob
        names(rmxEst) <- "prob"
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
            IF <- .getMBIF.binom(size = size, prob = prob)
            IF$radius <- r
        }else{
            LM <- .getLM.binom(r0 = r, prob = prob, size = size, 
                               aUp = aUp, cUp = cUp, delta = delta)
            IF <- .getOptIF.binom(size = size, prob = prob, radius = r, 
                                  A = LM$A, a = LM$A*LM$z, b = LM$A*LM$c0)
            IF$radius <- r
        }
        if(rlo == 0){
            b <- LM$A*LM$c0
            ineff <- (LM$A - b^2*r.as^2)/(prob*(1-prob))*size
        }else{
            if(rlo >= lcr){
                ineff <- 1
            }else{
                b <- LM$A*LM$c0
                ineff <- (LM$A - b^2*(r.as^2 - rlo^2))/LM$A
            }
        }
        rmxEst.all <- .kstep.binom(x = x, IF = IF, k = k, aUp = aUp, cUp = cUp,
                                   delta = delta)
        rmxEst <- rmxEst.all$prob
        names(rmxEst) <- "prob"
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

    RMX <- list(rmxEst = rmxEst, rmxIF = IF, initial.est = prob, 
                Infos = Info.matrix, n = length(x))
    class(RMX) <- "rmx"
    RMX
}

###############################################################################
## computation of radius-minimax IC
## using predefined functions included in "sysdata.rda"
###############################################################################
.getInterval.binom <- function(r, rlo, rup, prob, size, aUp, cUp, delta){
    LM <- .getLM.binom(r0 = r, prob = prob, size = size, aUp = aUp, 
                       cUp = cUp, delta = delta)
    A <- LM$A
    c0 <- LM$c0
    b <- A*c0
    
    lcr <- .lcr.binom(prob = prob, size = size)
    
    if(rup >= lcr){
        m0 <- qbinom(.5, size=size, prob=prob)
        p1 <- pbinom(m0, size=size, prob=prob)
        p2 <- pbinom(m0-1, size=size-1, prob=prob)
        
        inte <- m0*(2*p1 - 1) + size*prob*(1-2*p2)
        bmin <- prob*(1-prob)/inte
        
        effre <- b^2/bmin^2
    }
    else{
        LMre <- .getLM.binom(r0 = rup, prob = prob, size = size, aUp = 100*size, 
                             cUp = 1e4, delta = delta)
        effre <- (A - b^2*(r^2 - rup^2))/LMre$A
    }
    
    if(identical(all.equal(rlo, 0), TRUE))
        effli <- (A - b^2*r^2)/(prob*(1-prob))*size
    else{
        LMli <- .getLM.binom(r0 = rlo, prob = prob, size = size, aUp = 100*size, 
                             cUp = 1e4, delta = delta)
        effli <- (A - b^2*(r^2 - rlo^2))/LMli$A
    }
    
    effre-effli
}
###############################################################################
## computation of k-step construction
###############################################################################
.onestep.binom <- function(x, IF){
    prob <- IF$parameter["prob"]
    min(max(.Machine$double.eps, prob + mean(IF$IFun(x))), 1-.Machine$double.eps)
}
.updateIF.binom <- function(prob, IF, aUp, cUp, delta){
    r <- IF$radius
    if(r == 0){
        IF <- .getMLIF.binom(size = IF$parameter["size (known)"],
                             prob = prob)
        IF$radius <- r
        return(IF)
    }
    lcr <- .lcr.binom(prob = prob, size = IF$parameter["size (known)"])
    if(r >= lcr){
        IF <- .getMBIF.binom(size = IF$parameter["size (known)"],
                             prob = prob)
        IF$radius <- r
        return(IF)
    }
    
    LM <- .getLM.binom(r0 = r, prob = prob, size = IF$parameter["size (known)"], 
                       aUp = aUp, cUp = cUp, delta = delta)
    IF <- .getOptIF.binom(size = IF$parameter["size (known)"],
                          prob = prob, radius = r, 
                          A = LM$A, a = LM$A*LM$z, b = LM$A*LM$c0)
    IF$radius <- r
    IF
}
.kstep.binom <- function(x, IF, k, aUp, cUp, delta){
    prob <- .onestep.binom(x = x, IF = IF)
    if(k > 1){
        for(i in 2:k){
            IF <- .updateIF.binom(prob = prob, IF = IF, aUp = aUp, cUp = cUp,
                                  delta = delta)
            prob <- .onestep.binom(x = x, IF = IF)
        }
    }
    IF <- .updateIF.binom(prob = prob, IF = IF, aUp = aUp, cUp = cUp,
                          delta = delta)

    list(prob = prob, IF = IF)
}
