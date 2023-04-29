###############################################################################
## RMX estimator for probability of success of a binomial model
###############################################################################
rowRmx.pois <- function(x, eps.lower = 0, eps.upper, eps = NULL, initial.est = NULL, 
                        k = 3L, fsCor = FALSE, na.rm = TRUE, computeSE = FALSE, 
                        parallel = FALSE, ncores = NULL, aUp = 100*max(x), 
                        cUp = 1e4, delta = 1e-9){
    if(!is.null(eps)){
        if(eps == 0){
            rmxEst <- matrix(rowMeans(x), ncol = 1)
            colnames(rmxEst) <- "lambda"
            
            Info.matrix <- matrix(c("rowRmx", 
                                    paste("ML estimate")),
                                  ncol = 2, dimnames = list(NULL, c("method", "message")))
            
            if(computeSE){
                asSE <- matrix(sqrt(rmxEst/ncol(x)), ncol = 1)
                colnames(asSE) <- "SE.lambda"
            }else{
                asSE <- NA
            }
            
            RMX <- list(model = "pois", modelName = "Poisson mean",
                        rmxEst = rmxEst, asSE = asSE, radius = 0, 
                        Infos = Info.matrix)
            class(RMX) <- "RMXlist"
            return(RMX)
        }
    }
    if(ncol(x) <= 2){
        stop("A sample size of at least 3 is required!")
    }
    
    if(is.null(initial.est)){
        prob <- rowCVM(x, model = "pois", size = size, 
                       parallel = parallel, ncores = ncores)
    }else{
        stopifnot(is.numeric(initial.est))
        if(is.matrix(initial.est)){
            if(nrow(initial.est) != nrow(x) || ncol(initial.est) != 1)
                stop("'initial.est' has wrong dimension")
            lambda <- initial.est[,1]
        }else{
            if(length(initial.est) != nrow(x))
                stop("Length of 'initial.est' not equal to 'nrow(x)'")
            lambda <- initial.est
        }
    }
    
    lcr <- .lcr.pois(lambda = median(lambda))
    if(!is.null(eps)){
        r <- sqrt(ncol(x))*eps
        if(fsCor){ 
            r.as <- r
            if(r < lcr){
                r <- fsRadius.pois(r = r, n = ncol(x), 
                                   lambda = median(lambda),
                                   parallel = parallel, ncores = ncores)
            }
        }
    }else{
        sqrtn <- sqrt(ncol(x))
        rlo <- sqrtn*eps.lower
        rup <- sqrtn*eps.upper
        if(rlo >= lcr){
            r <- (rlo + rup)/2
            r.as <- r
        }else{
            r <- uniroot(.getInterval.pois, lower = rlo+1e-8, upper = rup, 
                         tol = .Machine$double.eps^0.25, rlo = rlo, 
                         rup = rup, prob = median(prob),
                         aUp = aUp, cUp = cUp, delta = delta)$root
            r.as <- r
        }
        if(fsCor){
            r.as <- r
            if(r < lcr){
                r <- fsRadius.pois(r = r, n = ncol(x), 
                                   lambda = median(lambda),
                                   parallel = parallel, ncores = ncores)
            }
        }
    }
    if(!is.null(eps)){
        if(r >= lcr){
            LM <- .getMBLM.pois(lambda = lambda)
        }else{
            LM <- .getLM.pois.vector(r0 = r, lambda = lambda, 
                                     parallel = parallel, ncores = ncores, 
                                     aUp = aUp, cUp = cUp, delta = delta)
        }
        rmxEst.all <- .kstep.pois.matrix(x = x, r = r, LM = LM, k = k, 
                                         lambda = lambda, na.rm = na.rm,
                                         parallel = parallel, ncores = ncores, 
                                         aUp = aUp, cUp = cUp, delta = delta)
        rmxEst <- rmxEst.all$lambda
        colnames(rmxEst) <- "lambda"
        if(computeSE){
            asVar <- rmxEst.all$LM$asVar
            asSE <- matrix(sqrt(asVar/ncol(x)), ncol = 1)
            colnames(asSE) <- "SE.lambda"
        }else{
            asSE <- NA
        }
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
    }else{
        if(r >= lcr){
            LM <- .getMBLM.pois(lambda = lambda)
        }else{
            LM <- .getLM.pois.vector(r0 = r, lambda = lambda,
                                     parallel = parallel, ncores = ncores, 
                                     aUp = aUp, cUp = cUp, delta = delta)
        }
        if(rlo == 0){
            ineff <- median((LM$A - LM$b^2*r.as^2)/lambda)
        }else{
            if(rlo >= lcr){
                ineff <- 1
            }else{
                ineff <- median((LM$A - LM$b^2*(r.as^2 - rlo^2))/LM$A)
            }
        }
        rmxEst.all <- .kstep.pois.matrix(x = x, r = r, LM = LM, k = k, 
                                         lambda = lambda, na.rm = na.rm,
                                         parallel = parallel, ncores = ncores, 
                                         aUp = aUp, cUp = cUp, delta = delta)
        rmxEst <- rmxEst.all$lambda
        colnames(rmxEst) <- "lambda"
        if(computeSE){
            asVar <- rmxEst.all$LM$A - r^2*rmxEst.all$LM$b^2
            asSE <- matrix(sqrt(asVar/ncol(x)), ncol = 1)
            colnames(asSE) <- "SE.lambda"
        }else{
            asSE <- NA
        }
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
    }
    RMX <- list(model = "pois", modelName = "Poisson mean", 
                rmxEst = rmxEst, asSE = asSE, radius = r, 
                Infos = Info.matrix)
    class(RMX) <- "RMXlist"
    RMX
}

##################################################################
## Lagrange multipliers of minimum bias estimator
##################################################################
.getMBLM.pois <- function(lambda){
    m0 <- qpois(.5,lambda)
    k0 <- m0/lambda
    m1 <- ceiling(m0)
    m2 <- floor(m0)
    p1 <- 1-ppois(m2-1,lambda)
    p2 <- 1-ppois(m2,lambda)
    p3 <- ppois(m1-1,lambda)
    p4 <- ppois(m1-2,lambda)
    
    inte <- p1 - p4 + k0*(p3-p2)
    b <- 1/inte
    names(b) <- NULL
    
    m0 <- qpois(.5,lambda)
    p <- ppois(m0-1,lambda)
    d <- dpois(m0,lambda)
    beta <- (2*p+d-1)/d
    
    A <- 1
    a <- m0/lambda - 1
    names(a) <- NULL
    
    list(A = A, a = a, b = b, beta = beta, asVar = b^2, MB = TRUE)
}
##################################################################
## Lagrange multipliers of minimum bias estimator
##################################################################
.getLM.pois.vector <- function(r0, lambda, parallel, ncores, 
                                aUp = 100*max(lambda), cUp = 1e4, delta = 1e-9){
    fun <- function(lambda, r0, aUp, cUp, delta){
        res <- .getLM.pois(r0 = r0, lambda = lambda, aUp = aUp, 
                            cUp = cUp, delta = delta)
        list(A = res$A, a = res$A*res$z, b = res$A*res$c0)
    }
    if(parallel){
        if(is.null(ncores)){
            cores <- detectCores()
            cl <- makeCluster(cores[1]-1)
        }else{
            cl <- makeCluster(ncores)
        }
        clusterExport(cl, list(".getLM.pois", ".getc.pois", ".geta.pois"),
                      envir = environment(fun = .getLM.pois.vector))
        LM <- parSapply(cl = cl, X = lambda, FUN = fun, r0 = r0, 
                        aUp = aUp, cUp = cUp, delta = delta)  
        stopCluster(cl)
    }else{
        LM <- sapply(lambda, fun, r0 = r0, aUp = aUp, cUp = cUp, 
                     delta = delta)
    }
    list(A = unlist(LM["A",]), a = unlist(LM["a",]), b = unlist(LM["b",]),
         asVar = unlist(LM["A",]) - r0^2*unlist(LM["b",])^2, MB = FALSE)
}
###############################################################################
## computation of k-step construction
###############################################################################
.onestep.pois.matrix <- function(x, LM, lambda, na.rm){
    if(LM$MB){
        M <- qpois(0.5, lambda = lambda)
        IFx <- rowMeans(LM$b*((x > M) - (x < M) + LM$beta*(x == M)), na.rm = na.rm)
    }else{
        Y <- LM$A*(x/lambda - 1) - LM$a
        ind <- LM$b <= abs(Y)
        IFx <- rowMeans(Y*(ind*LM$b/abs(Y) + !ind), na.rm = na.rm)
    }
        
    pmax(.Machine$double.eps, prob + IFx)
}
.updateLM.pois.matrix <- function(lambda, LM, r, parallel, ncores, 
                                  aUp = 100*max(lambda), cUp = 1e4, delta = 1e-9){
    if(r == 0){
        A <- lambda
        a <- rep(0, length(A))
        b <- rep(Inf, length(A))
        return(list(A = A, a = a, b = b, asVar = A))
    }
    lcr <- .lcr.pois(lambda = median(lambda))
    if(r >= lcr){
        LM <- .getMBLM.pois(lambda = lambda)
        return(LM)
    }
    LM <- .getLM.pois.vector(r0 = r, lambda = lambda, 
                             parallel = parallel, ncores = ncores, 
                             aUp = aUp, cUp = cUp, delta = delta)
    LM
}
.kstep.pois.matrix <- function(x, r, LM, lambda, na.rm, k, parallel, ncores, 
                                aUp = 100*max(lambda), cUp = 1e4, delta = 1e-9){
    lambda <- .onestep.pois.matrix(x = x, LM = LM, lambda = lambda, 
                                 na.rm = na.rm)
    if(k > 1){
        for(i in 2:k){
            LM <- .updateLM.pois.matrix(lambda = lambda, LM = LM, r = r,
                                        parallel = parallel, ncores = ncores, 
                                        aUp = aUp, cUp = cUp, delta = delta)
            prob <- .onestep.pois.matrix(x = x, LM = LM, lambda = lambda,
                                          na.rm = na.rm)
        }
    }
    LM <- .updateLM.pois.matrix(lambda = lambda, LM = LM, r = r,
                                 parallel = parallel, ncores = ncores, 
                                 aUp = aUp, cUp = cUp, delta = delta)

    list(lambda = matrix(lambda, ncol = 1), LM = LM)
}
