###############################################################################
## RMX estimator for probability of success of a binomial model
###############################################################################
rowRmx.binom <- function(x, eps.lower = 0, eps.upper, eps = NULL, initial.est = NULL, 
                         k = 3L, fsCor = FALSE, na.rm = TRUE, size, computeSE = FALSE, 
                         parallel = FALSE, ncores = NULL, aUp = 100*size, 
                         cUp = 1e4, delta = 1e-9){
    if(!is.null(eps)){
        if(eps == 0){
            rmxEst <- matrix(rowMeans(x)/size, ncol = 1)
            colnames(rmxEst) <- "prob"
            
            Info.matrix <- matrix(c("rowRmx", 
                                    paste("ML estimate")),
                                  ncol = 2, dimnames = list(NULL, c("method", "message")))
            
            if(computeSE){
                asSE <- matrix(sqrt((rmxEst*(1-rmxEst))/size)/sqrt(ncol(x)), 
                               ncol = 1)
                colnames(asSE) <- "SE.prob"
            }else{
                asSE <- NA
            }
            
            RMX <- list(model = "binom", modelName = "binomial probability",
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
        prob <- rowCVM(x, model = "binom", size = size, 
                       parallel = parallel, ncores = ncores)
    }else{
        stopifnot(is.numeric(initial.est))
        if(is.matrix(initial.est)){
            if(nrow(initial.est) != nrow(x) || ncol(initial.est) != 1)
                stop("'initial.est' has wrong dimension")
            prob <- initial.est[,1]
        }else{
            if(length(initial.est) != nrow(x))
                stop("Length of 'initial.est' not equal to 'nrow(x)'")
            prob <- initial.est
        }
    }
    
    lcr <- .lcr.binom(prob = median(prob), size = size)
    if(!is.null(eps)){
        r <- sqrt(ncol(x))*eps
        if(fsCor){ 
            r.as <- r
            if(r < lcr){
                r <- fsRadius.binom(r = r, n = ncol(x), 
                                    size = size, prob = median(prob),
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
            r <- uniroot(.getInterval.binom, lower = rlo+1e-8, upper = rup, 
                         tol = .Machine$double.eps^0.25, rlo = rlo, 
                         rup = rup, prob = median(prob), size = size, 
                         aUp = aUp, cUp = cUp, delta = delta)$root
            r.as <- r
        }
        if(fsCor){
            r.as <- r
            if(r < lcr){
                r <- fsRadius.binom(r = r, n = ncol(x), 
                                    size = size, prob = median(prob),
                                    parallel = parallel, ncores = ncores)
            }
        }
    }
    if(!is.null(eps)){
        if(r >= lcr){
            LM <- .getMBLM.binom(size = size, prob = prob)
        }else{
            LM <- .getLM.binom.vector(r0 = r, prob = prob, size = size, 
                                      parallel = parallel, ncores = ncores, 
                                      aUp = aUp, cUp = cUp, delta = delta)
        }
        rmxEst.all <- .kstep.binom.matrix(x = x, r = r, LM = LM, k = k, 
                                      prob = prob, na.rm = na.rm, size = size,
                                      parallel = parallel, ncores = ncores, 
                                      aUp = aUp, cUp = cUp, delta = delta)
        rmxEst <- rmxEst.all$prob
        colnames(rmxEst) <- "prob"
        if(computeSE){
            asVar <- rmxEst.all$LM$asVar
            asSE <- matrix(sqrt(asVar/ncol(x)), ncol = 1)
            colnames(asSE) <- "SE.prob"
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
            LM <- .getMBLM.binom(size = size, prob = prob)
        }else{
            LM <- .getLM.binom.vector(r0 = r, prob = prob, size = size, 
                                      parallel = parallel, ncores = ncores, 
                                      aUp = aUp, cUp = cUp, delta = delta)
        }
        if(rlo == 0){
            ineff <- median((LM$A - LM$b^2*r.as^2)/(prob*(1-prob))*size)
        }else{
            if(rlo >= lcr){
                ineff <- 1
            }else{
                ineff <- median((LM$A - LM$b^2*(r.as^2 - rlo^2))/LM$A)
            }
        }
        rmxEst.all <- .kstep.binom.matrix(x = x, r = r, LM = LM, k = k, 
                                      prob = prob, na.rm = na.rm, size = size,
                                      parallel = parallel, ncores = ncores, 
                                      aUp = aUp, cUp = cUp, delta = delta)
        rmxEst <- rmxEst.all$prob
        colnames(rmxEst) <- "prob"
        if(computeSE){
            asVar <- rmxEst.all$LM$A - r^2*rmxEst.all$LM$b^2
            asSE <- matrix(sqrt(asVar/ncol(x)), ncol = 1)
            colnames(asSE) <- "SE.prob"
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
    RMX <- list(model = "binom", modelName = "binomial probability", 
                rmxEst = rmxEst, asSE = asSE, radius = r, 
                Infos = Info.matrix)
    class(RMX) <- "RMXlist"
    RMX
}

##################################################################
## Lagrange multipliers of minimum bias estimator
##################################################################
.getMBLM.binom <- function(prob, size){
    m0 <- qbinom(0.5, size=size, prob=prob)
    p1 <- pbinom(m0, size=size, prob=prob)
    p2 <- pbinom(m0-1, size=size-1, prob=prob)
    
    inte <- m0*(2*p1 - 1) + size*prob*(1-2*p2)
    b <- prob*(1-prob)/inte
    
    p3 <- pbinom(m0-1, size = size, prob = prob)
    p4 <- pbinom(m0, size = size, prob = prob, lower.tail = FALSE)
    p5 <- dbinom(m0, size = size, prob = prob)
    beta <- (p3 - p4)/p5
    
    A <- rep(1, length(prob))
    a <- -size*prob/(prob*(1-prob))

    list(A = A, a = a, b = b, beta = beta, asVar = b^2, MB = TRUE)
}
##################################################################
## Lagrange multipliers of minimum bias estimator
##################################################################
.getLM.binom.vector <- function(r0, prob, size, parallel, ncores, 
                                aUp = 100*size, cUp = 1e4, delta = 1e-9){
    fun <- function(prob, r0, size, aUp, cUp, delta){
        res <- .getLM.binom(r0 = r0, prob = prob, size = size, aUp = aUp, 
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
        clusterExport(cl, list(".getLM.binom", ".getc.binom", ".geta.binom"),
                      envir = environment(fun = .getLM.binom.vector))
        LM <- parSapply(cl = cl, X = prob, FUN = fun, r0 = r0, size = size, 
                        aUp = aUp, cUp = cUp, delta = delta)  
        stopCluster(cl)
    }else{
        LM <- sapply(prob, fun, size = size, r0 = r0, aUp = aUp, cUp = cUp, 
                     delta = delta)
    }
    list(A = unlist(LM["A",]), a = unlist(LM["a",]), b = unlist(LM["b",]),
         asVar = unlist(LM["A",]) - r0^2*unlist(LM["b",])^2, MB = FALSE)
}
###############################################################################
## computation of k-step construction
###############################################################################
.onestep.binom.matrix <- function(x, LM, prob, size, na.rm){
    if(LM$MB){
        M <- qbinom(0.5, prob = prob, size = size)
        IFx <- rowMeans(LM$b*((x > M) - (x < M) + LM$beta*(x == M)), na.rm = na.rm)
    }else{
        Y <- LM$A*(x-size*prob)/(prob*(1-prob)) - LM$a
        ind <- LM$b <= abs(Y)
        IFx <- rowMeans(Y*(ind*LM$b/abs(Y) + !ind), na.rm = na.rm)
    }
        
    pmin(pmax(.Machine$double.eps, prob + IFx), 1-.Machine$double.eps)
}
.updateLM.binom.matrix <- function(prob, LM, size, r, parallel, ncores, 
                                   aUp = 100*size, cUp = 1e4, delta = 1e-9){
    if(r == 0){
        A <- (prob*(1-prob))/size
        a <- rep(0, length(A))
        b <- rep(Inf, length(A))
        return(list(A = A, a = a, b = b, asVar = A))
    }
    lcr <- .lcr.binom(size = size, prob = median(prob))
    if(r >= lcr){
        LM <- .getMBLM.binom(prob = prob, size = size)
        return(LM)
    }
    LM <- .getLM.binom.vector(r0 = r, prob = prob, size = size, 
                              parallel = parallel, ncores = ncores, 
                              aUp = aUp, cUp = cUp, delta = delta)
    LM
}
.kstep.binom.matrix <- function(x, r, LM, prob, na.rm, size, k, parallel, ncores, 
                                aUp = 100*size, cUp = 1e4, delta = 1e-9){
    prob <- .onestep.binom.matrix(x = x, LM = LM, prob = prob, size = size, 
                                  na.rm = na.rm)
    if(k > 1){
        for(i in 2:k){
            LM <- .updateLM.binom.matrix(prob = prob, LM = LM, size = size, r = r,
                                         parallel = parallel, ncores = ncores, 
                                         aUp = aUp, cUp = cUp, delta = delta)
            prob <- .onestep.binom.matrix(x = x, LM = LM, prob = prob, size = size,
                                          na.rm = na.rm)
        }
    }
    LM <- .updateLM.binom.matrix(prob = prob, LM = LM, size = size, r = r,
                                 parallel = parallel, ncores = ncores, 
                                 aUp = aUp, cUp = cUp, delta = delta)

    list(prob = matrix(prob, ncol = 1), LM = LM)
}
