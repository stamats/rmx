###############################################################################
## Function for finite-sample correction of the neighborhood radius
###############################################################################
fsRadius <- function(r, n, model = "norm", ...){
    stopifnot(is.numeric(r))
    stopifnot(length(r) == 1)
    stopifnot(r >= 0)
    if(r == 0) return(r)
    stopifnot(is.numeric(n))
    n <- as.integer(n)
    stopifnot(n > 0)
    if(n == 1) return(Inf)
    if(n == 2) return(Inf)

    if(model == "norm"){
        r <- fsRadius.norm(r = r, n = n)
        return(r)
    }
    if(model == "binom"){
        listDots <- list(...)
        if(!"prob" %in% names(listDots))
            stop("Parameter 'prob' must be specified!")
        prob <- listDots$prob
        if(!"size" %in% names(listDots))
            stop("Parameter 'size' must be specified!")
        size <- listDots$size
        M <- ifelse("M" %in% names(listDots), listDots$M, 1e4)
        parallel <- ifelse("parallel" %in% names(listDots), listDots$parallel, FALSE)
        if("ncores" %in% names(listDots)){
            ncores <- listDots$ncores
        }else{
            ncores <- NULL
        }
        r <- fsRadius.binom(r = r, n = n, prob = prob, size = size, M = M,
                            parallel = parallel, ncores = ncores)
        return(r)
    }
    warning("Finite-sample correction not yet implemented, asymptotic radius is used!")
    r
}
fsRadius.norm <- function(r, n){
    if(r >= 1.67) return(r)
    
    eps <- r/sqrt(n)
    ns <- c(3:50, seq(55, 100, by = 5), seq(110, 200, by = 10), 
            seq(250, 500, by = 50))
    ns <- as.integer(ns)
    epss <- c(seq(0.001, 0.01, by = 0.001), seq(0.02, to = 0.5, by = 0.01))
    if(n %in% ns){
        ind <- ns == n
        r <- max(r, approx(x = epss, y = .fsRadius.norm[,ind], 
                           xout = eps, rule = 2)$y)
    }else{
        if(n > 500){
            ind <- ns == 500L
            r1 <- approx(x = epss, y = .fsRadius.norm[,ind], 
                         xout = eps, rule = 2)$y
            r <- max(r, mean(c(r1, r)))
        }else{
            ind <- sort(order(abs(ns-n))[1:2])
            r1 <- approx(x = epss, y = .fsRadius.norm[,ind[1]], 
                         xout = eps, rule = 2)$y
            r2 <- approx(x = epss, y = .fsRadius.norm[,ind[2]], 
                         xout = eps, rule = 2)$y
            D <- ns[ind[2]] - ns[ind[1]]
            r <- max(r, (n-ns[ind[1]])/D*r1 + (ns[ind[2]]-n)/D*r2)
        }
    }
    r  
}
fsRadius.binom <- function(r, n, prob, size, M = 10000, parallel = FALSE, 
                           ncores = NULL){
    if(prob <= 0.5){
        D <- size  
    }else{
        D <- 0
    }
    
    lcr <- .lcr.binom(prob = prob, size = size)
    rs <- seq(from = r, to = lcr, length.out = 10)
    eps <- rs/sqrt(n)
    res <- numeric(M)
    RMX <- function(x, eps, size){
        rmx.binom(x, eps = eps, initial.est = NULL, 
                  k = 3, fsCor = FALSE, size = size)$rmxEst
    }
    mse <- numeric(length(eps))
    if(parallel){
        if(is.null(ncores)){
            cores <- detectCores()
            cl <- makeCluster(cores[1]-1)
        }else{
            cl <- makeCluster(ncores)
        }
        clusterExport(cl, list("rmx.binom", "cvm.binom", ".cvmdist", 
                               ".lcr.binom", ".getMBIF.binom", 
                               ".getLM.binom", ".getc.binom", 
                               ".geta.binom", ".getOptIF.binom",
                               ".kstep.binom", ".onestep.binom",
                               ".updateIF.binom"),
                      envir = environment(fun = fsRadius.binom))
    }
    for(i in seq_along(eps)){
        ind <- rbinom(n*M, size = 1, prob = eps[i])
        Mat <- matrix((1-ind)*rbinom(n*M, size = size, prob = prob) + ind*D, 
                      ncol = n)
        if(parallel){
            res <- parApply(cl = cl, X = Mat, MARGIN = 1, FUN = RMX, 
                            eps = eps[i], size = size)
        }else{
            res <- apply(Mat, 1, RMX, eps = eps[i], size = size)
        }
        mse[i] <- n*mean((res - prob)^2)
    }
    if(parallel) stopCluster(cl)
    r <- rs[which.min(mse)]
    r
}
