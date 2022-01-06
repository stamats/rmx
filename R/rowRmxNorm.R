###############################################################################
## Evaluate roblox on rows of a matrix
###############################################################################
rowRmx.norm <- function(x, eps.lower = 0, eps.upper, eps = NULL, initial.est = NULL, 
                        k = 3L, fsCor = TRUE, na.rm = TRUE, computeSE = TRUE){
    mad0 <- 1e-4
    if(!is.null(eps)){
        r <- sqrt(ncol(x))*eps # missing values per row are not considered!
        if(fsCor){ 
            r.as <- r
            r <- fsRadius.norm(r = r, n = ncol(x))
        }
    }else{
        sqrtn <- sqrt(ncol(x))
        rlo <- sqrtn*eps.lower
        rup <- sqrtn*eps.upper
        if(rlo > 10){
            r <- (rlo + rup)/2
            r.as <- r
        }else{
            r <- uniroot(.getInterval.norm, lower = rlo+1e-8, upper = rup, 
                         tol = .Machine$double.eps^0.25, rlo = rlo, rup = rup)$root
            r.as <- r
        }
        if(fsCor){
            r.as <- r
            r <- fsRadius.norm(r = r, n = ncol(x))
        }
    }
    if(ncol(x) <= 2){
        warning("Sample size <= 2! => Median and MAD are used for estimation.")
        MEAN <- rowMedians(x, na.rm = na.rm)
        SD <- rowMedians(abs(x-MEAN), na.rm = na.rm)/qnorm(0.75)
        rmxEst <- cbind(MEAN, SD)
        colnames(rmxEst) <- c("mean", "sd")

        Info.matrix <- matrix(c("rowRmx", 
                                paste("median and MAD")),
                              ncol = 2, dimnames = list(NULL, c("method", "message")))
        
        if(computeSE){
            b1 <- SD*sqrt(pi/2)
            b2 <- SD/(4*qnorm(0.75)*dnorm(qnorm(0.75)))
            SE1 <- b1/sqrt(ncol(x))
            SE2 <- b2/sqrt(ncol(x))
            asSE <- cbind(SE1, SE2)
            colnames(asSE) <- c("SE.mean", "SE.sd")
        }else{
            asSE <- NA
        }

        RMX <- list(model = "norm", modelName = "normal location and scale",
                    rmxEst = rmxEst, asSE = asSE, radius = r, 
                    Infos = Info.matrix)
        class(RMX) <- "RMXlist"
        return(RMX)
    }
    
    if(!is.null(eps)){
        if(eps == 0){
            MEAN <- rowMeans(x, na.rm = na.rm)
            n <- rowSums(!is.na(x))
            n[n < 1] <- NA
            SD <- sqrt(rowSums((x - MEAN)^2, na.rm = na.rm)/n)
            rmxEst <- cbind(MEAN, SD)

            colnames(rmxEst) <- c("mean", "sd")
            Info.matrix <- matrix(c("rowRmx", 
                                  paste("mean and sd")),
                                  ncol = 2, dimnames = list(NULL, c("method", "message")))
            
            if(computeSE){
                SE1 <- SD/sqrt(ncol(x))
                SE2 <- 0.5*SD/sqrt(ncol(x))
                asSE <- cbind(SE1, SE2)
                colnames(asSE) <- c("SE.mean", "SE.sd")
            }else{
                asSE <- NA
            }
            
            RMX <- list(model = "norm", modelName = "normal location and scale",
                        rmxEst = rmxEst, asSE = asSE, radius = r, 
                        Infos = Info.matrix)
            class(RMX) <- "RMXlist"
            return(RMX)
        }
    }

    if(is.null(initial.est)){
        MEAN <- rowMedians(x, na.rm = na.rm)
        SD <- rowMedians(abs(x-MEAN), na.rm = na.rm)/qnorm(0.75)
        if(any(SD == 0)){
            warning("'mad(x) = 0' for some samples => cannot compute a valid initial estimate.
                  To avoid division by zero 'mad0' is used. You could also specify 
                  a valid scale estimate via 'initial.est'. Note that you have to provide
                  a location and scale estimate.")
            SD[SD == 0] <- mad0
        }
        mean.sd <- cbind(mean = MEAN, sd = SD)
    }else{
        if(nrow(initial.est) != nrow(x) || ncol(initial.est) != 2)
            stop("'initial.est' has wrong dimension")
        MEAN <- initial.est[,1]
        SD <- initial.est[,2]
        if(any(SD <= 0))
            stop("initial estimate for scale <= 0 which is no valid scale estimate")
        mean.sd <- initial.est
    }

    if(!is.null(eps)){
        if(r > 10){
            b <- SD*1.618128043
            const <- 1.263094656
            A2 <- b^2*(1+r^2)/(1+const)
            A1 <- const*A2
            a <- -0.6277527697*A2/SD
            mse <- A1 + A2
        }else{
            A1 <- SD^2*.getA1.norm(r)
            A2 <- SD^2*.getA2.norm(r)
            a <- SD*.geta.norm(r)
            b <- SD*.getb.norm(r)
            mse <- A1 + A2
        }
        rmxEst <- .kstep.norm.matrix(x = x, initial.est = cbind(MEAN, SD), 
                                         A1 = A1, A2 = A2, a = a, b = b, k = k,
                                         na.rm = na.rm)
        colnames(rmxEst) <- c("mean", "sd")
        if(fsCor){
            Info.matrix <- matrix(c("rowRmx", 
                                    paste("fs-corrected estimate for 'eps' =", 
                                          round(eps, 3))),
                                  ncol = 2, dimnames = list(NULL, c("method", "message")))
        }else{
            Info.matrix <- matrix(c("rowRmx", 
                                    paste("asymptotic estimate for 'eps' =", 
                                          round(eps, 3))),
                                  ncol = 2, dimnames = list(NULL, c("method", "message")))
        }

        if(computeSE){
            asVar <- cbind(rmxEst[,2]^2 * .getAsVar.norm.approx(r)[1,1],
                           rmxEst[,2]^2 * .getAsVar.norm.approx(r)[2,2])
            asSE <- sqrt(asVar/ncol(x))
            colnames(asSE) <- c("SE.mean", "SE.sd")
        }else{
            asSE <- NA
        }
        
        RMX <- list(model = "norm", modelName = "normal location and scale",
                    rmxEst = rmxEst, asSE = asSE, radius = r, 
                    Infos = Info.matrix)
        class(RMX) <- "RMXlist"
    }else{
        if(r > 10){
            b <- SD*1.618128043
            const <- 1.263094656
            A2 <- b^2*(1+r^2)/(1+const)
            A1 <- const*A2
            a <- -0.6277527697*A2/SD
            mse <- A1 + A2
        }else{
            A1 <- SD^2*.getA1.norm(r)
            A2 <- SD^2*.getA2.norm(r)
            a <- SD*.geta.norm(r)
            b <- SD*.getb.norm(r)
            mse <- A1 + A2
        }
        if(rlo == 0){
            ineff <- (A1 + A2 - b^2*r.as^2)/(1.5*SD^2)
        }else{
            if(rlo > 10){
                ineff <- 1
            }else{
                A1lo <- SD^2*.getA1.norm(rlo)
                A2lo <- SD^2*.getA2.norm(rlo)
                ineff <- (A1 + A2 - b^2*(r.as^2 - rlo^2))/(A1lo + A2lo)
            }
        }
        rmxEst <- .kstep.norm.matrix(x = x, initial.est = cbind(MEAN, SD), 
                                     A1 = A1, A2 = A2, a = a, b = b, k = k,
                                     na.rm = na.rm)
        colnames(rmxEst) <- c("mean", "sd")
        if(fsCor){
            Info.matrix <- matrix(c(rep("rowRmx", 2), 
                                    paste("fs-corrected rmx estimate for eps in [", 
                                          round(eps.lower, 3), ", ", round(eps.upper, 3), "]", sep = ""),
                                    paste("least favorable (uncorrected) contamination: ", round(r.as/sqrtn, 3), sep = "")),
                                  ncol = 2, dimnames = list(NULL, c("method", "message")))
        }else{
            Info.matrix <- matrix(c(rep("rowRmx", 2), 
                                    paste("rmx estimate for eps in [", 
                                          round(eps.lower, 3), ", ", round(eps.upper, 3), "]", sep = ""),
                                    paste("least favorable contamination: ", round(r/sqrtn, 3), sep = "")),
                                  ncol = 2, dimnames = list(NULL, c("method", "message")))
        }
        if(computeSE){
            asVar <- cbind(rmxEst[,2]^2 * .getAsVar.norm.approx(r)[1,1],
                           rmxEst[,2]^2 * .getAsVar.norm.approx(r)[2,2])
            asSE <- sqrt(asVar/ncol(x))
            colnames(asSE) <- c("SE.mean", "SE.sd")
        }else{
            asSE <- NA
        }
        
        RMX <- list(model = "norm", modelName = "normal location and scale", 
                    rmxEst = rmxEst, asSE = asSE, radius = r, 
                    Infos = Info.matrix)
        class(RMX) <- "RMXlist"
    }
    RMX
}

###############################################################################
## computation of k-step construction in case a matrix x
###############################################################################
.onestep.norm.matrix <- function(x, initial.est, A1, A2, a, b, na.rm){
    MEAN <- initial.est[,1]
    SD <- initial.est[,2]
    u <- A1*(x-MEAN)/SD^2
    v <- A2*(((x-MEAN)/SD)^2-1)/SD - a
    ind <- b/sqrt(u^2 + v^2) <= 1
    IC1 <- rowMeans(u*(ind*b/sqrt(u^2 + v^2) + !ind), na.rm = na.rm)
    IC2 <- rowMeans(v*(ind*b/sqrt(u^2 + v^2) + !ind), na.rm = na.rm)
    IC <- cbind(IC1, IC2)
    initial.est + IC
}
.kstep.norm.matrix <- function(x, initial.est, A1, A2, a, b, k, na.rm){
    est <- .onestep.norm.matrix(x = x, initial.est = initial.est, 
                                A1 = A1, A2 = A2, a = a, b = b, na.rm = na.rm)
    if(k > 1){
        for(i in 2:k){
            A1 <- est[,2]^2*A1/initial.est[,2]^2
            A2 <- est[,2]^2*A2/initial.est[,2]^2
            a <- est[,2]*a/initial.est[,2]
            b <- est[,2]*b/initial.est[,2]
            initial.est <- est
            est <- .onestep.norm.matrix(x = x, initial.est = est,
                                        A1 = A1, A2 = A2, a = a, 
                                        b = b, na.rm = na.rm)
        }
    }
#    A1 <- est[,2]^2*A1/initial.est[,2]^2
#    A2 <- est[,2]^2*A2/initial.est[,2]^2
#    a <- est[,2]*a/initial.est[,2]
#    b <- est[,2]*b/initial.est[,2]
    
#    list(est = est, A1 = A1, A2 = A2, a = a, b = b)
    est
}


