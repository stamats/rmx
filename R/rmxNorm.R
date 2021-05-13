###############################################################################
## RMX estimator for normal location and scale
###############################################################################
rmx.norm <- function(x, eps.lower, eps.upper, eps, initial.est, k, fsCor){
    mad0 <- 1e-4
    if(!is.null(eps)){
        r <- sqrt(length(x))*eps
        if(fsCor){ 
            r.as <- r
            r <- fsRadius(r = r, n = length(x), model = "norm")
        }
    }else{
        sqrtn <- sqrt(length(x))
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
            r <- fsRadius(r = r, n = length(x), model = "norm")
        }
    }
    if(length(x) <= 2){
        warning("Sample size <= 2! => Median and MAD are used for estimation.")
        rmxEst <- c(median(x), mad(x))
        names(rmxEst) <- c("mean", "SD")
        Info.matrix <- matrix(c("rmx", "median and MAD"), ncol = 2, 
                              dimnames = list(NULL, c("method", "message")))
        AM <- rmxEst[1]
        SD <- rmxEst[2]
        b1 <- SD*sqrt(pi/2)
        A1 <- 1
        a1 <- 0
        
        b2 <- SD/(4*qnorm(0.75)*dnorm(qnorm(0.75)))
        A2 <- 1
        a2 <- (qnorm(0.75)^2 - 1)/SD
        A <- c(A1, A2)
        a <- c(a1, a2)
        b <- sqrt(b1^2 + b2^2)
        
        mse <- b1^2*(1+r^2) + b2^2*(1+r^2)
        names(mse) <- NULL
        bias <- r^2*(b1^2 + b2^2)
        names(bias) <- NULL
        V1 <- b1^2
        V2 <- b2^2
        asVar <- diag(c(V1, V2))
        rownames(asVar) <- c("mean", "SD")
        colnames(asVar) <- c("mean", "SD")
        
        param <- c(AM, SD)
        names(param) <- c("mean", "SD")
        IFun <- function(x){}
        body(IFun) <- substitute({ 
            z <- (x-AM)/sigma
            y1 <- b1*sign(z)
            y2 <- b2*sign((z^2 - 1)/sigma - a2)
            Y <- cbind(y1, y2)
            colnames(Y) <- c("IF for mean estimate", "IF for SD estimate")  
            Y
        }, list(AM = AM, sigma = SD, a2 = a2, b1 = b1, b2 = b2))
        range <- function(alpha, n = 501){} 
        body(range) <- substitute({
            rg <- qnorm(c(alpha/2, 1-alpha/2), mean = AM, sd = sigma)
            seq(from = rg[1], to = rg[2], length.out = n)
        }, list(AM = mean, sigma = sd))
        
        IF <- list(model = "norm", modelName = "normal location and scale", 
                   parameter = param, A = A, a = a, b = b, IFun = IFun,
                   range = range, asMSE = mse, asVar = asVar, asBias = bias,
                   radius = r)
        class(IF) <- "optIF"
        return(IF)
        
        RMX <- list(rmxEst = rmxEst, rmxIF = IF, IFinitial.est = NULL, 
                    Infos = Info.matrix)
        class(RMX) <- "rmx"
        return(RMX)
    }
    if(!is.null(eps)){
        if(eps == 0){
            n <- length(x)
            rmxEst <- c(mean(x), sqrt((n-1)/n)*sd(x))
            names(rmxEst) <- c("mean", "SD")
            Info.matrix <- matrix(c("rmx", "mean and SD"), ncol = 2, 
                                  dimnames = list(NULL, c("method", "message")))
            mean <- rmxEst[1]
            sd <- rmxEst[2]
            asVar <- sd^2*diag(c(1, 0.5))
            rownames(asVar) <- c("mean", "SD")
            colnames(asVar) <- c("mean", "SD")
            A <- sd^2*diag(c(1, 0.5))
            a <- c(0, 0)
            b <- Inf
            mse <- sd^2*(1 + 0.5)
            names(mse) <- NULL
            bias <- sqrt(mse - sum(diag(asVar)))
            names(bias) <- NULL
            param <- c(mean, sd)
            names(param) <- c("mean", "SD")
            IFun <- function(x){}
            body(IFun) <- substitute({ 
                z <- (x-AM)/sigma
                res <- sigma*cbind(z, 0.5*(z^2 - 1)) 
                colnames(res) <- c("IF for mean estimate", "IF for SD estimate") 
                res 
            }, list(AM = mean, sigma = sd))
            range <- function(alpha, n = 501){} 
            body(range) <- substitute({
                rg <- qnorm(c(alpha/2, 1-alpha/2), mean = AM, sd = sigma)
                seq(from = rg[1], to = rg[2], length.out = n)
            }, list(AM = mean, sigma = sd))
            
            IF <- list(model = "norm", modelName = "normal location and scale", 
                       parameter = param, A = A, a = a, b = b, IFun = IFun,
                       range = range, asMSE = mse, asVar = asVar, asBias = bias,
                       radius = r)
            class(IF) <- "optIF"

            RMX <- list(rmxEst = rmxEst, rmxIF = IF, initial.est = NULL, 
                        Infos = Info.matrix)
            class(RMX) <- "rmx"
            return(RMX)
        }
    }
    if(is.null(initial.est)){
        MEAN <- median(x, na.rm = TRUE)
        SD <- mad(x, center = MEAN, na.rm = TRUE)
        if(SD == 0){
            warning("'mad(x) = 0' => cannot compute a valid initial estimate. 
                     To avoid division by zero 'mad0' is used. You could also 
                     specify a valid initial estimate for scale via 'initial.est'. 
                     Note that you have to provide an initial estimate for mean and SD.")
            SD <- mad0
        }
    }else{
        stopifnot(is.numeric(initial.est))
        if(length(initial.est) != 2)
            stop("'initial.est' needs to be a numeric vector of length 2 or NULL")
        MEAN <- initial.est[1]
        SD <- initial.est[2]
        if(SD <= 0)
            stop("initial estimate for SD <= 0 which is no valid")
    }
    mean.sd <- c(MEAN, SD)
    names(mean.sd) <- c("mean","SD")
    if(!is.null(eps)){
        if(r > 10){
            b <- SD*1.618128043
            const <- 1.263094656
            A2 <- b^2*(1+r^2)/(1+const)
            A1 <- const*A2
            a <- c(0, -0.6277527697*A2/SD)
            mse <- A1 + A2
        }else{
            A1 <- SD^2*.getA1.norm(r)
            A2 <- SD^2*.getA2.norm(r)
            a <- SD*c(0, .geta.norm(r))
            b <- SD*.getb.norm(r)
            mse <- A1 + A2
        }
        rmxEst.all <- .kstep.norm(x = x, initial.est = c(MEAN, SD), 
                              A1 = A1, A2 = A2, a = a, b = b, k = k)
        rmxEst <- rmxEst.all$est
        names(rmxEst) <- c("mean", "SD")
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
        A1 <- rmxEst.all$A1 
        A2 <- rmxEst.all$A2
        A <- diag(c(A1, A2))
        a <- rmxEst.all$a 
        b <- rmxEst.all$b
        a1 <- A1/rmxEst[2]^2
        a3 <- A2/rmxEst[2]^2
        a2 <- a[2]/rmxEst[2]/a3 + 1
        
        asVar <- rmxEst.all$asvar
        rownames(asVar) <- c("mean", "SD")
        colnames(asVar) <- c("mean", "SD")
        mse <- rmxEst[2]^2*(a1 + a3)
        names(mse) <- NULL
        bias <- sqrt(mse - sum(diag(asVar)))
        names(bias) <- NULL
        param <- rmxEst
        IFun <- function(x){}
        body(IFun) <- substitute({ 
            z <- (x-AM)/sigma
            hvkt <- sqrt(a1^2*z^2/sigma^2 + (a3*(z^2-1)/sigma - a2)^2)
            ind1 <- (hvkt < b)
            w <- ind1 + (1-ind1)*b/hvkt 
            Y <- cbind(a1*z/sigma, a3*(z^2-1)/sigma-a2)
            res <- Y*w
            colnames(res) <- c("IF for mean estimate", "IF for SD estimate")  
            res 
        }, list(AM = rmxEst[1], sigma = rmxEst[2], 
                a1 = A[1,1], a2 = a[2], a3 = A[2,2], b = b))
        range <- function(alpha, n = 501){} 
        body(range) <- substitute({
            rg <- qnorm(c(alpha/2, 1-alpha/2), mean = AM, sd = sigma)
            seq(from = rg[1], to = rg[2], length.out = n)
        }, list(AM = rmxEst[1], sigma = rmxEst[2]))
        
        IF <- list(model = "norm", modelName = "normal location and scale", 
                   parameter = param, A = A, a = a, b = b, IFun = IFun,
                   range = range, asMSE = mse, asVar = asVar, asBias = bias,
                   radius = r)
        class(IF) <- "optIF"
    }else{
        if(r > 10){
            b <- SD*1.618128043
            const <- 1.263094656
            A2 <- b^2*(1+r^2)/(1+const)
            A1 <- const*A2
            a <- c(0, -0.6277527697*A2/SD)
            mse <- A1 + A2
        }else{
            A1 <- SD^2*.getA1.norm(r)
            A2 <- SD^2*.getA2.norm(r)
            a <- SD*c(0, .geta.norm(r))
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
        rmxEst.all <- .kstep.norm(x = x, initial.est = c(MEAN, SD), 
                                  A1 = A1, A2 = A2, a = a, b = b, k = k)
        rmxEst <- rmxEst.all$est
        names(rmxEst) <- c("mean", "SD")
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
        A1 <- rmxEst.all$A1 
        A2 <- rmxEst.all$A2 
        A <- diag(c(A1, A2))
        a <- rmxEst.all$a 
        b <- rmxEst.all$b
        a1 <- A1/rmxEst[2]^2
        a3 <- A2/rmxEst[2]^2
        a2 <- a[2]/rmxEst[2]/a3 + 1
        
        asVar <- rmxEst.all$asvar
        rownames(asVar) <- c("mean", "SD")
        colnames(asVar) <- c("mean", "SD")
        mse <- rmxEst[2]^2*(a1 + a3)
        names(mse) <- NULL
        bias <- sqrt(mse - sum(diag(asVar)))
        names(bias) <- NULL
        param <- rmxEst
        IFun <- function(x){}
        body(IFun) <- substitute({ 
            z <- (x-AM)/sigma
            hvkt <- sqrt(a1^2*z^2/sigma^2 + (a3*(z^2-1)/sigma - a2)^2)
            ind1 <- (hvkt < b)
            w <- ind1 + (1-ind1)*b/hvkt 
            Y <- cbind(a1*z/sigma, a3*(z^2-1)/sigma-a2)
            res <- Y*w
            colnames(res) <- c("IF for mean estimate", "IF for SD estimate")  
            res 
        }, list(AM = rmxEst[1], sigma = rmxEst[2], 
                a1 = A[1,1], a2 = a[2], a3 = A[2,2], b = b))
        range <- function(alpha, n = 501){} 
        body(range) <- substitute({
            rg <- qnorm(c(alpha/2, 1-alpha/2), mean = AM, sd = sigma)
            seq(from = rg[1], to = rg[2], length.out = n)
        }, list(AM = rmxEst[1], sigma = rmxEst[2]))
        
        IF <- list(model = "norm", modelName = "normal location and scale", 
                   parameter = param, A = A, a = a, b = b, IFun = IFun,
                   range = range, asMSE = mse, asVar = asVar, asBias = bias,
                   radius = r)
        class(IF) <- "optIF"
    }
    RMX <- list(rmxEst = rmxEst, rmxIF = IF, initial.est = c(MEAN, SD), 
                Infos = Info.matrix)
    class(RMX) <- "rmx"
    RMX
}

###############################################################################
## computation of radius-minimax IC
## using predefined functions included in "sysdata.rda"
###############################################################################
.getInterval.norm <- function(r, rlo, rup){
    if(r > 10){
        b <- 1.618128043
        const <- 1.263094656
        A2 <- b^2*(1+r^2)/(1+const)
        A1 <- const*A2
    }else{
        A1 <- .getA1.norm(r)
        A2 <- .getA2.norm(r)
        b <- .getb.norm(r)
    }
    
    if(rlo == 0){
        efflo <- (A1 + A2 - b^2*r^2)/1.5
    }else{
        A1lo <- .getA1.norm(rlo)
        A2lo <- .getA2.norm(rlo)
        efflo <- (A1 + A2 - b^2*(r^2 - rlo^2))/(A1lo + A2lo)
    }
    
    if(rup > 10){
        bup <- 1.618128043
        const.up <- 1.263094656
        A2up <- bup^2*(1+rup^2)/(1+const.up)
        A1up <- const.up*A2up
        effup <- (A1 + A2 - b^2*(r^2 - rup^2))/(A1up + A2up)
    }else{
        A1up <- .getA1.norm(rup)
        A2up <- .getA2.norm(rup)
        effup <- (A1 + A2 - b^2*(r^2 - rup^2))/(A1up + A2up)
    }
    
    return(effup-efflo)
}

###############################################################################
## computation of k-step construction
###############################################################################
.onestep.norm <- function(x, initial.est, A1, A2, a, b){
    MEAN <- initial.est[1]
    SD <- initial.est[2]
    u <- A1*(x-MEAN)/SD^2
    v <- A2*(((x-MEAN)/SD)^2-1)/SD - a[2]
    w <- pmin(1, b/sqrt(u^2 + v^2))
    IC <- c(mean(u*w, na.rm = TRUE), mean(v*w, na.rm = TRUE))
    return(initial.est + IC)
}
.kstep.norm <- function(x, initial.est, A1, A2, a, b, k){
    est <- .onestep.norm(x = x, initial.est = initial.est, 
                         A1 = A1, A2 = A2, a = a, b = b)
    if(k > 1){
        for(i in 2:k){
            A1 <- est[2]^2*A1/initial.est[2]^2
            A2 <- est[2]^2*A2/initial.est[2]^2
            a <- est[2]*a/initial.est[2]
            b <- est[2]*b/initial.est[2]
            initial.est <- est
            est <- .onestep.norm(x = x, initial.est = est,
                                 A1 = A1, A2 = A2, a = a, b = b)
        }
    }
    A1 <- est[2]^2*A1/initial.est[2]^2
    A2 <- est[2]^2*A2/initial.est[2]^2
    a <- est[2]*a/initial.est[2]
    b <- est[2]*b/initial.est[2]
    a1 <- A1/est[2]^2
    a3 <- A2/est[2]^2
    a2 <- a[2]/est[2]/a3 + 1
    asVar <- est[2]^2*.getAsVar.norm(b = b/est[2], a1 = a1, a2 = a2, a3 = a3)
    
    return(list(est = est, A1 = A1, A2 = A2, a = a, b = b, asvar = asVar))
}