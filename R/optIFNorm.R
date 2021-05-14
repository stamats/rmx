###############################################################################
## optimal IF for normal location and scale
###############################################################################
optIF.norm <- function(radius, mean = 0, sd = 1, A.loc.start = 1, a.sc.start = 0, 
                       A.sc.start = 0.5, bUp = 1000, delta = 1e-6, itmax = 100L, 
                       check = FALSE){
    if(radius == 0){
        asVar <- sd^2*diag(c(1, 0.5))
        rownames(asVar) <- c("mean", "sd")
        colnames(asVar) <- c("mean", "sd")
        A <- sd^2*diag(c(1, 0.5))
        a <- c(0, 0)
        b <- Inf
        mse <- sd^2*(1 + 0.5)
        names(mse) <- NULL
        bias <- sqrt(mse - sum(diag(asVar)))
        names(bias) <- NULL
        param <- c(mean, sd)
        names(param) <- c("mean", "sd")
        IFun <- function(x){}
        body(IFun) <- substitute({ 
                z <- (x-AM)/sigma
                res <- sigma*cbind(z, 0.5*(z^2 - 1)) 
                colnames(res) <- c("IF for mean", "IF for sd") 
                res 
            }, list(AM = mean, sigma = sd))
        range <- function(alpha, n = 501){} 
        body(range) <- substitute({
                rg <- qnorm(c(alpha/2, 1-alpha/2), mean = AM, sd = sigma)
                seq(from = rg[1], to = rg[2], length.out = n)
            }, list(AM = mean, sigma = sd))
        
        IF <- list(model = "norm", modelName = "normal location and scale", 
                   parameter = param, A = A, a = a, b = b, IFun = IFun,
                   range = range, asMSE = mse, asVar = asVar, asBias = bias)
        class(IF) <- "optIF"
        return(IF)
    }
    if(radius == Inf){
        b <- sd*1.618128043
        const <- 1.263094656
        a1 <- 1
        a3 <- 1/const
        A <- diag(c(a1, a3))
        a <- c(0, -0.6277527697/const/sd)
        mse <- Inf
        bias <- Inf
        ## fun1 <- function(x) {x^2 / (x^2 + (1/1.263094656*(x^2-1) + 0.6277527697/1.263094656)^2)*dnorm(x) }
        ## print(2*integrate(fun1, lower = 0, upper = Inf, rel.tol = 1e-14)$value, digits = 10)
        V1 <- 0.6347635562
        ## fun2 <- function(x) {(1/1.263094656*(x^2-1) + 0.6277527697/1.263094656)^2 / (x^2 + (1/1.263094656*(x^2-1) + 0.6277527697/1.263094656)^2)*dnorm(x) }
        ## print(2*integrate(fun2, lower = 0, upper = Inf, rel.tol = 1e-14)$value, digits = 10)
        V2 <- 0.3652364438
        asVar <- b^2*diag(c(V1, V2))
        rownames(asVar) <- c("mean", "sd")
        colnames(asVar) <- c("mean", "sd")
        
        param <- c(mean, sd)
        names(param) <- c("mean", "sd")
        IFun <- function(x){}
        body(IFun) <- substitute({ 
                z <- (x-AM)/sigma
                w <- 1/sqrt(z^2 + (a3*(z^2-1) - a2*sigma)^2)
                Y <- b*cbind(z, a3*(z^2-1)-a2*sigma)
                res <- Y*w
                colnames(res) <- c("IF for mean", "IF for sd") 
                res 
            }, list(AM = mean, sigma = sd, a2 = a[2], a3 = a3, b = b))
        range <- function(alpha, n = 501){} 
        body(range) <- substitute({
            rg <- qnorm(c(alpha/2, 1-alpha/2), mean = AM, sd = sigma)
            seq(from = rg[1], to = rg[2], length.out = n)
        }, list(AM = mean, sigma = sd))
        
        IF <- list(model = "norm", modelName = "normal location and scale", 
                   parameter = param, A = A, a = a, b = b, IFun = IFun,
                   range = range, asMSE = mse, asVar = asVar, asBias = bias)
        class(IF) <- "optIF"
        return(IF)
    }
    r <- radius
    a1 <- A.loc.start 
    a2 <- 1+a.sc.start 
    a3 <- A.sc.start
    b <- uniroot(.getr.norm, lower = 1e-4, upper = bUp, 
            tol = .Machine$double.eps^0.5, r = r, a1 = a1, a2 = a2, 
            a3 = a3)$root

    iter <- 0
    repeat{
        iter <- iter + 1
        if(iter > itmax){
            stop("Algorithm did not converge!\n", 
                 "=> increase itmax or try different starting values",
                 "'A.loc.start', 'a.sc.start' and 'A.sc.start'\n")
        }
        a1.old <- a1; a2.old <- a2; a3.old <- a3; b.old <- b

        a1a2a3 <- .geta1a2a3.norm(b = b, a1 = a1, a2 = a2, a3 = a3)
        a1 <- a1a2a3$a1
        a2 <- a1a2a3$a2
        a3 <- a1a2a3$a3

        b <- uniroot(.getr.norm, lower = 1e-4, upper = bUp, 
                tol = .Machine$double.eps^0.5, r = r, a1 = a1, a2 = a2, 
                a3 = a3)$root
        if(max(abs(a1.old-a1), abs(a2.old-a2), abs(a3.old-a3), abs(b.old-b))<delta)
            break
    }

    if(check){
        integrand1 <- function(x, b, a1, a2, a3){
            x^2*.getw.norm(x, b, a1, a2, a3)*dnorm(x)
        }
        Int1 <- 2*integrate(integrand1, lower = 0, upper = Inf, 
                        rel.tol = .Machine$double.eps^0.5, b = b, a1 = a1, 
                        a2 = a2, a3 = a3)$value
        ch1 <- a1*Int1

        integrand2 <- function(x, b, a1, a2, a3){
            (x^2 - a2)^2*.getw.norm(x, b, a1, a2, a3)*dnorm(x)
        }
        Int2 <- 2*integrate(integrand2, lower = 0, upper = Inf, 
                        rel.tol = .Machine$double.eps^0.5, b = b, a1 = a1, 
                        a2 = a2, a3 = a3)$value
        ch2 <- a3*Int2

        integrand3 <- function(x, b, a1, a2, a3){
            (x^2 - a2)*.getw.norm(x, b, a1, a2, a3)*dnorm(x)
        }
        Int3 <- 2*integrate(integrand3, lower=0, upper=Inf, 
                        rel.tol = .Machine$double.eps^0.5, b = b, a1 = a1, 
                        a2 = a2, a3 = a3)$value
        ch3 <- a3*Int3

        ch4 <- .getr.norm(b = b, r = r, a1 = a1, a2 = a2, a3 = a3)

        message("Fisher consistency of eta.loc:\t", ch1-1)
        message("Fisher consistency of eta.sc:\t", ch2-1)
        message("centering of eta.sc:\t", ch3)
        message("MSE equation:\t", ch4)
    }

    asVar <- sd^2*.getAsVar.norm(b = b, a1 = a1, a2 = a2, a3 = a3)
    rownames(asVar) <- c("mean", "sd")
    colnames(asVar) <- c("mean", "sd")
    A <- sd^2*diag(c(a1, a3))
    a <- sd*c(0, a3*(a2-1))
    b <- sd*b
    mse <- sd^2*(a1 + a3)
    names(mse) <- NULL
    bias <- sqrt(mse - sum(diag(asVar)))
    names(bias) <- NULL
    param <- c(mean, sd)
    names(param) <- c("mean", "sd")
    IFun <- function(x){}
    body(IFun) <- substitute({ 
        z <- (x-AM)/sigma
        hvkt <- sqrt(a1^2*z^2/sigma^2 + (a3*(z^2-1)/sigma - a2)^2)
        ind1 <- (hvkt < b)
        w <- ind1 + (1-ind1)*b/hvkt 
        Y <- cbind(a1*z/sigma, a3*(z^2-1)/sigma-a2)
        res <- Y*w
        colnames(res) <- c("IF for mean", "IF for sd") 
        res 
        }, list(AM = mean, sigma = sd, 
                a1 = A[1,1], a2 = a[2], a3 = A[2,2], b = b))
    range <- function(alpha, n = 501){} 
    body(range) <- substitute({
        rg <- qnorm(c(alpha/2, 1-alpha/2), mean = AM, sd = sigma)
        seq(from = rg[1], to = rg[2], length.out = n)
    }, list(AM = mean, sigma = sd))
    
    IF <- list(model = "norm", modelName = "normal location and scale", 
               parameter = param, A = A, a = a, b = b, IFun = IFun,
               range = range, asMSE = mse, asVar = asVar, asBias = bias)
    class(IF) <- "optIF"
    IF
}

###############################################################################
## weight function
###############################################################################
.getw.norm <- function(x, b, a1, a2, a3){
    hvkt <- sqrt(a3^2*x^4 + (a1^2 - 2*a2*a3^2)*x^2 + a2^2*a3^2)
    ind1 <- (hvkt < b)
    
    return(ind1 + (1-ind1)*b/hvkt)
}

###############################################################################
## computation of r
###############################################################################
.getr.norm <- function(b, r, a1, a2, a3){
    integrandr <- function(x, b, a1, a2, a3){
        hvkt <- sqrt(a3^2*x^4 + (a1^2 - 2*a2*a3^2)*x^2 + a2^2*a3^2)/b - 1
        return((hvkt > 0)*hvkt*dnorm(x))
    }
    Int <- integrate(integrandr, lower = 0, upper = Inf, 
                     rel.tol = .Machine$double.eps^0.5, a1 = a1, a2 = a2, 
                     a3 = a3, b = b)$value
    
    return(r-sqrt(2*Int))
}


###############################################################################
## computation of a1, a2 and a3
###############################################################################
.geta1a2a3.norm <- function(b, a1, a2, a3){
    integrand1 <- function(x, b, a1, a2, a3){ 
        x^2*.getw.norm(x, b, a1, a2, a3)*dnorm(x)
    }
    Int1 <- 2*integrate(integrand1, lower = 0, upper = Inf, 
                        rel.tol = .Machine$double.eps^0.5, b = b, a1 = a1, 
                        a2 = a2, a3 = a3)$value
    a1 <- 1/Int1
    
    integrand2 <- function(x, b, a1, a2, a3){
        .getw.norm(x, b, a1, a2, a3)*dnorm(x)
    }
    Int2 <- 2*integrate(integrand2, lower = 0, upper = Inf, 
                        rel.tol = .Machine$double.eps^0.5, b = b, a1 = a1, 
                        a2 = a2, a3 = a3)$value
    a2 <- Int1/Int2
    
    integrand3 <- function(x, b, a1, a2, a3){
        (x^2 - a2)^2*.getw.norm(x, b, a1, a2, a3)*dnorm(x)
    }
    Int3 <- 2*integrate(integrand3, lower = 0, upper = Inf, 
                        rel.tol = .Machine$double.eps^0.5, b = b, a1 = a1, 
                        a2 = a2, a3 = a3)$value
    a3 <- 1/Int3
    
    return(list(a1=a1, a2=a2, a3=a3))
}

###############################################################################
## computation of asymptotic (co)variance
###############################################################################
.getAsVar.norm <- function(b, a1, a2, a3){
    integrand1 <- function(x, b, a1, a2, a3){ 
        x^2*.getw.norm(x, b, a1, a2, a3)^2*dnorm(x)
    }
    Int1 <- 2*integrate(integrand1, lower = 0, upper = Inf, 
                        rel.tol = .Machine$double.eps^0.5, b = b, a1 = a1, 
                        a2 = a2, a3 = a3)$value
    V1 <- a1^2*Int1
    
    integrand2 <- function(x, b, a1, a2, a3){
        (x^2 - a2)^2*.getw.norm(x, b, a1, a2, a3)^2*dnorm(x)
    }
    Int2 <- 2*integrate(integrand2, lower = 0, upper = Inf, 
                        rel.tol = .Machine$double.eps^0.5, b = b, a1 = a1, 
                        a2 = a2, a3 = a3)$value
    V2 <- a3^2*Int2
    
    return(diag(c(V1, V2)))
}
