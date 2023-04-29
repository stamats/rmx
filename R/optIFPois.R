###############################################################################
## optimal IF for normal location and scale
###############################################################################
optIF.pois <- function(radius, lambda = 1, aUp = 100*lambda, cUp = 1e4, 
                        delta = 1e-9, check = FALSE){
    if(radius == 0){ # IF of ML estimator
        IF <- .getMLIF.pois(lambda = lambda)
        IF$radius <- 0
        return(IF)
    }
    lcr <- .lcr.pois(lambda = lambda)
    if(radius >= lcr){ # IF of minimum bias estimator
        IF <- .getMBIF.pois(lambda = lambda)
        IF$asBias <- radius*IF$b
        IF$asMSE <- IF$asVar + IF$asBias^2
        IF$radius <- radius
        
        return(IF)
    }
    
    ## compute Lagrange multipliers
    LM <- .getLM.pois(r0 = radius, lambda = lambda, aUp = aUp, 
                      cUp = cUp, delta = delta)
    A <- LM$A
    z <- LM$z
    a <- A*z
    c0 <- LM$c0
    b <- A*c0
    
    if(check){
        supp <- seq(from = 0, to = qpois(1-1e-15, lambda = lambda), by = 1)
        WS <- dpois(supp, lambda = lambda)
        Y0 <- supp/lambda - 1
        M <- pmin(abs(Y0 - z), c0)
        ch1 <- sum(A*abs(Y0-z)*pmin(abs(Y0 - z), c0)*WS)
        
        ind2 <- abs(Y0-z) < c0
        res <- ind2*(Y0-z) + (1-ind2)*c0*sign(Y0-z)
        ch2 <- sum(res*WS)
        
        ch3 <- sum(pmax(abs(Y0-z)-c0, 0)*WS) - radius^2*c0

        message("Fisher consistency of eta:\t", ch1-1)
        message("centering of eta:\t", ch2)
        message("MSE equation:\t", ch3)
    }

    IF <- .getOptIF.pois(lambda = lambda, radius = radius, A = A, a = a, b = b)
    IF$radius <- radius
    IF
}

###############################################################################
## lower case radius
###############################################################################
.lcr.pois <- function(lambda){
    m0 <- qpois(0.5, lambda = lambda)
    wsm <- dpois(m0, lambda = lambda)
    
    if(wsm > 0){
        supp <- seq(from = 0, to = qpois(1-1e-15, lambda = lambda), by = 1)
        gam <- min(abs(supp[supp != m0] - m0))/lambda
        if(gam > 0){
            Int <- sum(abs(supp - m0)*dpois(supp, lambda = lambda))/lambda
            p1 <- ppois(m0-1, lambda = lambda)
            p2 <- ppois(m0, lambda = lambda, lower.tail = FALSE)
            beta <- (p1 - p2)/wsm
            
            if(p1 == 0 || p2 == 0){
                M <- (1 + abs(beta))/gam
            }else{
                gam1 <- min(supp[supp > m0] - m0)/lambda
                gam2 <- max(supp[supp < m0] - m0)/lambda
                M <- max((1-beta)/gam1, (1+beta)/(-gam2))
            }
            rad <- sqrt(max(M*Int - (1-wsm) - beta^2*wsm, 0))
        }else{
            rad <- Inf
        }
    }else{
        gam1 <- min(supp[supp > m0] - m0)/lambda
        gam2 <- max(supp[supp < m0] - m0)/lambda
        if(gam1 > 0 || gam2 < 0){
            Int <- sum(abs(supp - m0)*dpois(supp, lambda = lambda))/lambda
            M <- 2/(gam1-gam2)
            rad <- sqrt(max(M*Int - 1, 0))
        }else{
            rad <- Inf
        }
    }
    rad
}

##################################################################
## centering
##################################################################
.geta.pois <- function(a, c0, lambda){
    c1 <- ceiling(lambda*(a - c0))
    c2 <- floor(lambda*(a + c0))
    
    p11 <- ppois(c1-1, lambda=lambda)
    p12 <- ppois(c1-2, lambda=lambda)
    p21 <- ppois(c2, lambda=lambda)
    p22 <- ppois(c2-1, lambda=lambda)
    
    c0 + (a-c0)*p11 + (p22-p12)*(p22>=p12) - (a+c0)*p21
}
##################################################################
## clipping
##################################################################
.getc.pois <- function(c0, r0, lambda, aUp, delta){
    a <- uniroot(.geta.pois, lower = 0, upper = aUp, tol = delta, c0 = c0, 
                 lambda = lambda)$root 
    
    c1 <- floor(lambda*(a-c0))
    c2 <- ceiling(lambda*(a+c0))
    
    p11 <- ppois(c1,lambda=lambda)
    p12 <- ppois(c1-1,lambda=lambda)
    p21 <- 1-ppois(c2-1,lambda=lambda)
    p22 <- 1-ppois(c2-2,lambda=lambda)
    
    p22 - p12 - (a+c0)*p21 + (a-c0)*p11 - r0^2*c0
}
##################################################################
## Lagrange multipliers
##################################################################
.getLM.pois <- function(r0, lambda, aUp = 100*lambda, cUp = 1e4, delta = 1e-9){
    c0 <- uniroot(.getc.pois, lower = 1e-10, upper = cUp, 
                  tol = delta, r0 = r0, lambda = lambda, 
                  aUp = aUp, delta = delta)$root
    a <- uniroot(.geta.pois, lower = 0, upper = aUp, tol = delta, c0 = c0, 
                 lambda = lambda)$root
    
    c1 <- ceiling(lambda*(a - c0))
    c2 <- floor(lambda*(a + c0))
    
    p11 <- ppois(c1-1,lambda=lambda)
    p12 <- ppois(c1-2,lambda=lambda)
    p13 <- ppois(c1-3,lambda=lambda)
    p21 <- ppois(c2,lambda=lambda)
    p22 <- ppois(c2-1,lambda=lambda)
    p23 <- ppois(c2-2,lambda=lambda)
    
    aa <- c0*a*p11 - c0*p12 + (p23 - p13 + (p22-p12)*(1/lambda-2*a) + 
                                   a^2*(p21-p11))*(c2>=c1)
    aa <- aa + c0*(1-p22) - c0*a*(1-p21)
    
    A <- 1/aa
    z <- a - 1
    
    list(A = A, z = z, c0 = c0)    
}
##################################################################
## IF of ML estimator
##################################################################
.getMLIF.pois <- function(lambda){
    asVar <- lambda
    A <- asVar
    names(asVar) <- "lambda"
    a <- 0
    b <- Inf
    mse <- asVar
    names(mse) <- NULL
    bias <- sqrt(mse - asVar)
    names(bias) <- NULL
    param <- lambda
    names(param) <- "lambda"
    IFun <- function(x){}
    body(IFun) <- substitute({ 
        res <- matrix(x - lambda, ncol = 1)
        colnames(res) <- c("IF for lambda") 
        res 
    }, list(lambda = lambda))
    range <- function(alpha){} 
    body(range) <- substitute({
        rg <- qpois(c(alpha/2, 1-alpha/2), lambda = lambda)
        seq(from = rg[1], to = rg[2], by = 1)
    }, list(lambda = lambda))
    
    IF <- list(model = "pois", modelName = "Poisson mean", 
               parameter = param, A = A, a = a, b = b, IFun = IFun,
               range = range, asMSE = mse, asVar = asVar, asBias = bias)
    class(IF) <- "optIF"
    IF
}
##################################################################
## IF of minimum bias estimator
##################################################################
.getMBIF.pois <- function(lambda){
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
    mse <- Inf
    bias <- Inf
    asVar <- b^2
    names(asVar) <- "lambda"
    
    param <- lambda
    names(param) <- "lambda"
    IFun <- function(x){}
    body(IFun) <- substitute({
        res <- matrix(bmin*((x > M) - (x < M) + beta*(x == M)), ncol = 1)
        colnames(res) <- "IF for lambda"
        res 
    }, list(bmin = b, M = m0, beta = beta))
    range <- function(alpha){} 
    body(range) <- substitute({
        rg <- qpois(c(alpha/2, 1-alpha/2), lambda = lambda)
        seq(from = rg[1], to = rg[2], by = 1)
    }, list(lambda = lambda))
    
    IF <- list(model = "pois", modelName = "Poisson mean", 
               parameter = param, A = A, a = a, b = b, lowerCase = beta, 
               IFun = IFun, range = range, asMSE = mse, asVar = asVar, 
               asBias = bias)
    class(IF) <- "optIF"
    IF
}
##################################################################
## IF of optimally robust estimator
##################################################################
.getOptIF.pois <- function(lambda, radius, A, a, b){
    asVar <- A - radius^2*b^2
    names(asVar) <- "lambda"
    mse <- A
    names(mse) <- NULL
    bias <- radius*b
    names(bias) <- NULL
    param <- lambda
    names(param) <- c("lambda")
    IFun <- function(x){}
    body(IFun) <- substitute({ 
        Y0 <- x/lambda - 1
        Y <- A*Y0 - a
        ind1 <- abs(Y) < b
        res <- matrix(ind1*Y + (1-ind1)*b*sign(Y), ncol = 1)
        colnames(res) <- "IF for lambda"
        res 
    }, list(lambda = lambda, A = A, a = a, b = b))
    range <- function(alpha){} 
    body(range) <- substitute({
        rg <- qpois(c(alpha/2, 1-alpha/2), lambda = lambda)
        seq(from = rg[1], to = rg[2], by = 1)
    }, list(lambda = lambda))
    
    IF <- list(model = "pois", modelName = "Poisson mean", 
               parameter = param, A = A, a = a, b = b, IFun = IFun,
               range = range, asMSE = mse, asVar = asVar, asBias = bias)
    class(IF) <- "optIF"
    IF
}
