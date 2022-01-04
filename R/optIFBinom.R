###############################################################################
## optimal IF for normal location and scale
###############################################################################
optIF.binom <- function(radius, size = 1, prob = 0.5, aUp = 100*size, cUp = 1e4, 
                        delta = 1e-9, check = FALSE){
    if(radius == 0){ # IF of ML estimator
        IF <- .getMLIF.binom(size = size, prob = prob)
        IF$radius <- 0
        return(IF)
    }
    lcr <- .lcr.binom(prob = prob, size = size)
    if(radius >= lcr){ # IF of minimum bias estimator
        IF <- .getMBIF.binom(size = size, prob = prob)
        IF$radius <- radius
        
        return(IF)
    }
    
    ## compute Lagrange multipliers
    LM <- .getLM.binom(r0 = radius, prob = prob, size = size, aUp = aUp, 
                       cUp = cUp, delta = delta)
    A <- LM$A
    z <- LM$z
    a <- A*z
    c0 <- LM$c0
    b <- A*c0
    
    if(check){
        supp <- seq(from = 0, to = size, by = 1)
        WS <- dbinom(supp, size = size, prob = prob)
        Y0 <- (supp - size*prob)/(prob*(1-prob))
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

    IF <- .getOptIF.binom(size = size, prob = prob, radius = radius, 
                          A = A, a = a, b = b)
    IF$radius <- radius
    IF
}

###############################################################################
## lower case radius
###############################################################################
.lcr.binom <- function(prob, size){
    m0 <- qbinom(0.5, size = size, prob = prob)
    wsm <- dbinom(m0, size = size, prob = prob)
    
    if(wsm > 0){
        supp <- seq(from = 0, to = size, by = 1)
        gam <- min(abs(supp[supp != m0] - m0))/(prob*(1-prob))
        if(gam > 0){
            Int <- sum(abs(supp - m0)*dbinom(supp, size = size, prob = prob))/(prob*(1-prob))
            p1 <- pbinom(m0-1, size = size, prob = prob)
            p2 <- pbinom(m0, size = size, prob = prob, lower.tail = FALSE)
            beta <- (p1 - p2)/wsm
            
            if(p1 == 0 || p2 == 0){
                M <- (1 + abs(beta))/gam
            }else{
                gam1 <- min(supp[supp > m0] - m0)/(prob*(1-prob))
                gam2 <- max(supp[supp < m0] - m0)/(prob*(1-prob))
                M <- max((1-beta)/gam1, (1-beta)/(-gam2))
            }
            rad <- sqrt(max(M*Int - (1-wsm) - beta^2*wsm, 0))
        }else{
            rad <- Inf
        }
    }else{
        gam1 <- min(supp[supp > m0] - m0)/(prob*(1-prob))
        gam2 <- max(supp[supp < m0] - m0)/(prob*(1-prob))
        if(gam1 > 0 || gam2 < 0){
            Int <- sum(abs(supp - m0)*dbinom(supp, size = size, prob = prob))/(prob*(1-prob))
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
.geta.binom <- function(a, c0, prob, size){
    c1 <- ceiling(prob*(1-prob)*(a - c0))
    c2 <- floor(prob*(1-prob)*(a + c0))
    
    if(size > 1){
        p11 <- pbinom(c1-1, size=size, prob=prob)
        p12 <- pbinom(c1-2, size=size-1, prob=prob)
        p21 <- pbinom(c2, size=size, prob=prob)
        p22 <- pbinom(c2-1, size=size-1, prob=prob)
        
        fa <- c0 + (a-c0)*p11 + size/(1-prob)*(p22-p12)*(c2>=c1) - (a+c0)*p21
    }
    else{
        fa <- - a*min(1,c0/a)*(1-prob) + (1/(prob*(1-prob))-a)*min(1,c0/abs(1/(prob*(1-prob))-a))*prob
    }
    fa    
}
##################################################################
## clipping
##################################################################
.getc.binom <- function(c0, r0, prob, size, aUp, delta){
    a <- uniroot(.geta.binom, lower = 0, upper = aUp, 
                 tol = delta, c0 = c0, prob = prob, size = size)$root 

    c1 <- floor(prob*(1-prob)*(a - c0))
    c2 <- ceiling(prob*(1-prob)*(a + c0))
    
    if(size > 1){
        p11 <- pbinom(c1, size=size, prob=prob)
        p12 <- pbinom(c1-1, size=size-1, prob=prob)
        p21 <- pbinom(c2-1, size=size, prob=prob, lower.tail = FALSE)
        p22 <- pbinom(c2-2, size=size-1, prob=prob, lower.tail = FALSE)
        
        erg <- size/(1-prob)*(p22-p12) - (a+c0)*p21 + (a-c0)*p11 - r0^2*c0
    }
    else{
        erg <- max(0, a-c0)*(1-prob) + max(0, abs(1/(prob*(1-prob))-a)-c0)*prob - r0^2*c0
    }
    
    erg
}
##################################################################
## Lagrange multipliers
##################################################################
.getLM.binom <- function(r0, prob, size, aUp = 100*size, cUp = 1e4, delta = 1e-9){
    c0 <- uniroot(.getc.binom, lower = 1e-10, upper = cUp, 
                  tol = delta, r0 = r0, prob = prob, size = size, 
                  aUp = aUp, delta = delta)$root
    a <- uniroot(.geta.binom, lower = 0, upper = aUp, 
                 tol = delta, c0 = c0, prob = prob, size = size)$root
    
    c1 <- ceiling(prob*(1-prob)*(a - c0))
    c2 <- floor(prob*(1-prob)*(a + c0))
    
    if(size > 2){
        p11 <- pbinom(c1-2, size=size-1, prob=prob)
        p12 <- pbinom(c1-3, size=size-2, prob=prob)
        p21 <- pbinom(c2-1, size=size-1, prob=prob)
        p22 <- pbinom(c2-2, size=size-2, prob=prob)
        
        aa <- size/(1-prob)*((a-c0)*p11 - (a+c0)*p21 + c0) 
        aa <- aa + size/(1-prob)^2*((p21-p11)/prob + (size-1)*(p22-p12))*(c2>=c1)
    }
    else{
        if(size > 1){
            p11 <- pbinom(c1-2, size=size-1, prob=prob)
            p21 <- pbinom(c2-1, size=size-1, prob=prob)
            
            aa <- size/(1-prob)*((a-c0)*p11 - (a+c0)*p21 + c0) 
            aa <- aa + size/(1-prob)^2*((p21-p11)/prob)*(c2>=c1)
        }
        else{
            aa <- (1/(prob*(1-prob))-a)/(prob*(1-prob))*min(1,c0/abs(1/(prob*(1-prob))-a))*prob
        }
    }
    
    A <- 1/aa
    z <- a - size/(1-prob)

    list(A = A, z = z, c0 = c0)    
}
##################################################################
## IF of ML estimator
##################################################################
.getMLIF.binom <- function(size, prob){
    asVar <- (prob*(1-prob))/size
    names(asVar) <- "prob"
    A <- asVar
    a <- 0
    b <- Inf
    mse <- asVar
    names(mse) <- NULL
    bias <- sqrt(mse - asVar)
    names(bias) <- NULL
    param <- c(prob, size)
    names(param) <- c("prob", "size (known)")
    IFun <- function(x){}
    body(IFun) <- substitute({ 
        res <- matrix(x/size - prob, ncol = 1)
        colnames(res) <- c("IF for prob") 
        res 
    }, list(size = size, prob = prob))
    range <- function(alpha){} 
    body(range) <- substitute({
        rg <- qbinom(c(alpha/2, 1-alpha/2), size = size, prob = prob)
        seq(from = rg[1], to = rg[2], by = 1)
    }, list(size = size, prob = prob))
    
    IF <- list(model = "binom", modelName = "binomial probability", 
               parameter = param, A = A, a = a, b = b, IFun = IFun,
               range = range, asMSE = mse, asVar = asVar, asBias = bias)
    class(IF) <- "optIF"
    IF
}
##################################################################
## IF of minimum bias estimator
##################################################################
.getMBIF.binom <- function(size, prob){
    m0 <- qbinom(0.5, size=size, prob=prob)
    p1 <- pbinom(m0, size=size, prob=prob)
    p2 <- pbinom(m0-1, size=size-1, prob=prob)
    
    inte <- m0*(2*p1 - 1) + size*prob*(1-2*p2)
    b <- prob*(1-prob)/inte
    names(b) <- NULL
    
    p3 <- pbinom(m0-1, size = size, prob = prob)
    p4 <- pbinom(m0, size = size, prob = prob, lower.tail = FALSE)
    p5 <- dbinom(m0, size = size, prob = prob)
    beta <- (p3 - p4)/p5
    
    A <- 1
    a <- -size*prob/(prob*(1-prob))
    names(a) <- NULL
    mse <- Inf
    bias <- Inf
    asVar <- b^2
    names(asVar) <- "prob"
    
    param <- c(prob, size)
    names(param) <- c("prob", "size (known)")
    IFun <- function(x){}
    body(IFun) <- substitute({
        res <- matrix(bmin*((x > M) - (x < M) + beta*(x == M)), ncol = 1)
        colnames(res) <- "IF for prob"
        res 
    }, list(bmin = b, M = m0, beta = beta))
    range <- function(alpha){} 
    body(range) <- substitute({
        rg <- qbinom(c(alpha/2, 1-alpha/2), size = size, prob = prob)
        seq(from = rg[1], to = rg[2], by = 1)
    }, list(size = size, prob = prob))
    
    IF <- list(model = "binom", modelName = "binomial probability", 
               parameter = param, A = A, a = a, b = b, lowerCase = beta, 
               IFun = IFun, range = range, asMSE = mse, asVar = asVar, 
               asBias = bias)
    class(IF) <- "optIF"
    IF
}
##################################################################
## IF of optimally robust estimator
##################################################################
.getOptIF.binom <- function(size, prob, radius, A, a, b){
    asVar <- A - radius^2*b^2
    names(asVar) <- "prob"
    mse <- A
    names(mse) <- NULL
    bias <- radius*b
    names(bias) <- NULL
    param <- c(prob, size)
    names(param) <- c("prob", "size (known)")
    IFun <- function(x){}
    body(IFun) <- substitute({ 
        Y0 <- (x - size*prob)/(prob*(1-prob))
        Y <- A*Y0 - a
        ind1 <- abs(Y) < b
        res <- matrix(ind1*Y + (1-ind1)*b*sign(Y), ncol = 1)
        colnames(res) <- "IF for prob"
        res 
    }, list(size = size, prob = prob, A = A, a = a, b = b))
    range <- function(alpha){} 
    body(range) <- substitute({
        rg <- qbinom(c(alpha/2, 1-alpha/2), size = size, prob = prob)
        seq(from = rg[1], to = rg[2], by = 1)
    }, list(size = size, prob = prob))
    
    IF <- list(model = "binom", modelName = "binomial probability", 
               parameter = param, A = A, a = a, b = b, IFun = IFun,
               range = range, asMSE = mse, asVar = asVar, asBias = bias)
    class(IF) <- "optIF"
    IF
}