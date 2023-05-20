###############################################################################
## optimal IF for exponential scale
###############################################################################
optIF.exp <- function(radius, scale = 1, aUp = 1, cUp = 0.97, delta = 1e-9){
    if(radius == 0){ # IF of ML estimator
        IF <- .getMLIF.exp(scale = scale)
        IF$radius <- 0
        return(IF)
    }
    if(radius >= 50){ # IF of minimum bias estimator
        IF <- .getMBIF.exp(scale = scale)
        IF$asBias <- radius*IF$b
        IF$asMSE <- IF$asVar + IF$asBias^2
        IF$radius <- radius
        
        return(IF)
    }
    
    ## compute Lagrange multipliers
    LM <- .getLM.exp(r0 = radius, aUp = aUp, cUp = cUp, delta = delta)
    A <- LM$A
    a <- LM$a
    b <- LM$b
    
    IF <- .getOptIF.exp(scale = scale, radius = radius, A = A, a = a, b = b)
    IF$radius <- radius
    IF
}

##################################################################
## centering
##################################################################
.geta.exp <- function(a, c0){
    c1 <- max(a-c0,0)
    c2 <- c0 + a
    (1 + c1 - (a-c0))*exp(-c1) - exp(-c2) - c0
}
##################################################################
## clipping
##################################################################
.getc.exp <- function(c0, r0, aUp, delta){
    c0 <- c0^2/(1-c0^2)
    a <- uniroot(.geta.exp, lower = 0, upper = aUp, tol = delta, c0=c0)$root 
    
    c1 <- max(a-c0,0)
    c2 <- c0 + a
    r1 <- (1 + c1 - (a-c0))*exp(-c1) + exp(-c2) + (a-c0) - 1
    r <- sqrt(r1/c0)
    
    r0-r
}
##################################################################
## Lagrange multipliers
##################################################################
.getLM.exp <- function(r0, aUp = 1, cUp = 0.97, delta = 1e-9){
    c0 <- uniroot(.getc.exp, lower = 0.01, upper = cUp, 
                  tol = delta, r0 = r0, aUp = aUp, delta = delta)$root
    c0 <- c0^2/(1-c0^2)
    a <- uniroot(.geta.exp, lower = 0, upper = aUp, tol = delta, c0 = c0)$root
    
    c1 <- max(a-c0,0)
    c2 <- c0 + a
    
    aa <- (c1*(2 + c1 - (a-c0)) + 2 - (a-c0))*exp(-c1) - (2+c2)*exp(-c2) - c0
    A <- 1/aa

    list(a=A*(a-1), A=A, b=A*c0)
}
##################################################################
## IF of ML estimator
##################################################################
.getMLIF.exp <- function(scale){
    asVar <- scale^2
    A <- asVar
    names(asVar) <- "lambda"
    a <- 0
    b <- Inf
    mse <- asVar
    names(mse) <- NULL
    bias <- sqrt(mse - asVar)
    names(bias) <- NULL
    param <- scale
    names(param) <- "scale"
    IFun <- function(x){}
    body(IFun) <- substitute({ 
        res <- matrix(x - scale, ncol = 1)
        colnames(res) <- c("IF for scale") 
        res 
    }, list(scale = scale))
    range <- function(alpha, n = 501){} 
    body(range) <- substitute({
        rg <- qexp(c(alpha/2, 1-alpha/2), rate = 1/scale)
        seq(from = rg[1], to = rg[2], length.out = n)
    }, list(scale = scale))
    
    IF <- list(model = "exp", modelName = "Exponential scale", 
               parameter = param, A = A, a = a, b = b, IFun = IFun,
               range = range, asMSE = mse, asVar = asVar, asBias = bias)
    class(IF) <- "optIF"
    IF
}
##################################################################
## IF of minimum bias estimator
##################################################################
.getMBIF.exp <- function(scale){
    b <- scale/log(2)
    m0 <- scale*log(2)

    A <- 1
    a <- (m0/scale - 1)/scale
    names(a) <- NULL
    mse <- Inf
    bias <- Inf
    asVar <- b^2
    names(asVar) <- "scale"
    
    param <- scale
    names(param) <- "scale"
    IFun <- function(x){}
    body(IFun) <- substitute({
        res <- matrix(bmin*((x > M) - (x < M)), ncol = 1)
        colnames(res) <- "IF for scale"
        res 
    }, list(bmin = b, M = m0))
    range <- function(alpha, n = 501){} 
    body(range) <- substitute({
        rg <- qexp(c(alpha/2, 1-alpha/2), rate = 1/scale)
        seq(from = rg[1], to = rg[2], length.out = n)
    }, list(scale = scale))
    
    IF <- list(model = "exp", modelName = "Exponential scale", 
               parameter = param, A = A, a = a, b = b, 
               IFun = IFun, range = range, asMSE = mse, asVar = asVar, 
               asBias = bias)
    class(IF) <- "optIF"
    IF
}
##################################################################
## IF of optimally robust estimator
##################################################################
.getOptIF.exp <- function(scale, radius, A, a, b){
    A <- scale^2*A
    a <- scale*a
    b <- scale*b
    asVar <- A - radius^2*b^2
    names(asVar) <- "scale"
    mse <- A
    names(mse) <- NULL
    bias <- radius*b
    names(bias) <- NULL
    param <- scale
    names(param) <- c("scale")
    IFun <- function(x){}
    body(IFun) <- substitute({ 
        Y0 <- (x/scale - 1)/scale
        Y <- A*Y0 - a
        ind1 <- abs(Y) < b
        res <- matrix(ind1*Y + (1-ind1)*b*sign(Y), ncol = 1)
        colnames(res) <- "IF for scale"
        res 
    }, list(scale = scale, A = A, a = a, b = b))
    range <- function(alpha, n = 501){} 
    body(range) <- substitute({
        rg <- qexp(c(alpha/2, 1-alpha/2), rate = 1/scale)
        seq(from = rg[1], to = rg[2], length.out = n)
    }, list(scale = scale))
    
    IF <- list(model = "exp", modelName = "Exponential scale", 
               parameter = param, A = A, a = a, b = b, IFun = IFun,
               range = range, asMSE = mse, asVar = asVar, asBias = bias)
    class(IF) <- "optIF"
    IF
}
