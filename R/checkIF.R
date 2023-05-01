checkIF <- function(x, ...){
  UseMethod("checkIF")
}
checkIF.optIF <- function(x, rel.tol = .Machine$double.eps^0.5, ...){
  if(x$model == "norm"){
    chs <- .checkIF.norm(IF = x, rel.tol = rel.tol)
  }
  if(x$model == "binom"){
    chs <- .checkIF.binom(IF = x)
  }
  if(x$model == "pois"){
    chs <- .checkIF.pois(IF = x)
  }
  
  check <- list(Fisher = chs$Fisher, center = chs$center, mse.eq = chs$mse.eq)
  class(check) <- "checkIF"
  check
}
checkIF.rmx <- function(x, rel.tol = .Machine$double.eps^0.5, ...){
  checkIF(x$rmxIF, rel.tol = rel.tol)
}
print.checkIF <- function(x, digits = getOption("digits"), prefix = " ", ...){
  cat("\n")
  cat(strwrap("Check of IF properties", prefix = prefix), sep = "\n")
  cat("\n")
  cat(format("Fisher consistency:", width = 20L, justify = "right"), "\n")
  print(x$Fisher, digits = digits, ...)
  cat("\n")
  cat(format("Centering:", width = 11L, justify = "right"), "\n")
  print(x$center, digits = digits, ...)
  cat("\n")
  cat(format("MSE equation:", width = 14L, justify = "right"), "\n")
  print(x$mse.eq, digits = digits, ...)
}
.checkIF.norm <- function(IF, rel.tol = .Machine$double.eps^0.5){
  a1a3 <- diag(IF$A)/IF$parameter["sd"]^2
  a1 <- a1a3[1]
  a3 <- a1a3[2]
  a2 <- IF$a[2]/IF$parameter["sd"]/a3 + 1
  b <- IF$b/IF$parameter["sd"]
  r <- IF$radius
  
  if(r == 0 || is.infinite(r)){
    IF.mean.x <- function(x, mean, sd){ 
      (x-mean)/sd*IF$IFun(x)[,1]*dnorm(x, mean = mean, sd = sd)/sd
    }
    IF.sd.x <- function(x, mean, sd){ 
      (((x-mean)/sd)^2-1)*IF$IFun(x)[,2]*dnorm(x, mean = mean, sd = sd)/sd
    }
    IF.mean <- function(x, mean, sd){ 
      IF$IFun(x)[,1]*dnorm(x, mean = mean, sd = sd) 
    }
    IF.sd <- function(x, mean, sd){ 
      IF$IFun(x)[,2]*dnorm(x, mean = mean, sd = sd) 
    }
    ch1 <- integrate(IF.mean.x, lower = -Inf, upper = Inf, rel.tol = rel.tol,
                     mean = IF$parameter["mean"], sd = IF$parameter["sd"])$value
    ch2 <- integrate(IF.sd.x, lower = -Inf, upper = Inf, rel.tol = rel.tol,
                     mean = IF$parameter["mean"], sd = IF$parameter["sd"])$value
    ch3 <- integrate(IF.mean, lower = -Inf, upper = Inf, rel.tol = rel.tol,
                     mean = IF$parameter["mean"], sd = IF$parameter["sd"])$value
    ch4 <- integrate(IF.sd, lower = -Inf, upper = Inf, rel.tol = rel.tol,
                     mean = IF$parameter["mean"], sd = IF$parameter["sd"])$value
    center <- c(ch3, ch4)
    names(center) <- c("mean", "sd")
    mse.eq <- NULL
  }else{
    integrand1 <- function(x, b, a1, a2, a3){
      x^2*.getw.norm(x, b, a1, a2, a3)*dnorm(x)
    }
    Int1 <- 2*integrate(integrand1, lower = 0, upper = Inf, 
                        rel.tol = rel.tol, b = b, a1 = a1, 
                        a2 = a2, a3 = a3)$value
    ch1 <- a1*Int1
    
    integrand2 <- function(x, b, a1, a2, a3){
      (x^2 - a2)^2*.getw.norm(x, b, a1, a2, a3)*dnorm(x)
    }
    Int2 <- 2*integrate(integrand2, lower = 0, upper = Inf, 
                        rel.tol = rel.tol, b = b, a1 = a1, 
                        a2 = a2, a3 = a3)$value
    ch2 <- a3*Int2
    
    integrand3 <- function(x, b, a1, a2, a3){
      (x^2 - a2)*.getw.norm(x, b, a1, a2, a3)*dnorm(x)
    }
    Int3 <- 2*integrate(integrand3, lower=0, upper=Inf, 
                        rel.tol = rel.tol, b = b, a1 = a1, 
                        a2 = a2, a3 = a3)$value
    ch3 <- a3*Int3
    
    ch4 <- .getr.norm(b = b, r = r, a1 = a1, a2 = a2, a3 = a3)
    
    center <- c(0, ch3)
    names(center) <- c("mean", "sd")
    mse.eq <- ch4
  }
  Fisher <- matrix(0, ncol = 2, nrow = 2)
  Fisher[1,1] <- ch1 - 1
  Fisher[2,2] <- ch2 - 1
  colnames(Fisher) <- rownames(Fisher) <- c("mean", "sd")
  
  list(Fisher = Fisher, center = center, mse.eq = mse.eq)
}
.checkIF.binom <- function(IF){
  prob <- IF$parameter["prob"]
  size <- IF$parameter["size (known)"]
  A <- IF$A
  a <- IF$a
  b <- IF$b
  r <- IF$radius
  
  supp <- seq(from = 0, to = size, by = 1)
  WS <- dbinom(supp, size = size, prob = prob)
  Y0 <- (supp - size*prob)/(prob*(1-prob))
  
  ch1 <- sum(Y0*IF$IFun(supp)*WS)
  ch2 <- sum(IF$IFun(supp)*WS)
  
  if(r == 0 || is.infinite(r) || ("lowerCase" %in% names(IF))){
    mse.eq <- NULL
  }else{
    ch3 <- sum(pmax(abs(A*Y0-a)-b, 0)*WS) - r^2*b
    mse.eq <- ch3
  }  

  Fisher <- ch1 - 1
  names(Fisher) <- "prob"
  center <- ch2
  names(center) <- "prob"
  
  list(Fisher = Fisher, center = center, mse.eq = mse.eq)
}
.checkIF.pois <- function(IF){
  lambda <- IF$parameter
  A <- IF$A
  a <- IF$a
  b <- IF$b
  r <- IF$radius
  
  supp <- seq(from = 0, to = qpois(1-1e-15, lambda = lambda), by = 1)
  WS <- dpois(supp, lambda = lambda)
  Y0 <- supp/lambda - 1
  
  ch1 <- sum(Y0*IF$IFun(supp)*WS)
  ch2 <- sum(IF$IFun(supp)*WS)
  
  if(r == 0 || is.infinite(r) || ("lowerCase" %in% names(IF))){
    mse.eq <- NULL
  }else{
    ch3 <- sum(pmax(abs(A*Y0-a)-b, 0)*WS) - r^2*b
    mse.eq <- ch3
  }  
  
  Fisher <- ch1 - 1
  names(Fisher) <- "prob"
  center <- ch2
  names(center) <- "prob"
  
  list(Fisher = Fisher, center = center, mse.eq = mse.eq)
}