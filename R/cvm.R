###############################################################################
## Compute Cramer von Mises minimum distance estimators for various models
###############################################################################
## main function
cvm <- function(x, model = "norm", mu = "model", na.rm = TRUE, ...){
  if(missing(x))
    stop("'x' is missing with no default")
  stopifnot(is.numeric(x))
  
  stopifnot(is.character(model))
  stopifnot(length(model) == 1)

  stopifnot(is.character(mu))
  stopifnot(length(mu) == 1)
  stopifnot(mu %in% c("model", "data"))

  stopifnot(length(na.rm) == 1)
  stopifnot(is.logical(na.rm))
  
  completecases <- !is.na(x)
  x.org <- x
  if(na.rm) x <- x[completecases] 
  
  if(model == "norm"){# normal distribution
    listDots <- list(...)
    if("startPar" %in% names(listDots)){
      startPar <- listDots$startPar
    }else{
      startPar <- NULL
    }
    est <- cvm.norm(x = x, mu = mu, startPar = startPar)
  }
  if(model == "binom"){# binomial distribution
    listDots <- list(...)
    if(!("size" %in% names(listDots)))
      stop("Parameter 'size' must be specified!")
    size <- listDots$size
    est <- cvm.binom(x = x, mu = mu, size = size)
  }
  if(model == "pois"){# Poisson distribution
    est <- cvm.pois(x = x, mu = mu)
  }
  est
}
## helper function to compute Cramer von Mises distance
.cvmdist <- function(x, pfun, dfun, mu, abscont, supp){
  if(mu == "data"){
    ecdf.x <- ecdf(x)
    Dist <- sqrt(mean((ecdf.x(x) - pfun(x))^2))
  }
  if(mu == "model"){
    if(abscont){
      ## code adapted from by P. Ruckdeschel, file CvMDist.R in package distrEx
      x1 <- sort(x)
      p <- pfun
      p0 <- pfun(x1)
      p1 <- 1-p0
      p2 <- 1-p0^2
      n <- length(x1)
      i1 <- 2*(1:n)-1
      Dist <- sqrt(mean(i1*p1)/n-mean(p2)+1/3)
    }else{
      ecdf.x <- ecdf(x)
      Dist <- sqrt(sum((ecdf.x(supp) - pfun(supp))^2*dfun(supp)))
    }
  }
  Dist
}
cvm.norm <- function(x, mu, startPar = NULL){
  if(is.null(startPar)) startPar <- c(median(x), mad(x))
  fn <- function(par, Data, mu){
    if(par[2] <= 0) par[2] <- abs(par[2]) + .Machine$double.eps
    M <- S <- NULL
    pfun <- function(q){
      pnorm(q, mean = M, sd = S)
    }
    body(pfun) <- substitute({ pnorm(q, mean = M, sd = S) },
                             list(M = par[1], S = par[2]))
    dfun <- function(x){
      dnorm(x, mean = M, sd = S)
    }
    body(dfun) <- substitute({ dnorm(x, mean = M, sd = S) },
                             list(M = par[1], S = par[2]))
    .cvmdist(x = Data, pfun = pfun, dfun = dfun, mu = mu, 
             abscont = TRUE, supp = NULL)
  }
  est <- optim(par = startPar, fn = fn, Data = x, mu = mu)$par
  names(est) <- c("mean", "sd")
  est
}
cvm.binom <- function(x, mu, size){
  f <- function(par, Data, mu){
    S <- P <- NULL
    pfun <- function(q){
      pbinom(q, size = S, prob = P)
    }
    body(pfun) <- substitute({ pbinom(q, size = S, prob = P) },
                             list(S = size, P = par))
    dfun <- function(x){
      dbinom(x, size = S, prob = P)
    }
    body(dfun) <- substitute({ dbinom(x, size = S, prob = P) },
                             list(S = size, P = par))
    .cvmdist(x = Data, pfun = pfun, dfun = dfun, mu = mu, 
             abscont = FALSE, supp = seq(from = 0, to = size, by = 1))
  }
  est <- optimize(f = f, interval = c(.Machine$double.eps, 1 - .Machine$double.eps),
                  Data = x, mu = mu)$minimum
  names(est) <- "prob"
  est
}
cvm.pois <- function(x, mu){
  f <- function(par, Data, mu){
    L <- NULL
    pfun <- function(q){
      ppois(q, lambda = L)
    }
    body(pfun) <- substitute({ ppois(q, lambda = L) }, list(L = par))
    dfun <- function(x){
      dpois(x, lambda = L)
    }
    body(dfun) <- substitute({ dpois(x, lambda = L) }, list(L = par))
    .cvmdist(x = Data, pfun = pfun, dfun = dfun, mu = mu, 
             abscont = FALSE, supp = seq(from = 0, 
                                         to = qpois(1-.Machine$double.eps, 
                                                    lambda = par)))
  }
  est <- optimize(f = f, 
                  interval = c(.Machine$double.eps, max(x)),
                  Data = x, mu = mu)$minimum
  names(est) <- "lambda"
  est
}
rowCVM <- function(x, model, mu = "model", na.rm = TRUE, parallel = FALSE, 
                   ncores = NULL, ...){
  if(parallel){
    if(is.null(ncores)){
      cores <- detectCores()
      cl <- makeCluster(cores[1]-1)
    }else{
      cl <- makeCluster(ncores)
    }
    clusterExport(cl, list(".cvmdist", "cvm.norm", "cvm.binom", "cvm.pois"), 
                  envir = environment(fun = rowCVM))
    res <- parApply(cl = cl, X = x, MARGIN = 1, FUN = cvm, model = model, 
                     mu = mu, na.rm = na.rm, ...)
    stopCluster(cl)
  }else{
    res <- apply(x, 1, cvm, model = model, mu = mu, na.rm = na.rm, ...)
  }
  res
}