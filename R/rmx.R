###############################################################################
## Compute RMX estimators for various models
###############################################################################
rmx <- function(x, model = "norm", eps.lower=0, eps.upper=NULL, eps=NULL, k = 3L, 
                initial.est=NULL, fsCor = TRUE, na.rm = TRUE, message = TRUE){
  es.call <- match.call()
  
  if(missing(x))
    stop("'x' is missing with no default")
  stopifnot(is.numeric(x))

  stopifnot(is.character(model))
  stopifnot(length(model) == 1)
  if(!is.na(pmatch(model, "norm")))  model <- "norm"
  MODELS <- c("norm", "binom", "pois", "gamma", "weibull")
  model <- pmatch(model, MODELS)
  if(is.null(eps) && is.null(eps.upper)){
    res0 <- rmx(x = x, model = MODELS[model], eps.upper = 0.5, k = k, 
                initial.est = initial.est, fsCor = fsCor, na.rm = na.rm,
                message = message)
    eps.lower <- 0
    eps.upper <- outlier(res0)$prop.outlier
    if(eps.upper == 0){
      eps <- 0
    }
  }
  
  if(is.null(eps)){
    if(length(eps.lower) != 1 || length(eps.upper) != 1)
      stop("'eps.lower' and 'eps.upper' have to be of length 1")
    if(!is.numeric(eps.lower) || !is.numeric(eps.upper) || eps.lower >= eps.upper) 
      stop("'eps.lower' < 'eps.upper' is not fulfilled")
    if((eps.lower < 0) || (eps.upper > 0.5))
      stop("'eps.lower' and 'eps.upper' have to be in [0, 0.5]")
  }else{
    if(length(eps) != 1)
      stop("'eps' has to be of length 1")
    if((eps < 0) || (eps > 0.5))
      stop("'eps' has to be in (0, 0.5]")
    if(eps == 0 && message){
      message("'eps=0': Maximum likelihood estimator is computed.")
    }
  }
  
  stopifnot(is.numeric(k))
  if(!is.integer(k))  k <- as.integer(k)
  if(k < 1){
    stop("'k' has to be some positive integer value")
  }
  stopifnot(length(k) == 1)
  stopifnot(length(fsCor) == 1)
  stopifnot(is.logical(fsCor))
  stopifnot(length(na.rm) == 1)
  stopifnot(is.logical(na.rm))
  
  completecases <- !is.na(x)
  x.org <- x
  if(na.rm) x <- x[completecases] 
  
  if(model == 1){ # normal distribution
    RMXest <- rmx.norm(x, eps.lower=eps.lower, eps.upper=eps.upper, eps=eps, 
                       k = k, initial.est=initial.est, fsCor = fsCor)
  }
  if(model != 1){
    stop("Given 'model' not yet implemented")
  }
  if(length(x) < length(x.org)){
    RMXest$Infos <- rbind(RMXest$Infos, 
                          c("rmx", "'NA' values have been removed before the computation"))
  }

  RMXest$rmxIF$call <- es.call
  RMXest$x <- x.org
  RMXest$n <- length(x)
  RMXest$eps.lower <- ifelse(is.null(eps), eps.lower, NA)
  RMXest$eps.upper <- ifelse(is.null(eps), eps.upper, NA)
  RMXest$eps <- ifelse(is.null(eps), NA, eps)
  RMXest$fsCor <- fsCor
  RMXest$k <- k
  RMXest$call <- es.call
  RMXest
}

## adapted from MASS:::print.fitdistr
print.rmx <- function(x, digits = getOption("digits"), ...){
  cat("\n")
  cat(strwrap(paste0("RMX estimator for ", x$rmxIF$modelName), 
              prefix = " "), sep = "\n")
  cat("\n")
  SD <- sqrt(diag(x$rmxIF$asVar))/sqrt(x$n)
  ans <- format(rbind(x$rmxEst, SD), digits = digits)
  ans[1L, ] <- sapply(ans[1L, ], function(x) paste("", x))
  ans[2L, ] <- sapply(ans[2L, ], function(x) paste("(", x, ")", sep = ""))
  dn <- dimnames(ans)
  dn[[1L]] <- rep("", 2L)
  dn[[2L]] <- paste(substring("      ", 1L, (nchar(ans[2L,]) - nchar(dn[[2L]]))%/%2), dn[[2L]])
  dn[[2L]] <- paste(dn[[2L]], substring("      ", 1L, (nchar(ans[2L,]) - nchar(dn[[2L]]))%/%2))
  dimnames(ans) <- dn
  print(ans, quote = FALSE)
  cat("\n NOTE: asymptotic standard errors are shown\n")
  cat("\n Call:\n ", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  invisible(x)
}

summary.rmx <- function(object, digits = getOption("digits"), ...){
  cat("\n")
  cat(strwrap(paste0("RMX estimator for ", object$rmxIF$modelName), 
              prefix = " "), sep = "\n")
  cat("\n")
  SD <- sqrt(diag(object$rmxIF$asVar))/sqrt(object$n)
  ans <- format(rbind(object$rmxEst, SD), digits = digits)
  ans[1L, ] <- sapply(ans[1L, ], function(x) paste("", x))
  ans[2L, ] <- sapply(ans[2L, ], function(x) paste("(", x, ")", sep = ""))
  dn <- dimnames(ans)
  dn[[1L]] <- rep("", 2L)
  dn[[2L]] <- paste(substring("      ", 1L, (nchar(ans[2L,]) - nchar(dn[[2L]]))%/%2), dn[[2L]])
  dn[[2L]] <- paste(dn[[2L]], substring("      ", 1L, (nchar(ans[2L,]) - nchar(dn[[2L]]))%/%2))
  dimnames(ans) <- dn
  print(ans, quote = FALSE)
  cat("\n NOTE: asymptotic standard errors are shown\n")
  cat("\n")
  cat(paste(format("Sample size", width = 24L, justify = "right"), 
            format(object$n, digits = digits), sep = " = "), sep = "\n")
  if(is.na(object$eps)){
    cat(paste(format("Amount gross-errors", width = 24L, justify = "right"), 
              paste0(100*signif(object$eps.lower, digits = 3), " - ",
                     100*signif(object$eps.upper, digits = 3), " %"), sep = " = "), 
        sep = "\n")
  }else{
    cat(paste(format("Amount gross-errors", width = 24L, justify = "right"), 
              paste0(100*signif(object$eps, digits = 3), " %"), sep = " = "), 
              sep = "\n")
  }
  if(object$fsCor){
    cat(paste(format("FS-corrected radius", width = 24L, justify = "right"), 
              format(object$rmxIF$radius, digits = 3), sep = " = "), sep = "\n")
  }else{
    cat(paste(format("Infinitesimal radius", width = 24L, justify = "right"), 
              format(object$rmxIF$radius, digits = 3), sep = " = "), sep = "\n")
  }
  cat(paste(format("Maximum asymptotic MSE", width = 24L, justify = "right"), 
            format(object$rmxIF$asMSE, digits = 3), sep = " = "), sep = "\n")
  cat(paste(format("Maximum asymptotic bias", width = 24L, justify = "right"), 
            format(object$rmxIF$asBias, digits = 3), sep = " = "), sep = "\n")
  cat(format("Asymptotic variance:", width = 24L, justify = "right"), "\n")
  print(object$rmxIF$asVar, digits = 3, ...)
  cat("\n Call:\n ", paste(deparse(object$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  invisible(object)
}

coef.rmx <- function(object, complete = TRUE, ...){
  cf <- object$rmxEst
  if (complete) 
    cf
  else cf[!is.na(cf)]
}
vcov.rmx <- function(object, ...){
    object$rmxIF$asVar/object$n
}
## from stats:::format.perc
.format.perc <- function (probs, digits){
  paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits), "%")
}
confint.rmx <- function (object, parm, level = 0.95, method = "as", R = 9999, 
                         parallel = FALSE, ncores = NULL, ...){
  if(method == "as"){
    Method <- "Asymptotic (LAN-based) confidence interval"
    ci <- confint.default(object)
  }
  if(method == "as.bias"){
    Method <- "Asymptotic (LAN-based), uniform (bias-aware) confidence interval"
    cf <- coef(object)
    pnames <- names(cf)
    if (missing(parm)) 
      parm <- pnames
    else if (is.numeric(parm)) 
      parm <- pnames[parm]
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    pct <- .format.perc(a, 3)
    fac <- qnorm(a, mean = c(-object$rmxIF$asBias, object$rmxIF$asBias))
    ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, pct))
    ses <- sqrt(diag(vcov(object)))[parm]
    ci[] <- cf[parm] + ses %o% fac
  }
  if(method == "boot"){
    Method <- "Bootstrap confidence interval"
    n <- length(object$x)
    t0 <- c(object$rmxEst, diag(object$rmxIF$asVar))
    seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    X <- matrix(sample(object$x, size = R*n, replace = TRUE), nrow = R)
    boot.res <- rowRmx(X, model = object$rmxIF$model, computeSE = TRUE,
                       parallel = parallel, ncores = ncores)
    t <- cbind(boot.res$rmxEst, n*boot.res$asSE^2)
    boot.out <- list(t0 = t0, t = t, R = R, data = object$x, seed = seed, 
                     statistic = function(x, i){}, sim = "ordinary", 
                     call = boot.res$call, stype = "i", 
                     strata = rep(1, length(object$x)),
                     weights = rep(1/length(object$x), length(object$x)))
    class(boot.out) <- "boot"
    attr(boot.out, "boot_type") <- "boot"
    if(object$rmxIF$model %in% c("norm", "gamma")){
      ci <- list(boot.ci(boot.out, index = c(1,3)),
                 boot.ci(boot.out, index = c(2,4)))
      names(ci) <- names(object$rmxIF$parameter)
    }
  }
  attr(ci, "conf.level") <- level
  CI <- list(method = Method, conf.int = ci, rmxEst = object$rmxEst)
  class(CI) <- "rmxCI"
  CI
}
## adapted from MKinfer:::print.confint
print.rmxCI <- function (x, digits = getOption("digits"), prefix = " ", ...){
  cat("\n")
  cat(strwrap(x$method, prefix = prefix), sep = "\n")
  cat("\n")
  out <- x$conf.int
  attr(out, "conf.level") <- NULL
  if (is.list(out) | inherits(out, "bootci")) {
    print(out, digits = digits, ...)
  }
  else {
    if (nrow(x$conf.int) > 1) {
      cat(format(100 * attr(x$conf.int, "conf.level")), 
          " percent confidence intervals:\n", sep = "")
    }
    else {
      cat(format(100 * attr(x$conf.int, "conf.level")), 
          " percent confidence interval:\n", sep = "")
    }
    print(out, digits = digits, ...)
  }
  if (!is.null(x$rmxEst)) {
    cat("\n")
    if (length(x$rmxEst) == 1) 
      cat("RMX estimate:\n")
    if (length(x$rmxEst) > 1) 
      cat("RMX estimates:\n")
    print(x$rmxEst, digits = digits, ...)
  }
  cat("\n")
  invisible(x)
}
plot.rmx <- function(x, which = 1, 
                     control = list(ifPlot = NULL, qqPlot = NULL,
                                    ppPlot = NULL, dPlot = NULL,
                                    aiPlot = NULL, riPlot = NULL,
                                    iiPlot = NULL), ...){
  show <- rep(FALSE, 7)
  show[which] <- TRUE
  if(show[1]){
    do.call("ifPlot", args = c(list(x = x), control$ifPlot))
  }
  if(show[2]){
    do.call("qqPlot", args = c(list(x = x), control$qqPlot))
  }
  if(show[3]){
    do.call("ppPlot", args = c(list(x = x), control$ppPlot))
  }
  if(show[4]){
    do.call("dPlot", args = c(list(x = x), control$dPlot))
  }
  if(show[5]){
    do.call("aiPlot", args = c(list(x = x), control$aiPlot))
  }
  if(show[6]){
    do.call("riPlot", args = c(list(x = x), control$riPlot))
  }
  if(show[7]){
    do.call("iiPlot", args = c(list(x = x), control$iiPlot))
  }
}
