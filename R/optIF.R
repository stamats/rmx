optIF <- function(model = "norm", radius = NULL, check = FALSE, 
                  delta = 1e-6, itmax = 100L, ...){
  es.call <- match.call()
  
  stopifnot(is.character(model))
  stopifnot(length(model) == 1)
  if(!is.na(pmatch(model, "norm")))  model <- "norm"
  MODELS <- c("norm", "binom", "pois", "gamma", "weibull")
  model <- pmatch(model, MODELS)
  
  if(is.null(radius))
    stop("Please specify 'radius'")
  if(length(radius) != 1)
    stop("'radius' has to be of length 1")
  if(radius < 0)
    stop("'radius' has to be in (0, Inf]")
  if(radius == 0){
    message("'radius=0': IF of the maximum likelihood estimator is computed.")
  }
  if(radius == Inf)
    message("'radius=Inf': IF of the minimum bias estimator is computed.")
  
  stopifnot(length(check) == 1)
  stopifnot(is.logical(check))
  
  stopifnot(length(delta) == 1)
  stopifnot(is.numeric(delta))
  if(delta <= 0)
    stop("'delta' must be positive")
  if(delta > 0.1)
    warning("'delta' is expected to be small or very small.")
  
  stopifnot(is.numeric(itmax))
  if(!is.integer(itmax))  itmax <- as.integer(itmax)
  if(itmax < 1){
    stop("'itmax' has to be some positive integer value")
  }
  
  if(model == 1){ # normal distribution
    IF <- optIF.norm(radius, check = check, delta = delta, itmax = itmax, ...)
  }
  if(model != 1){
    stop("Given 'model' not yet implemented")
  }
  
  IF$radius <- radius
  IF$call <- es.call
  IF
}

print.optIF <- function(x, digits = getOption("digits"), prefix = " ", ...){
  cat("\n")
  if(x$radius == 0){
    cat(strwrap(paste0("IF of ML estimator for ", x$modelName), 
                prefix = prefix), sep = "\n")
  } 
  if(x$radius == Inf){
    cat(strwrap(paste0("IF of minimum bias estimator for ", x$modelName), 
                prefix = prefix), sep = "\n")
  }
  if(x$radius > 0 && is.finite(x$radius)){
    cat(strwrap(paste0("IF of optimally-robust estimator for ", x$modelName), 
                prefix = prefix), sep = "\n")
  }
  cat("\n")
  for(i in 1:length(x$parameter)){
    cat(paste(format(names(x$parameter)[i], width = 20L, justify = "right"), 
              format(x$parameter[i], digits = digits), sep = " = "), sep = "\n")
  }
  cat(paste(format("Radius", width = 20L, justify = "right"), 
            format(x$radius, digits = digits), sep = " = "), sep = "\n")
  cat(paste(format("Clipping constant", width = 20L, justify = "right"), 
            format(x$b, digits = digits), sep = " = "), sep = "\n")
  cat("\n")
  cat(format("Centering vector:", width = 20L, justify = "right"), "\n")
  print(x$a, digits = digits, ...)
  cat("\n")
  cat(format("Standardising matrix:", width = 20L, justify = "right"), "\n")
  print(x$A, digits = digits, ...)
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
}

summary.optIF <- function(object, digits = getOption("digits"), prefix = " ", ...){
  cat("\n")
  if(object$radius == 0){
    cat(strwrap(paste0("IF of ML estimator for ", object$modelName), 
                prefix = prefix), sep = "\n")
  } 
  if(object$radius == Inf){
    cat(strwrap(paste0("IF of minimum bias estimator for ", object$modelName), 
                prefix = prefix), sep = "\n")
  }
  if(object$radius > 0 && is.finite(object$radius)){
    cat(strwrap(paste0("IF of optimally-robust estimator for ", object$modelName), 
                prefix = prefix), sep = "\n")
  }
  cat("\n")
  for(i in 1:length(object$parameter)){
    cat(paste(format(names(object$parameter)[i], width = 24L, justify = "right"), 
              format(object$parameter[i], digits = digits), sep = " = "), sep = "\n")
  }
  cat(paste(format("Infinitesimal radius", width = 24L, justify = "right"), 
            format(object$radius, digits = digits), sep = " = "), sep = "\n")
  cat(paste(format("Maximum asymptotic MSE", width = 24L, justify = "right"), 
            format(object$asMSE, digits = digits), sep = " = "), sep = "\n")
  cat(paste(format("Maximum asymptotic bias", width = 24L, justify = "right"), 
            format(object$asBias, digits = digits), sep = " = "), sep = "\n")
  cat(format("Asymptotic variance:", width = 24L, justify = "right"), "\n")
  print(object$asVar, digits = digits, ...)
  cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
}

plot.optIF <- function(x, alpha = 1e-6, digits = 2, ...){
  y <- x$range(alpha, ...)
  IF <- x$IFun(y)
  IFmin <- min(IF)
  IFmax <- max(IF)
  IFnames <- colnames(IF)
  DF <- data.frame(y, IF)
  if(x$model %in% c("norm", "gamma")){
    if(ncol(DF) > 2){
      gg <- vector(mode = "list", length = ncol(DF)-1)
      Param <- paste(paste(names(x$parameter), signif(x$parameter, digits), 
                           sep = " = "), collapse = ", ")
      for(i in 1:(ncol(DF)-1)){
        gg[[i]] <- ggplot(DF, aes_string(x = "y", y = names(DF)[i+1])) +
          geom_line() + xlab("x") + ylab("IF(x)") + ylim(c(IFmin, IFmax)) +
          ggtitle(paste0(IFnames[i], " (", Param, ")"))
      }
      grid.draw(arrangeGrob(grobs = gg, ncol = ncol(DF)-1, nrow = 1))
    }else{
      Param <- paste(names(x$parameter), signif(x$parameter, digits), 
                     sep = " = ")
      gg <- ggplot(DF, aes_string(x = "y", y = names(DF)[2])) +
        geom_line() + xlab("x") + ylab("IF(x)") + ylim(c(IFmin, IFmax)) +
        ggtitle(paste0(names(DF)[2], " (", Param, ")"))
      print(gg)
    }
  }
  invisible(gg)
}
