optIF <- function(model = "norm", radius = NULL, ...){
  es.call <- match.call()
  
  stopifnot(is.character(model))
  stopifnot(length(model) == 1)

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
  
  listDots <- list(...)
  if(model == "norm"){ # normal distribution
    mean <- ifelse("mean" %in% names(listDots), listDots$mean, 0)
    sd <- ifelse("sd" %in% names(listDots), listDots$sd, 1)
    A.loc.start <- ifelse("A.loc.start" %in% names(listDots), 
                          listDots$A.loc.start, 1)
    A.sc.start <- ifelse("A.sc.start" %in% names(listDots), 
                          listDots$A.sc.start, 0.5)
    a.sc.start <- ifelse("a.sc.start" %in% names(listDots), 
                         listDots$a.sc.start, 0)
    bUp <- ifelse("bUp" %in% names(listDots), listDots$bUp, 1000)
    delta <- ifelse("delta" %in% names(listDots), listDots$delta, 1e-6)
    itmax <- ifelse("itmax" %in% names(listDots), listDots$itmax, 100L)
    IF <- optIF.norm(radius = radius, mean = mean, sd = sd, A.loc.start = A.loc.start, 
                     A.sc.start = A.sc.start, a.sc.start = a.sc.start, bUp = bUp, 
                     delta = delta, itmax = itmax)
  }
  if(model == "binom"){# binomial distribution
    size <- ifelse("size" %in% names(listDots), listDots$size, 1)
    prob <- ifelse("prob" %in% names(listDots), listDots$prob, 0.5) 
    aUp <- ifelse("aUp" %in% names(listDots), listDots$aUp, 100*size)
    cUp <- ifelse("cUp" %in% names(listDots), listDots$cUp, 1e4)
    delta <- ifelse("delta" %in% names(listDots), listDots$delta, 1e-9)
    IF <- optIF.binom(radius = radius, size = size, prob = prob, 
                      aUp = aUp, cUp = cUp, delta = delta)
  }
  if(model == "pois"){# Poisson distribution
    lambda <- ifelse("lambda" %in% names(listDots), listDots$lambda, 1) 
    aUp <- ifelse("aUp" %in% names(listDots), listDots$aUp, 100*lambda)
    cUp <- ifelse("cUp" %in% names(listDots), listDots$cUp, 1e4)
    delta <- ifelse("delta" %in% names(listDots), listDots$delta, 1e-9)
    IF <- optIF.pois(radius = radius, lambda = lambda, aUp = aUp, cUp = cUp, 
                     delta = delta)
  }
  if(model == "exp"){# Exponential distribution
    scale <- ifelse("scale" %in% names(listDots), listDots$scale, 1) 
    aUp <- ifelse("aUp" %in% names(listDots), listDots$aUp, 1)
    cUp <- ifelse("cUp" %in% names(listDots), listDots$cUp, 0.97)
    delta <- ifelse("delta" %in% names(listDots), listDots$delta, 1e-9)
    IF <- optIF.exp(radius = radius, scale = scale, aUp = aUp, cUp = cUp, 
                    delta = delta)
  }
  if(!model %in% c("norm", "binom", "pois", "exp")){
    stop("Given 'model' not yet implemented")
  }
  
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
  if("lowerCase" %in% names(x)){
    cat(format("Lower case beta:", width = 20L, justify = "right"), "\n")
    print(x$lowerCase, digits = digits, ...)
  }
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

plot.optIF <- function(x, alpha = 1e-6, digits = 2, plot = TRUE, n = 501, ...){
  if(x$model %in% c("norm", "gamma")){
    y <- x$range(alpha, n = n)
    IF <- x$IFun(y)
    IFmin <- min(IF)
    IFmax <- max(IF)
    IFnames <- colnames(IF)
    DF <- data.frame(y, IF)
    if(ncol(DF) > 2){
      gg <- vector(mode = "list", length = ncol(DF)-1)
      Param <- paste(paste(names(x$parameter), signif(x$parameter, digits), 
                           sep = " = "), collapse = ", ")
      for(i in 1:(ncol(DF)-1)){
        gg[[i]] <- ggplot(DF, aes(x = y, y = .data[[names(DF)[i+1]]])) +
          geom_line() + xlab("x") + ylab("IF(x)") + ylim(c(IFmin, IFmax)) +
          ggtitle(paste0(IFnames[i], " (", Param, ")"))
      }
      if(plot){
        grid.newpage()
        grid.draw(arrangeGrob(grobs = gg, ncol = ncol(DF)-1, nrow = 1))
      }
    }else{
      Param <- paste(names(x$parameter), signif(x$parameter, digits), 
                     sep = " = ")
      gg <- ggplot(DF, aes(x = y, y = .data[[names(DF)[2]]])) +
        geom_line() + xlab("x") + ylab("IF(x)") + ylim(c(IFmin, IFmax)) +
        ggtitle(paste0(IFnames, " (", Param, ")"))
    }
  }
  if(x$model %in% c("binom", "pois")){
    if(x$model == "binom"){
      y <- x$range(alpha = alpha)
    }
    if(x$model == "pois"){
      y <- x$range(alpha = alpha)
    }
    IF <- x$IFun(y)
    IFmin <- min(IF)
    IFmax <- max(IF)
    IFnames <- colnames(IF)
    DF <- data.frame(y, IF)
    Param <- paste(names(x$parameter), signif(x$parameter, digits), 
                   sep = " = ")
    gg <- ggplot(DF, aes(x = y, y = .data[[names(DF)[2]]])) +
      geom_point() + geom_line() + xlab("x") + ylab("IF(x)") + 
      ylim(c(IFmin, IFmax)) + ggtitle(paste0(IFnames, " (", Param, ")"))
    if(plot) print(gg)
  }
  if(x$model == "exp"){
    y <- x$range(alpha = alpha)
    IF <- x$IFun(y)
    IFmin <- min(IF)
    IFmax <- max(IF)
    IFnames <- colnames(IF)
    DF <- data.frame(y, IF)
    Param <- paste(names(x$parameter), signif(x$parameter, digits), 
                   sep = " = ")
    gg <- ggplot(DF, aes(x = y, y = .data[[names(DF)[2]]])) +
      geom_line() + xlab("x") + ylab("IF(x)") + 
      ylim(c(IFmin, IFmax)) + ggtitle(paste0(IFnames, " (", Param, ")"))
    if(plot) print(gg)
  }
  invisible(gg)
}
