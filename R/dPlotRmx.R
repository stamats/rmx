dPlot <- function(x, ...){
  UseMethod("dPlot")
}
dPlot.rmx <- function(x, param.digits = 3, ggplot.xlab = "x", 
                      ggplot.ylab = "Density", ggplot.ggtitle = NULL, 
                      density.col = "#0072B5", density.lwd = 1, ...){
  stopifnot(length(ggplot.xlab) == 1)
  stopifnot(length(ggplot.ylab) == 1)
  
  Dname <- x$rmxIF$model
  param <- sapply(as.list(x$rmxIF$parameter), signif, digits = param.digits)
  
  if(is.null(ggplot.ggtitle)){
    Param <- paste(paste(names(param), param, sep = " = "), collapse = ", ")
    ggt <- ggtitle(paste0("Density-Plot for ", Dname, "(", Param, ")"))
  }else{
    stopifnot((length(ggplot.ggtitle) == 1))
    ggt <- ggtitle(ggplot.ggtitle)
  }
  if(x$rmxIF$model == "norm"){
    ggd <- stat_function(fun = dnorm, args = param, lwd = density.lwd)
  }
  
  DF <- data.frame(x = x$x)
  gg <- ggplot(DF, aes(x = x)) + 
    geom_rug() + ggd +
    geom_density(color = density.col, lwd = density.lwd) + 
    xlab(ggplot.xlab) + ylab(ggplot.ylab) + ggt
  print(gg)
  invisible(gg)
}
