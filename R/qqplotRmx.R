qqplotRmx <- function(x, param.digits = 3, ggplot.xlab = "Theoretical Quantiles", 
                      ggplot.ylab = "Sample Quantiles", ggplot.ggtitle = NULL,
                      point.col = "#0072B5", point.alpha = 1){
  stopifnot(inherits(x, "rmx"))
  stopifnot(length(ggplot.xlab) == 1)
  stopifnot(length(ggplot.ylab) == 1)
  
  Dname <- x$rmxIF$model
  param <- sapply(as.list(x$rmxIF$parameter), signif, digits = param.digits)
  
  if(is.null(ggplot.ggtitle)){
    Param <- paste(paste(names(param), param, sep = " = "), collapse = ", ")
    ggt <- ggtitle(paste0("qq-Plot for ", Dname, "(", Param, ")"))
  }else{
    stopifnot((length(ggplot.ggtitle) == 1))
    ggt <- ggtitle(ggplot.ggtitle)
  }
  
  DF <- data.frame(x = x$x)
  gg <- ggplot(DF, aes(sample = x)) + 
    stat_qq_band(dparams = param, distribution = Dname, identity = TRUE) + 
    stat_qq_point(dparams = param, distribution = Dname, 
                  color = point.col, alpha = point.alpha) + 
    stat_qq_line(dparams = param, distribution = Dname, identity = TRUE) + 
    xlab(ggplot.xlab) + ylab(ggplot.ylab) + ggt
  print(gg)
  invisible(gg)
}
