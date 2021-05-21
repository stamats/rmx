iiPlot <- function(x, ...){
  UseMethod("iiPlot")
}
iiPlot.rmx <- function(x, param.digits = 2, ggplot.xlim = NULL,
                       ggplot.xlab = "Absolute Information of ML", 
                       ggplot.ylab = "Absolute Information of RMX",
                       ggplot.ggtitle = NULL,
                       point.col = "#0072B5", point.alpha = 0.4, ...){
  stopifnot(length(ggplot.xlab) == 1)
  stopifnot(length(ggplot.ylab) == 1)
  
  ML <- rmx(x$x, model = x$rmxIF$model, eps = 0)#, message = FALSE)
  IFx.ML <- ML$rmxIF$IFun(x$x)
  IFx <- x$rmxIF$IFun(x$x)
  if(ncol(IFx) == 1){
    info <- IFx^2
    info.ML <- IFx.ML^2
  }else{
    info <- rowSums(IFx^2)
    info.ML <- rowSums(IFx.ML^2)
  }
  DF <- data.frame(x = info.ML, y = info)
  if(length(x$rmxIF$parameter > 1)){
    Param <- paste(paste(names(x$rmxIF$parameter), 
                         signif(x$rmxIF$parameter, param.digits), 
                         sep = " = "), collapse = ", ")
  }else{
    Param <- paste(names(x$rmxIF$parameter), 
                   signif(x$rmxIF$parameter, param.digits), 
                   sep = " = ")
  }
  if(is.null(ggplot.ggtitle)){
    ggt <- ggtitle(paste0("Absolute Information of RMX vs ML for ", x$rmxIF$model))
  }else{
    stopifnot(length(ggplot.ggtitle) == 1)
    ggt <- ggtitle(ggplot.ggtitle)
  }
  rg <- c(min(DF$x, DF$y), max(DF$y)*1.1)
  DFline <- data.frame(x = rg, y = rg)
  gg <- ggplot(DF, aes_string(x = "x", y = "y")) +
    xlab(ggplot.xlab) + ylab(ggplot.ylab) + 
    geom_vline(xintercept = x$rmxIF$b^2) +
    geom_point(color = point.col, alpha = point.alpha) + 
    geom_abline(slope = 1, intercept = 0)
    ggt
  if(!is.null(ggplot.xlim)) gg <- gg + xlim(ggplot.xlim)
  print(gg)
  invisible(gg)
}

