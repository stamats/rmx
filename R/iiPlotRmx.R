iiPlot <- function(x, ...){
  UseMethod("iiPlot")
}
iiPlot.rmx <- function(x, param.digits = 2, ggplot.ylim = NULL,
                       ggplot.xlab = "Absolute Information of RMX", 
                       ggplot.ylab = "Absolute Information of ML",
                       ggplot.ggtitle = NULL, color.line = "#E18727",
                       point.col = "#0072B5", point.alpha = 0.4, ...){
  stopifnot(length(ggplot.xlab) == 1)
  stopifnot(length(ggplot.ylab) == 1)
  
  if(x$rmxIF$model %in% c("norm", "pois")){
    ML <- rmx(x$x, model = x$rmxIF$model, eps = 0, message = FALSE)
  }
  if(x$rmxIF$model == "binom"){
    ML <- rmx(x$x, model = x$rmxIF$model, eps = 0, message = FALSE,
              size = x$rmxIF$parameter["size (known)"])
  }
  
  IFx.ML <- ML$rmxIF$IFun(x$x)
  IFx <- x$rmxIF$IFun(x$x)
  if(ncol(IFx) == 1){
    info <- IFx^2
    info.ML <- IFx.ML^2
  }else{
    info <- rowSums(IFx^2)
    info.ML <- rowSums(IFx.ML^2)
  }
  DF <- data.frame(info.RMX = info, info.ML = info.ML)
  names(DF) <- c("info.RMX", "info.ML")
  if(length(x$rmxEst) > 1){
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
  gg <- ggplot(DF, aes(x = .data$info.RMX, y = .data$info.ML)) +
    xlab(ggplot.xlab) + ylab(ggplot.ylab) + 
    geom_abline(slope = 1, intercept = 0) +
    geom_point(color = point.col, alpha = point.alpha) + 
    geom_vline(xintercept = x$rmxIF$b^2, color = color.line) +
    ggt
  if(!is.null(ggplot.ylim)) gg <- gg + ylim(ggplot.ylim)
  gg
}
