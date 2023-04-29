aiPlot <- function(x, ...){
  UseMethod("aiPlot")
}
aiPlot.rmx <- function(x, range.alpha = 1e-6, range.n = 501, 
                       param.digits = 2, ggplot.xlab = "x", 
                       ggplot.ylab = expression(paste(abs(IF(x))^2)),
                       ggplot.ggtitle = NULL,
                       point.col = "#0072B5", point.alpha = 0.4, ...){
  stopifnot(length(range.alpha) == 1)
  stopifnot(is.numeric(range.alpha))
  stopifnot(range.alpha > 0 & range.alpha < 0.5)
  stopifnot(length(ggplot.xlab) == 1)
  stopifnot(length(ggplot.ylab) == 1)
  
  if(x$rmxIF$model %in% c("norm", "gamma")){
    rg <- x$rmxIF$range(alpha = range.alpha, n = 2)
    y <- c(seq(from = min(rg[1], x$x), to = max(rg[2], x$x), length.out = range.n), 
           x$x)
    y <- sort(unique(y))
  }
  
  if(x$rmxIF$model == "binom"){
    y <- x$rmxIF$range(alpha = 0)
  }
  if(x$rmxIF$model == "pois"){
    y <- x$rmxIF$range(alpha = 1e-15)
  }
  
  IF <- x$rmxIF$IFun(y)
  IFx <- x$rmxIF$IFun(x$x)
  if(ncol(IFx) == 1){
    info.x <- IFx^2
    info <- IF^2
  }else{
    info.x <- rowSums(IFx^2)
    info <- rowSums(IF^2)
  }
  IFmin <- min(info)
  IFmax <- max(info)
  DF <- data.frame(y, info)
  names(DF) <- c("y", "info")
  DFx <- data.frame(x$x, info.x)
  names(DFx) <- c("x", "info")
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
    ggt <- ggtitle(paste0("Absolute Information of ", x$rmxIF$model, 
                          " (", Param, ")"))
  }else{
    stopifnot(length(ggplot.ggtitle) == 1)
    ggt <- ggtitle(ggplot.ggtitle)
  }
  gg <- ggplot(DF, aes_string(x = "y", y = "info")) +
    geom_line() + xlab(ggplot.xlab) + ylab(ggplot.ylab) + 
    ylim(c(IFmin, IFmax)) +
    geom_point(data = DFx, aes_string(x = "x", y = "info"), 
               inherit.aes = FALSE, color = point.col, alpha = point.alpha) +
    ggt
  gg
}

