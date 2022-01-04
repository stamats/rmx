ppPlot <- function(x, ...){
  UseMethod("ppPlot")
}
ppPlot.rmx <- function(x, param.digits = 3, 
                       ggplot.xlab = "Theoretical Cumulative Probabilities", 
                       ggplot.ylab = "Empirical Cumulative Probabilities", 
                       ggplot.ggtitle = NULL,
                       point.col = "#0072B5", point.alpha = 1, ...){
  stopifnot(length(ggplot.xlab) == 1)
  stopifnot(length(ggplot.ylab) == 1)
  
  Dname <- x$rmxIF$model
  if(Dname %in% c("binom", "pois"))
    stop("'ppPlot' is only implemented for continous models.")
    
  param <- sapply(as.list(x$rmxIF$parameter), signif, digits = param.digits) 

  if(is.null(ggplot.ggtitle)){
    Param <- paste(paste(names(param), param, sep = " = "), collapse = ", ")
    ggt <- ggtitle(paste0("pp-Plot for ", Dname, "(", Param, ")"))
  }else{
    stopifnot((length(ggplot.ggtitle) == 1))
    ggt <- ggtitle(ggplot.ggtitle)
  }
  
  DF <- data.frame(x = x$x)
  gg <- ggplot(DF, aes(sample = x)) + 
    stat_pp_band(dparams = param, distribution = Dname) + 
    stat_pp_point(dparams = param, distribution = Dname, 
                  color = point.col, alpha = point.alpha) + 
    stat_pp_line() + xlab(ggplot.xlab) + ylab(ggplot.ylab) + ggt
  gg
}
