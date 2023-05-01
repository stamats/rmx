riPlot <- function(x, ...){
  UseMethod("riPlot")
}
riPlot.rmx <- function(x, range.alpha = 1e-6, range.n = 501, 
                       info.digits = 2, param.digits = 2, 
                       ggplot.xlab = "x", ggplot.ylab = NULL, 
                       ggplot.ggtitle = NULL,
                       point.col = "#0072B5", point.alpha = 0.4, 
                       point.length.out = 5, point.range = c(1,7), 
                       plot = TRUE, ...){
  if(length(x$rmxEst) == 1)
    stop(paste0("Relative Information is only useful for models, ", 
                "where at least two parameters have to be estimated."))

  stopifnot(length(range.alpha) == 1)
  stopifnot(is.numeric(range.alpha))
  stopifnot(range.alpha > 0 & range.alpha < 0.5)

  if(x$rmxIF$model %in% c("norm", "gamma")){
    rg <- x$rmxIF$range(alpha = range.alpha, n = 2) 
  }
  y <- c(seq(from = min(rg[1], x$x), to = max(rg[2], x$x), 
             length.out = range.n), x$x)
  y <- sort(unique(y))
  IF <- x$rmxIF$IFun(y)
  IFx <- x$rmxIF$IFun(x$x)
  info.x <- rowSums(IFx^2)
  info <- rowSums(IF^2)
  relInfo <- IF^2/info
  relInfo.x <- IFx^2/info.x
  IFmin <- min(relInfo)
  IFmax <- max(relInfo)
  DF <- data.frame(y, relInfo)
  DFx <- data.frame(x = x$x, relInfo.x, 
                    info = signif(x = info.x, digits = info.digits))
  if(x$rmxIF$model %in% c("norm", "gamma")){
    gg <- vector(mode = "list", length = ncol(DF)-1)
    Param <- paste(paste(names(x$rmxIF$parameter), 
                         signif(x$rmxIF$parameter, param.digits), 
                         sep = " = "), collapse = ", ")
    if(length(ggplot.xlab) == 1){
      ggplot.xlab <- rep(ggplot.xlab, ncol(DF)-1)
    }
    if(length(ggplot.xlab) != (ncol(DF)-1)){
      stop("'ggplot.xlab' must have length 1 or equal to number of parameters")
    }
    ggplot.ylab.org <- ggplot.ylab
    ggplot.ylab <- character(length(x$rmxIF$parameter))
    for(i in 1:(ncol(DF)-1)){
      if(is.null(ggplot.ylab.org)){
        ggplot.ylab[i] <- paste0("Relative information for ", names(x$rmxIF$parameter)[i])
      }else{
        if(length(ggplot.ylab) == 1){
          ggplot.ylab <- rep(ggplot.ylab, ncol(DF)-1)
        }
        if(length(ggplot.ylab) != (ncol(DF)-1)){
          stop("'ggplot.ylab' must have length 1 or equal to number of parameters")
        }
      }
      if(is.null(ggplot.ggtitle)){
        ggt <- ggtitle(paste0(x$rmxIF$model, " (", Param, ")"))
      }else{
        if(length(ggplot.ggtitle) == 1){
          ggplot.ggtitle <- rep(ggplot.ggtitle, ncol(DF)-1)
        }
        if(length(ggplot.ggtitle) != (ncol(DF)-1)){
          stop("'ggplot.ggtitle' must have length 1 or equal to number of parameters")
        }
        ggt <- ggtitle(ggplot.ggtitle[i])
      }
      gg[[i]] <- ggplot(DF, aes(x = .data$y, y = .data[[names(DF)[i+1]]])) +
        geom_line() + xlab(ggplot.xlab[i]) + ylab(ggplot.ylab[i]) + 
        ylim(c(IFmin, IFmax)) +
        geom_point(data = DFx, aes(x = .data$x, y = .data[[names(DFx)[i+1]]], 
                                   size = .data$info), 
                   inherit.aes = FALSE, color = point.col, alpha = point.alpha) +
        scale_size(breaks = seq(from = min(DFx$info), to = max(DFx$info), 
                                length.out = point.length.out), range = point.range) +
        ggt
    }
    if(plot){
      grid.newpage()
      grid.draw(arrangeGrob(grobs = gg, ncol = ncol(DF)-1, nrow = 1)) 
    }
  }
  invisible(gg)
}

