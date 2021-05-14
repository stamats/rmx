plotIF <- function(x, range.alpha = 1e-6, range.n = 501, 
                   info.digits = 2, param.digits = 2, 
                   ggplot.xlab = "x", ggplot.ylab = "IF(x)",
                   ggplot.ggtitle = NULL,
                   point.col = "#0072B5", point.alpha = 0.4, 
                   point.length.out = 5, point.range = c(1,7)){
  stopifnot(inherits(x, "rmx"))
  stopifnot(length(range.alpha) == 1)
  stopifnot(is.numeric(range.alpha))
  stopifnot(range.alpha > 0 & range.alpha < 1)

  y <- sort(unique(c(x$x, x$rmxIF$range(alpha = range.alpha, n = range.n))))
  IF <- x$rmxIF$IFun(y)
  IFmin <- min(IF)
  IFmax <- max(IF)
  IFnames <- colnames(IF)
  DF <- data.frame(y, IF)
  IFx <- x$rmxIF$IFun(x$x)
  if(ncol(IFx) == 1){
    info <- IFx
  }else{
    info <- sqrt(rowSums(IFx^2))
  }
  DFx <- data.frame(x = x$x, IFx, info = signif(info, info.digits))
  if(x$rmxIF$model %in% c("norm", "gamma")){
    if(ncol(DF) > 2){
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
      if(length(ggplot.ylab) == 1){
        ggplot.ylab <- rep(ggplot.ylab, ncol(DF)-1)
      }
      if(length(ggplot.ylab) != (ncol(DF)-1)){
        stop("'ggplot.ylab' must have length 1 or equal to number of parameters")
      }
      for(i in 1:(ncol(DF)-1)){
        if(is.null(ggplot.ggtitle)){
          ggt <- ggtitle(paste0(IFnames[i], " (", Param, ")"))
        }else{
          if(length(ggplot.ggtitle) == 1){
            ggplot.ggtitle <- rep(ggplot.ggtitle, ncol(DF)-1)
          }
          if(length(ggplot.ggtitle) != (ncol(DF)-1)){
            stop("'ggplot.ggtitle' must have length 1 or equal to number of parameters")
          }
          ggt <- ggtitle(ggplot.ggtitle[i])
        }
        gg[[i]] <- ggplot(DF, aes_string(x = "y", y = names(DF)[i+1])) +
          geom_line() + xlab(ggplot.xlab[i]) + ylab(ggplot.ylab[i]) + 
          ylim(c(IFmin, IFmax)) +
          geom_point(data = DFx, aes_string(x = "x", y = names(DFx)[i+1], 
                                            size = "info"), 
                     inherit.aes = FALSE, color = point.col, alpha = point.alpha) +
          scale_size(breaks = seq(from = min(DFx$info), to = max(DFx$info), 
                                  length.out = point.length.out), range = point.range) +
          ggt
      }
      grid.newpage()
      grid.draw(arrangeGrob(grobs = gg, ncol = ncol(DF)-1, nrow = 1))
    }else{
      Param <- paste(names(x$rmxIF$parameter), 
                     signif(x$rmxIF$parameter, param.digits), 
                     sep = " = ")
      if(is.null(ggplot.ggtitle)){
        ggt <- ggtitle(paste0(names(DF)[2], " (", Param, ")"))
      }else{
        if(length(ggplot.ggtitle) != 1){
          stop("'ggplot.ggtitle' must have length 1")
        }
        ggt <- ggtitle(ggplot.ggtitle)
      }
      gg <- ggplot(DF, aes_string(x = "y", y = names(DF)[2])) +
        geom_line() + xlab(ggplot.xlab) + ylab(ggplot.ylab) + 
        ylim(c(IFmin, IFmax)) +
        geom_point(data = DFx, aes_string(x = "x", y = names(DFx)[2], 
                                          size = "info"), 
                   inherit.aes = FALSE, color = point.col, alpha = point.alpha) +
        scale_size(breaks = c(min(DFx$info), max(DFx$info)), range = point.range) +
        ggt
      print(gg)
    }
  }
  invisible(gg)
}

