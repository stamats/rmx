dPlot <- function(x, ...){
  UseMethod("dPlot")
}
dPlot.rmx <- function(x, param.digits = 3, ggplot.xlab = "x", 
                      ggplot.ylab = NULL, ggplot.ggtitle = NULL, 
                      density.col = "#0072B5", density.lwd = 1, 
                      density.n = 501, ...){
  stopifnot(length(ggplot.xlab) == 1)
  if(!is.null(ggplot.ylab)) stopifnot(length(ggplot.ylab) == 1)
  
  Dname <- x$rmxIF$model
  param <- sapply(as.list(x$rmxIF$parameter), signif, digits = param.digits)
  if(x$rmxIF$model == "binom") names(param)[2] <- "size"
  
  if(is.null(ggplot.ggtitle)){
    Param <- paste(paste(names(param), param, sep = " = "), collapse = ", ")
    ggt <- ggtitle(paste0("Density-Plot for ", Dname, "(", Param, ")"))
  }else{
    stopifnot((length(ggplot.ggtitle) == 1))
    ggt <- ggtitle(ggplot.ggtitle)
  }
  if(x$rmxIF$model == "norm"){
    if(is.null(ggplot.ylab)) ggplot.ylab <- "Density"
    ggd <- stat_function(fun = dnorm, args = param, lwd = density.lwd,
                         n = density.n)
    ggempD <- geom_density(color = density.col, lwd = density.lwd, bw = "SJ")
    DF <- data.frame(x = x$x)
    gg <- ggplot(DF, aes(x = x)) + geom_rug(color = density.col) + 
      ggempD + ggd + xlab(ggplot.xlab) + ylab(ggplot.ylab) + ggt +
      labs(caption = "Empirical density with rug plot") +
      theme(plot.caption = element_text(face = "bold", color = density.col))
  }
  if(x$rmxIF$model == "binom"){
    if(is.null(ggplot.ylab)) ggplot.ylab <- "Relative Frequency / Probability"
    size <- x$rmxIF$parameter["size (known)"]
    supp <- seq(from = 0, to = size, by = 1)
    y <- NULL
    DFsupp <- data.frame(x = supp, 
                         y = dbinom(supp, prob = x$rmxEst, size = size))
    ggd <- geom_bar(data = DFsupp, aes(x = x, y = y), inherit.aes = FALSE, 
                    stat = "identity")
    ggempD <- geom_point(color = density.col)
    DF <- data.frame(x = sort(unique(x$x)), 
                     y = as.vector(table(x$x)/length(x$x)))
    gg <- ggplot(DF, aes(x = x, y = y)) + ggd + ggempD + 
      geom_segment(aes(x=x, xend=x, y=0, yend=y), 
                   color = density.col, lwd = density.lwd) +
      xlab(ggplot.xlab) + ylab(ggplot.ylab) + ggt +
      labs(caption = "Observed relative frequency") +
      theme(plot.caption = element_text(face = "bold", color = density.col)) +
      scale_x_continuous(breaks = supp)
  }
  if(x$rmxIF$model == "pois"){
    if(is.null(ggplot.ylab)) ggplot.ylab <- "Relative Frequency / Probability"
    size <- qpois(1-1e-15, lambda = x$rmxEst)
    supp <- seq(from = 0, to = size, by = 1)
    y <- NULL
    DFsupp <- data.frame(x = supp, y = dpois(supp, lambda = x$rmxEst))
    ggd <- geom_bar(data = DFsupp, aes(x = x, y = y), inherit.aes = FALSE, 
                    stat = "identity")
    ggempD <- geom_point(color = density.col)
    DF <- data.frame(x = sort(unique(x$x)), 
                     y = as.vector(table(x$x)/length(x$x)))
    gg <- ggplot(DF, aes(x = x, y = y)) + ggd + ggempD + 
      geom_segment(aes(x=x, xend=x, y=0, yend=y), 
                   color = density.col, lwd = density.lwd) +
      xlab(ggplot.xlab) + ylab(ggplot.ylab) + ggt +
      labs(caption = "Observed relative frequency") +
      theme(plot.caption = element_text(face = "bold", color = density.col)) +
      scale_x_continuous(breaks = supp)
  }
  
  gg
}
