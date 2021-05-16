cniper <- function(x, ...){
  UseMethod("cniper")
}
cniper.rmx <- function(x, range.alpha = 1e-6, ...){
  stopifnot(length(range.alpha) == 1)
  stopifnot(is.numeric(range.alpha))
  stopifnot((range.alpha > 0) && (range.alpha < 0.5))
  if(x$rmxIF$model == "norm"){
    SD <- x$rmxEst["sd"]
    tr.invF <- 1.5*SD^2
    maxMSE <- x$rmxIF$asMSE
    Delta <- sqrt(maxMSE - tr.invF)/x$rmxIF$radius
    fun <- function(y){
      z <- (y - x$rmxEst["mean"])/x$rmxEst["sd"]
      Y <- x$rmxEst["sd"]*cbind(z, 0.5*(z^2 - 1)) 
      sqrt(as.vector(Y %*% t(Y))) - Delta
    }
    rg <- x$rmxIF$range(alpha = range.alpha, n = 2)
    pt.lo <- uniroot(fun, lower = rg[1], upper = x$rmxEst["mean"])$root
    d.lo <- x$rmxEst["mean"] - pt.lo
    pt.up <- x$rmxEst["mean"] + d.lo
    names(pt.up) <- NULL
    prop.lo <- sum(x$x < pt.lo)/x$n
    prop.up <- sum(x$x > pt.up)/x$n
    p.lo <- pnorm(pt.lo, mean = x$rmxEst["mean"], sd = x$rmxEst["sd"])
    p.up <- p.lo
    res <- list(rmx = x,
                lower = pt.lo, upper = pt.up, prop.cniper = prop.lo+prop.up, 
                p.cniper = 2*p.lo, prop.lower = prop.lo, prop.upper = prop.up, 
                p.lower = p.lo, p.upper = p.up)
    class(res) <- "cniper"
  }else{
    stop("not yet implemented")
  }
  res
}
getCnipers <- function(x, ...){
  if(!inherits(x, "cniper")) x <- cniper(x, ...)
  if(x$rmx$rmxIF$model == "norm"){
    out <- c(which(x$rmx$x < x$lower), which(x$rmx$x > x$upper))
    val <- x$rmx$x[out]
  }
  ind <- order(val)
  list(values = val[ind], indices = out[ind])
}
print.cniper <- function(x, digits = 3, ...){
  cat("\n")
  cat("     ", strwrap(paste0("Cniper region for ", x$rmx$rmxIF$modelName), 
              prefix = " "), "\n")
  cat("\n")
  cat(paste(format("Parameter:", width = 30L, justify = "right"), 
            paste0(paste(names(x$rmx$rmxIF$parameter), 
                   signif(x$rmx$rmxIF$parameter, digits = digits), 
                   sep = " = "), collapse = ", ")), "\n")
  if(x$rmx$rmxIF$model == "norm"){
    cat(paste(format("Cniper region:", width = 30L, justify = "right"), 
              paste0("(-Inf, ", signif(x$lower, digits = digits), 
                     ") or (", signif(x$upper, digits = digits), ", Inf)")), 
        "\n")
    cat(paste(format("Data in cniper region:", width = 30L, justify = "right"),
              paste0(100*signif(x$prop.lower, digits = digits), " % + ",
                     100*signif(x$prop.upper, digits = digits), " % = ",
                     100*signif(x$prop.lower+x$prop.upper, digits = digits), " %")), 
        "\n")
    cat(paste(format("Prob. of cniper region:", width = 30L, justify = "right"),
              paste0(100*signif(x$p.lower, digits = digits), " % + ",
                     100*signif(x$p.upper, digits = digits), " % = ",
                     100*signif(x$p.lower+x$p.upper, digits = digits), " %")), 
        "\n\n")
  }
  if(is.na(x$rmx$eps)){
    cat(paste(format("Gross-errors for RMX", width = 29L, justify = "right"), 
              paste0(100*signif(x$rmx$eps.lower, digits = digits), " - ",
                     100*signif(x$rmx$eps.upper, digits = digits), " %"), sep = " = "), 
        "\n\n")
  }else{
    cat(paste(format("Gross-errors for RMX", width = 29L, justify = "right"), 
              paste0(100*signif(x$rmx$eps, digits = digits), " %"), sep = " = "), 
        "\n\n")
  }
  invisible(x)
}
plot.cniper <- function(x, add.data = TRUE, color.data = "#0072B5",
                        alpha.data = 0.4, range.alpha = 1e-6, range.n = 501, 
                        color.vline = "#E18727", 
                        ggplot.ggtitle = "Cniper Contamination", 
                        ggplot.xlab = "contamination point", 
                        ggplot.ylab = "asMSE(ML) - asMSE(RMX)", ...){
  if(x$rmx$rmxIF$model == "norm"){
    SD <- x$rmx$rmxEst["sd"]
    tr.invF <- 1.5*SD^2
    maxMSE <- x$rmx$rmxIF$asMSE
    Delta <- sqrt(maxMSE - tr.invF)/x$rmx$rmxIF$radius
    fun <- function(y){
      z <- (y - x$rmx$rmxEst["mean"])/x$rmx$rmxEst["sd"]
      Y <- x$rmx$rmxEst["sd"]*cbind(z, 0.5*(z^2 - 1)) 
      sqrt(as.vector(Y %*% t(Y))) - Delta
    }
    rg <- x$rmx$rmxIF$range(alpha = range.alpha, n = 2)
    y <- c(seq(from = min(rg[1], x$rmx$x), to = max(rg[2], x$rmx$x), 
               length.out = range.n), x$rmx$x)
    y <- sort(unique(y))
    Risk.delta <- sapply(y, fun)
    DF <- data.frame(x = y, y = Risk.delta)
    gg <- ggplot(DF, aes(x = x, y = y)) +
      geom_line() + geom_hline(yintercept = 0) + 
      xlab(ggplot.xlab) + ylab(ggplot.ylab) +
      geom_vline(xintercept = c(x$lower, x$upper), color = color.vline) +
      ggtitle(ggplot.ggtitle)
    if(abs(min(y) - x$rmx$rmxEst["mean"]) >= abs(max(y)- x$rmx$rmxEst["mean"])){
      gg <- gg + 
        annotate(geom = "text", x = min(y), y = 0, hjust = 0, vjust = -1, 
                 label = "RMX better") +
        annotate(geom = "text", x = min(y), y = 0, hjust = 0, vjust = 2, 
                 label = "ML better")
    }else{
      gg <- gg + 
        annotate(geom = "text", x = max(y), y = 0, hjust = 1, vjust = -1, 
                 label = "RMX better") +
        annotate(geom = "text", x = max(y), y = 0, hjust = 1, vjust = 2, 
                 label = "ML better")
    }
    if(add.data){
      DFx <- data.frame(x = x$rmx$x, y = sapply(x$rmx$x, fun))
      gg <- gg + geom_point(data = DFx, aes(x = x, y = y), color = color.data,
                            alpha = alpha.data, inherit.aes = FALSE)
    }
    print(gg)
  }else{
    stop("not yet implemented")
  }
  invisible(gg)
}

