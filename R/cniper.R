cniper <- function(x, ...){
  UseMethod("cniper")
}
cniper.rmx <- function(x, range.alpha = 1e-6){
  stopifnot(length(range.alpha) == 1)
  stopifnot(is.numeric(range.alpha))
  stopifnot((range.alpha > 0) && (range.alpha < 0.5))
  if(x$rmxIF$model == "norm"){
    res <- .cniper.norm(x = x, range.alpha = range.alpha)
  }
  if(x$rmxIF$model == "binom"){
    res <- .cniper.binom(x = x)
  }
  if(!x$rmxIF$model %in% c("norm", "binom"))
    stop("'cniper' not yet implemented for given model.")
  res
}
.cniper.norm <- function(x, range.alpha){
  SD <- x$rmxEst["sd"]
  tr.invF <- 1.5*SD^2
  maxMSE <- x$rmxIF$asMSE
  r <- x$rmxIF$radius
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
  res
}
.cniper.binom <- function(x){
  prob <- x$rmxEst["prob"]
  size <- x$rmxIF$parameter["size (known)"]
  tr.invF <- (prob*(1-prob))/size
  maxMSE <- x$rmxIF$asMSE
  
  Delta <- sqrt(maxMSE - tr.invF)/x$rmxIF$radius
  supp <- seq(from = 0, to = size, by = 1)
  Diff <- abs(supp/size - prob) - Delta
  if(all(Diff < 0)){
    pt.lo <- -1
    pt.up <- size + 1
    names(pt.up) <- NULL
    prop.lo <- 0
    prop.up <- 0
    p.lo <- 0
    p.up <- 0
  }else{
    M <- trunc(size*prob)
    if(any(Diff[supp >= M] > 0)){
      ind.up <- min(which(Diff[supp >= M] > 0))
      pt.up <- supp[supp >= M][ind.up]
    }else{
      pt.up <- size + 1
    }
    if(any(Diff[supp < M] > 0)){
      ind.lo <- max(which(Diff[supp < M] > 0))
      pt.lo <- supp[supp < M][ind.lo]
    }else{
      pt.lo <- -1
    }
  }
  
  prop.lo <- sum(x$x <= pt.lo)/x$n
  prop.up <- sum(x$x >= pt.up)/x$n
  p.lo <- pbinom(pt.lo, prob = prob, size = size)
  p.up <- pbinom(pt.up-1, prob = prob, size = size, lower.tail = FALSE)
  res <- list(rmx = x,
              lower = pt.lo, upper = pt.up, prop.cniper = prop.lo+prop.up, 
              p.cniper = p.lo+p.up, prop.lower = prop.lo, prop.upper = prop.up, 
              p.lower = p.lo, p.upper = p.up)
  class(res) <- "cniper"
  res
}
getCnipers <- function(x, ...){
  if(!inherits(x, "cniper")) x <- cniper(x, ...)
  out <- c(which(x$rmx$x <= x$lower), which(x$rmx$x >= x$upper))
  val <- x$rmx$x[out]
  ind <- order(val)
  list(values = val[ind], indices = out[ind])
}
print.cniper <- function(x, digits = 3){
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
  if(x$rmx$rmxIF$model == "binom"){
    supp <- seq(from = 0, to = x$rmx$rmxIF$parameter["size (known)"], by = 1)
    cniper.lo <- supp[supp <= x$lower] 
    cniper.up <- supp[supp >= x$upper]
    if(length(cniper.lo) > 5){
      text.lo <- paste0("{", paste0(cniper.lo[1:3], collapse = ", "), 
                        ", ..., ", cniper.lo[length(cniper.lo)], "}")
    }
    else{
      if(length(cniper.lo) > 0){
        text.lo <- paste0("{", paste0(cniper.lo, collapse = ", "), "}")
      }else{
        text.lo <- NULL
      }
    }
      
    if(length(cniper.up) > 5){
      text.up <- paste0("{", paste0(cniper.up[1:3], collapse = ", "), 
                        ", ..., ", cniper.up[length(cniper.up)], "}")
    }else{
      if(length(cniper.up) > 0){
        text.up <- paste0("{", paste0(cniper.up, collapse = ", "), "}")
      }else{
        text.up <- NULL
      }
    }
    if(is.null(text.lo) && is.null(text.up)){
      text.region <- "{}"
    }else{
      if(is.null(text.lo)){
        text.region <- text.up
      }
      if(is.null(text.up)){
        text.region <- text.lo
      }
      if(!is.null(text.lo) && !is.null(text.up)){
        text.region <- paste0(text.lo, " or ", text.up) 
      }
    }
    cat(paste(format("Cniper region:", width = 30L, justify = "right"), 
                text.region), 
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
                        ggplot.ylab = "asMSE(ML) - asMSE(RMX)"){
  if(x$rmx$rmxIF$model == "norm"){
    SD <- x$rmx$rmxEst["sd"]
    tr.invF <- 1.5*SD^2
    maxMSE <- x$rmx$rmxIF$asMSE
    r <- x$rmx$rmxIF$radius
    fun <- function(y){
      z <- (y - x$rmx$rmxEst["mean"])/x$rmx$rmxEst["sd"]
      Y <- x$rmx$rmxEst["sd"]*cbind(z, 0.5*(z^2 - 1)) 
      tr.invF + r^2*as.vector(Y %*% t(Y)) - maxMSE
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
  }
  if(x$rmx$rmxIF$model == "binom"){
    prob <- x$rmx$rmxEst["prob"]
    size <- x$rmx$rmxIF$parameter["size (known)"]
    tr.invF <- (prob*(1-prob))/size
    maxMSE <- x$rmx$rmxIF$asMSE
    r <- x$rmx$rmxIF$radius
    
    supp <- seq(from = 0, to = size, by = 1)
    Risk.delta <- tr.invF + r^2*(supp/size - prob)^2 - maxMSE
    DF <- data.frame(x = supp, y = Risk.delta)
    gg <- ggplot(DF, aes(x = x, y = y)) +
      geom_line() + geom_hline(yintercept = 0) + 
      xlab(ggplot.xlab) + ylab(ggplot.ylab) +
      geom_vline(xintercept = c(x$lower, x$upper), color = color.vline) +
      ggtitle(ggplot.ggtitle)
    if(x$rmx$rmxEst["prob"] >= 0.5){
      gg <- gg + 
        annotate(geom = "text", x = -2, y = 0, hjust = 0, vjust = -1, 
                 label = "RMX better") +
        annotate(geom = "text", x = -2, y = 0, hjust = 0, vjust = 2, 
                 label = "ML better")
    }else{
      gg <- gg + 
        annotate(geom = "text", x = size+2, y = 0, hjust = 1, vjust = -1, 
                 label = "RMX better") +
        annotate(geom = "text", x = size+2, y = 0, hjust = 1, vjust = 2, 
                 label = "ML better")
    }
    if(add.data){
      DFx <- data.frame(x = x$rmx$x, y = tr.invF + r^2*(x$rmx$x/size - prob)^2 - maxMSE)
      gg <- gg + geom_point(data = DFx, aes(x = x, y = y), color = color.data,
                            alpha = alpha.data, inherit.aes = FALSE)
    }
  }
  if(!x$rmx$rmxIF$model %in% c("norm", "binom"))
    stop("'plot' not yet implemented for given model.")
  gg
}
