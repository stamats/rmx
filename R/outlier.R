outlier <- function(x, ...){
  UseMethod("outlier")
}
outlier.rmx <- function(x, prob = 0.001){
  stopifnot(length(prob) == 1)
  stopifnot(is.numeric(prob))
  stopifnot((prob > 0) && (prob < 0.5))
  
  if(x$rmxIF$model %in% c("norm", "gamma")){
    pt.out <- x$rmxIF$range(alpha = prob, n = 2)
    prop.lo <- sum(x$x < pt.out[1])/x$n
    prop.up <- sum(x$x > pt.out[2])/x$n
    p.lower <- p.upper <- prob
  }
  if(x$rmxIF$model %in% c("binom", "pois")){
    pt.out <- x$rmxIF$range(alpha = prob)
    pt.out <- pt.out[c(1, length(pt.out))]
    prop.lo <- sum(x$x <= pt.out[1])/x$n
    prop.up <- sum(x$x >= pt.out[2])/x$n
    p.lower <- pbinom(pt.out[1], size = x$rmxIF$parameter["size (known)"],
                      prob = x$rmxEst)
    p.upper <- pbinom(pt.out[2]-1, size = x$rmxIF$parameter["size (known)"],
                      prob = x$rmxEst, lower.tail = FALSE)
  }
  
  res <- list(rmx = x,
              lower = pt.out[1], upper = pt.out[2], 
              prop.outlier = prop.lo+prop.up, 
              p.outlier = p.lower+p.upper, prop.lower = prop.lo, 
              prop.upper = prop.up, p.lower = p.lower, p.upper = p.upper)
  class(res) <- "outlier"
  res
}
getOutliers <- function(x, ...){
  if(!inherits(x, "outlier")) x <- outlier(x, ...)
  out <- c(which(x$rmx$x <= x$lower), which(x$rmx$x >= x$upper))
  val <- x$rmx$x[out]
  ind <- order(val)
  list(values = val[ind], indices = out[ind])
}
print.outlier <- function(x, digits = 3){
  cat("\n")
  cat("     ", strwrap(paste0("Outlier region for ", x$rmx$rmxIF$modelName), 
                       prefix = " "), "\n")
  cat("\n")
  cat(paste(format("Parameter:", width = 30L, justify = "right"), 
            paste0(paste(names(x$rmx$rmxIF$parameter), 
                         signif(x$rmx$rmxIF$parameter, digits = digits), 
                         sep = " = "), collapse = ", ")), "\n")
  if(x$rmx$rmxIF$model == "norm"){
    cat(paste(format("Outlier region:", width = 30L, justify = "right"), 
              paste0("(-Inf, ", signif(x$lower, digits = digits), 
                     ") or (", signif(x$upper, digits = digits), ", Inf)")), 
        "\n")
    cat(paste(format("Data in oulier region:", width = 30L, justify = "right"),
              paste0(100*signif(x$prop.lower, digits = digits), " % + ",
                     100*signif(x$prop.upper, digits = digits), " % = ",
                     100*signif(x$prop.lower+x$prop.upper, digits = digits), " %")), 
        "\n")
    cat(paste(format("Prob. of outlier region:", width = 30L, justify = "right"),
              paste0(100*signif(x$p.lower, digits = digits), " % + ",
                     100*signif(x$p.upper, digits = digits), " % = ",
                     100*signif(x$p.lower+x$p.upper, digits = digits), " %")), 
        "\n\n")
  }
  if(x$rmx$rmxIF$model == "binom"){
    supp <- seq(from = 0, to = x$rmx$rmxIF$parameter["size (known)"], by = 1)
    outlier.lo <- supp[supp <= x$lower] 
    outlier.up <- supp[supp >= x$upper]
    if(length(outlier.lo) > 5){
      text.lo <- paste0("{", paste0(outlier.lo[1:3], collapse = ", "), 
                        ", ..., ", outlier.lo[length(outlier.lo)], "}")
    }
    else{
      if(length(outlier.lo) > 0){
        text.lo <- paste0("{", paste0(outlier.lo, collapse = ", "), "}")
      }else{
        text.lo <- NULL
      }
    }
    
    if(length(outlier.up) > 5){
      text.up <- paste0("{", paste0(outlier.up[1:3], collapse = ", "), 
                        ", ..., ", outlier.up[length(outlier.up)], "}")
    }else{
      if(length(outlier.up) > 0){
        text.up <- paste0("{", paste0(outlier.up, collapse = ", "), "}")
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
    cat(paste(format("Outlier region:", width = 30L, justify = "right"), 
              text.region), 
        "\n")
    cat(paste(format("Data in outlier region:", width = 30L, justify = "right"),
              paste0(100*signif(x$prop.lower, digits = digits), " % + ",
                     100*signif(x$prop.upper, digits = digits), " % = ",
                     100*signif(x$prop.lower+x$prop.upper, digits = digits), " %")), 
        "\n")
    cat(paste(format("Prob. of outlier region:", width = 30L, justify = "right"),
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
