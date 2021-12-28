###############################################################################
## Evaluate rmx on rows of a matrix
###############################################################################
rowRmx <- function(x, model = "norm", eps.lower=0, eps.upper=0.5, eps=NULL, 
                   k = 3L, initial.est=NULL, fsCor = TRUE, na.rm = TRUE, 
                   computeSE = FALSE, parallel = FALSE, ncores = NULL,
                   message = TRUE){
    es.call <- match.call()
    if(missing(x))
        stop("'x' is missing with no default")
    if(is.data.frame(x))
        x <- data.matrix(x)
    else
        x <- as.matrix(x)
    if(!is.matrix(x))
        stop("'x' has to be a matrix resp. convertable to a matrix by 'as.matrix'
              or 'data.matrix'")
    
    stopifnot(is.character(model))
    stopifnot(length(model) == 1)
    if(!is.na(pmatch(model, "norm")))  model <- "norm"
    MODELS <- c("norm", "binom", "pois", "gamma")
    model <- pmatch(model, MODELS)
    
    if(is.null(eps)){
        if(length(eps.lower) != 1 || length(eps.upper) != 1)
            stop("'eps.lower' and 'eps.upper' have to be of length 1")
        if(!is.numeric(eps.lower) || !is.numeric(eps.upper) || eps.lower >= eps.upper) 
            stop("'eps.lower' < 'eps.upper' is not fulfilled")
        if((eps.lower < 0) || (eps.upper > 0.5))
            stop("'eps.lower' and 'eps.upper' have to be in [0, 0.5]")
    }else{
        if(length(eps) != 1)
            stop("'eps' has to be of length 1")
        if((eps < 0) || (eps > 0.5))
            stop("'eps' has to be in (0, 0.5]")
        if(eps == 0 && message){
            message("'eps=0': Maximum likelihood estimator is computed.")
        }
    }
    
    stopifnot(is.numeric(k))
    if(!is.integer(k))  k <- as.integer(k)
    if(k < 1){
        stop("'k' has to be some positive integer value")
    }
    stopifnot(length(k) == 1)
    stopifnot(length(fsCor) == 1)
    stopifnot(is.logical(fsCor))
    stopifnot(length(na.rm) == 1)
    stopifnot(is.logical(na.rm))
    
    if(model == 1){ # normal distribution
        RMXest <- rowRmx.norm(x, eps.lower = eps.lower, eps.upper = eps.upper, 
                              eps = eps, k = k, initial.est=initial.est, 
                              fsCor = fsCor, na.rm = na.rm, computeSE = computeSE, 
                              parallel = parallel, ncores = ncores)
    }
    if(model != 1){
        stop("Given 'model' not yet implemented")
    }
    completecases <- rowSums(!is.na(x))
    if(any(completecases < ncol(x))){
        RMXest$Infos <- rbind(RMXest$Infos, 
                              c("rowRmx", 
                                "'NA' values have been removed before the computation"))
    }

    RMXest$n <- ncol(x)
    RMXest$eps.lower <- ifelse(is.null(eps), eps.lower, NA)
    RMXest$eps.upper <- ifelse(is.null(eps), eps.upper, NA)
    RMXest$eps <- ifelse(is.null(eps), NA, eps)
    RMXest$fsCor <- fsCor
    RMXest$k <- k
    RMXest$call <- es.call
    RMXest
}
print.RMXlist <- function (x, digits = getOption("digits"), prefix = " ", head.n = 6L, ...){
    cat("\n")
    cat(strwrap(paste0("RMX estimator for ", x$modelName), 
                prefix = " "), sep = "\n")
    cat("\n")
    cat(paste(format("Sample size", width = 24L, justify = "right"), 
              format(x$n, digits = digits), sep = " = "), sep = "\n")
    if(is.na(x$eps)){
        cat(paste(format("Amount gross-errors", width = 24L, justify = "right"), 
                  paste0(100*signif(x$eps.lower, digits = 3), " - ",
                         100*signif(x$eps.upper, digits = 3), " %"), sep = " = "), 
            sep = "\n")
    }else{
        cat(paste(format("Amount gross-errors", width = 24L, justify = "right"), 
                  paste0(100*signif(x$eps, digits = 3), " %"), sep = " = "), 
            sep = "\n")
    }
    if(x$fsCor){
        cat(paste(format("FS-corrected radius", width = 24L, justify = "right"), 
                  format(x$radius, digits = digits), sep = " = "), sep = "\n")
    }else{
        cat(paste(format("Infinitesimal radius", width = 24L, justify = "right"), 
                  format(x$radius, digits = digits), sep = " = "), sep = "\n")
    }
    cat("\n Estimates:\n")
    print(head(x$rmxEst, n = head.n), digits = digits, ...)
    if(head.n > nrow(x$rmxEst)) cat(paste0("[", head.n+1, ",]  ...", collapse = ""), "\n")
    if(any(!is.na(x$asSE))) {
        cat("\n Asymptotic standard errors:\n")
        print(head(x$asSE, n = head.n), digits = digits, ...)
    }
    if(head.n < nrow(x$rmxEst)) cat(paste0("[", head.n+1, ",]  ...", collapse = ""), "\n")
    if(nrow(x$rmxEst) > head.n){
        cat("\n NOTE:", nrow(x$rmxEst)-head.n, "rows omitted\n")
    }
    cat("\n Call:\n ", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    invisible(x)
}
