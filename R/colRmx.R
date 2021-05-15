###############################################################################
## Evaluate rmx on columns of a matrix
###############################################################################
colRmx <- function(x, model = "norm", eps.lower=0, eps.upper=0.5, eps=NULL, 
                   k = 3L, initial.est=NULL, fsCor = TRUE, na.rm = TRUE, 
                   computeSE = FALSE, parallel = FALSE, ncores = NULL){
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

    res <- rowRmx(x = t(x), model = model, 
                  eps.lower = eps.lower, eps.upper = eps.upper, eps = eps, 
                  k = k, initial.est = initial.est, fsCor = fsCor, 
                  na.rm = na.rm, computeSE = computeSE, parallel = parallel,
                  ncores = ncores)
    res$call <- es.call
    res
}
