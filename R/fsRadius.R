###############################################################################
## Function for finite-sample correction of the neighborhood radius
###############################################################################
fsRadius <- function(r, n, model = "norm"){
    stopifnot(is.numeric(r))
    stopifnot(length(r) == 1)
    stopifnot(r >= 0)
    if(r == 0) return(r)
    stopifnot(is.numeric(n))
    n <- as.integer(n)
    stopifnot(n > 0)
    if(n == 1) return(Inf)
    if(n == 2) return(Inf)
    
    if(model == "norm"){
        if(model == "norm" & r >= 1.67) return(r)
        
        eps <- r/sqrt(n)
        ns <- c(3:50, seq(55, 100, by = 5), seq(110, 200, by = 10), 
                seq(250, 500, by = 50))
        ns <- as.integer(ns)
        epss <- c(seq(0.001, 0.01, by = 0.001), seq(0.02, to = 0.5, by = 0.01))
        if(n %in% ns){
            ind <- ns == n
            r <- max(r, approx(x = epss, y = .fsRadius.norm[,ind], 
                               xout = eps, rule = 2)$y)
            return(r)
        }else{
            if(n > 500){
                ind <- ns == 500L
                r1 <- approx(x = epss, y = .fsRadius.norm[,ind], 
                             xout = eps, rule = 2)$y
                r <- max(r, mean(c(r1, r)))
                return(r)
            }else{
                ind <- sort(order(abs(ns-n))[1:2])
                r1 <- approx(x = epss, y = .fsRadius.norm[,ind[1]], 
                             xout = eps, rule = 2)$y
                r2 <- approx(x = epss, y = .fsRadius.norm[,ind[2]], 
                             xout = eps, rule = 2)$y
                D <- ns[ind[2]] - ns[ind[1]]
                r <- max(r, (n-ns[ind[1]])/D*r1 + (ns[ind[2]]-n)/D*r2)
                return(r)
            }
        }
    }else{
        stop("argument 'model' has to be 'norm'")
    }
}