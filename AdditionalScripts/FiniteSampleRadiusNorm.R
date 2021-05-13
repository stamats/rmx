###############################################################################
## Find finite-sample correction factor for asymptotic radius
###############################################################################

library(robustbase)

## source rowRmxNorm.R
## source ZZZinternal.R

## in combination with sysdata.rda of package rmx
rowRMX <- function(x, r, k = 3L){
    mean <- rowMedians(x, na.rm = TRUE)
    sd <- rowMedians(abs(x-mean), na.rm = TRUE)/qnorm(0.75)
    if(r > 10){
        b <- sd*1.618128043
        const <- 1.263094656
        A2 <- b^2*(1+r^2)/(1+const)
        A1 <- const*A2
        a <- -0.6277527697*A2/sd
        mse <- A1 + A2
    }else{
        A1 <- sd^2*.getA1.norm(r)
        A2 <- sd^2*.getA2.norm(r)
        a <- sd*.geta.norm(r)
        b <- sd*.getb.norm(r)
        mse <- A1 + A2
    }
    robEst <- .kstep.norm.matrix(x = x, initial.est = cbind(mean, sd), 
                                 A1 = A1, A2 = A2, a = a, b = b, k = k,
                                 na.rm = TRUE)
    colnames(robEst$est) <- c("mean", "sd")
    robEst$est
}

## find out at which contamination is the maximum finite-sample risk attained
# n <- 10
# M <- 1e5
# eps <- 0.01
# D <- 0.1
# fun <- function(r, x, n){
#     RadMinmax <- rowRMX(x, r = r)
#     n*(mean(RadMinmax[,1]^2) + mean((RadMinmax[,2]-1)^2))
# }
# 
# r <- rbinom(n*M, prob = eps, size = 1)
# Mid <- rnorm(n*M)
# Mcont <- rep(D, n*M)
# Mre <- matrix((1-r)*Mid + r*Mcont, ncol = n)
# ind <- rowSums(matrix(r, ncol = n)) >= n/2
# while(any(ind)){
#     M1 <- sum(ind)
#     cat("M1:\t", M1, "\n")
#     r <- rbinom(n*M1, prob = eps, size = 1)
#     Mid <- rnorm(n*M1)
#     Mcont <- r(contD)(n*M1)
#     Mre[ind,] <- (1-r)*Mid + r*Mcont
#     ind[ind] <- rowSums(matrix(r, ncol = n)) >= n/2
# }
# 
# fun1 <- function(D){
#     Mcont <- rep(D, n*M)
#     Mre <- matrix((1-r)*Mid + r*Mcont, ncol = n)
#     fun(r = 1, x = Mre, n = n)
# }
# sapply(c(seq(0.1, 10, length = 20), 20, 50, 100, 1000, 1e4, 1e6), fun1)


## finite-sample optimal radius
## n at least 3, for n = 2 not possible to have less than 50% contamination
n <- c(3:50, seq(55, 100, by = 5), seq(110, 200, by = 10), seq(250, 500, by = 50))
eps <- c(seq(0.001, 0.01, by = 0.001), seq(0.02, to = 0.5, by = 0.01))
M <- 1e5
D <- 1e6

r.fi <- matrix(NA, nrow = length(eps), ncol = length(n))
colnames(r.fi) <- n
rownames(r.fi) <- eps
r.as <- r.fi
for(j in seq(along = n)){
    ptm <- proc.time()
    cat("aktuelles n:\t", n[j], "\n")
    i <- 0
    repeat{
        i <- i + 1
        cat("aktuelles eps:\t", eps[i], "\n")
        r <- rbinom(n[j]*M, prob = eps[i], size = 1)
        Mid <- rnorm(n[j]*M)
        Mcont <- rep(D, n[j]*M)
        Mre <- matrix((1-r)*Mid + r*Mcont, ncol = n[j])
        rm(Mid, Mcont)
        gc()
        ind <- rowSums(matrix(r, ncol = n[j])) >= n[j]/2
        rm(r)
        gc()
        while(any(ind)){
            M1 <- sum(ind)
            cat("M1:\t", M1, "\n")
            r <- rbinom(n[j]*M1, prob = eps[i], size = 1)
            Mid <- rnorm(n[j]*M1)
            Mcont <- rep(D, n[j]*M1)
            Mre[ind,] <- (1-r)*Mid + r*Mcont
            ind[ind] <- rowSums(matrix(r, ncol = n[j])) >= n[j]/2
            rm(Mid, Mcont, r)
            gc()
        }
        fun <- function(r, x, n){
            RadMinmax <- rowRMX(x, r = r)
            n*(mean(RadMinmax[,1]^2) + mean((RadMinmax[,2]-1)^2))
        }
        r.fi[i,j] <- optimize(fun, interval = c(eps[i], min(max(2, n[j]*eps[i]*25), 10)), x = Mre, n = n[j])$minimum
        r.as[i,j] <- sqrt(n[j])*eps[i]
        cat("finit:\t", r.fi[i,j], "\t asympt:\t", r.as[i,j], "\n")
        rm(Mre)
        gc()
        if(round(r.fi[i,j], 2) == 1.74 | i == length(eps)) break
    }
    save.image(file = "FiniteSampleRadiusNorm.RData")
    cat("Dauer:\t", proc.time() - ptm, "\n")
}

r.as <- outer(eps, sqrt(n))

## keine Konvergenz der Optimierung bzw. "falsche" Optima
sum(r.fi > 1.744, na.rm = TRUE)
r.fi[r.fi > 1.744 & !is.na(r.fi)]

r.finit <- r.fi
## setze auf NA
r.finit[r.fi > 1.744 & !is.na(r.fi)] <- NA
r.finit[r.finit > 1.744 & !is.na(r.finit)]
sum(r.finit > r.as & !is.na(r.finit), na.rm = TRUE)
sum(!is.na(r.finit))

## r.finit muss größer oder gleich r.as sein
r.finit[r.finit <= r.as & !is.na(r.finit)]
r.as[r.finit <= r.as & !is.na(r.finit)]
r.finit[r.finit <= r.as & !is.na(r.finit)] <- r.as[r.finit <= r.as & !is.na(r.finit)]
r.finit[r.finit < r.as & !is.na(r.finit)]

for(i in 1:nrow(r.finit)){
    cat("====================================================================\n")
    cat(eps[i], "\n")
    ind <- r.finit[i,] > 1.744
    fit <- lm(r.finit[i,!ind] ~ I(1/sqrt(n[!ind])) + I(1/n[!ind]) + I(1/n[!ind]^(3/2)))
    print(summary(fit))
#    fit <- lm(r.finit[i,!ind] ~ I(1/n[!ind]^(3/2)))
#    print(summary(fit))
}

## ab ca. 1.655 keine Änderung mehr sichtbar
r.finit[18,]
r.finit2 <- r.finit
## minimales b, ab ca. r = 1.67
r.finit2[is.na(r.finit) | r.finit > 1.655] <- 1.67

plot(n, r.finit[19,])
plot(n, r.finit2[19,])

for(i in c(1:11)){
    r.finit2[i, r.finit2[i,] > 1.655] <- approx(n[r.finit2[i,] <= 1.655], 
                                                r.finit2[i,r.finit2[i,] <= 1.655], 
                                                n[r.finit2[i,] > 1.655], rule = 2)$y
}
r.finit2[12,2] <- 1.67

r.finit2 <- round(r.finit2, 4)

#load("sysdata.rda")
.fsRadius.norm <- r.finit2
save(.fsRadius.norm, .radius.gitter.norm, .A1.norm, .A2.norm, 
     .a.norm, .b.norm, file = "sysdata.rda")
