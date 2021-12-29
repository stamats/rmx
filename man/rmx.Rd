\name{rmx}
\alias{rmx}
\title{Optimally Robust RMX Estimator}
\description{
  The function \code{rmx} computes optimally robust rmx estimators. 
  The definition of these estimators can be found in Kohl (2005) and 
  Rieder et al. (2008), respectively.
}
\usage{
rmx(x, model = "norm", eps.lower=0, eps.upper=NULL, eps=NULL, k = 3L, 
    initial.est=NULL, fsCor = TRUE, na.rm = TRUE, message = TRUE, \dots)
}
\arguments{
  \item{x}{ numeric vector \code{x} of data values. }
  \item{model}{ character: short name of the model/distribution 
                (default = \code{"norm"}); see also details. }
  \item{eps.lower}{ positive real (0 <= \code{eps.lower} <= \code{eps.upper}): 
        lower bound for the amount of gross errors; see details below. }
  \item{eps.upper}{ positive real (\code{eps.lower} <= \code{eps.upper} <= 0.5): 
        upper bound for the amount of gross errors; see details below. }
  \item{eps}{ positive real (0 < \code{eps} <= 0.5): amount of gross errors. 
        See details below. }
  \item{k}{ positive integer: k-step is used to compute the optimally robust estimator. }
  \item{initial.est}{ initial estimate for \code{mean} and \code{sd}. If missing 
        median and MAD are used. }
  \item{fsCor}{ logical: perform finite-sample correction; see function \code{\link{fsRadius}}. }
  \item{na.rm}{ logical: if \code{TRUE}, \code{NA} values are removed before the estimator is evaluated.}
  \item{message}{ logical: if \code{FALSE}, messages are suppressed.}
  \item{\dots}{ further arguments passed through; e.g., known parameters as 
                for instance \code{size} in case of the binomial model. }
}
\details{
  If the amount of gross errors (contamination) is known, it can be 
  specified by \code{eps}. The radius of the corresponding infinitesimal 
  contamination neighborhood is obtained by multiplying \code{eps} 
  by the square root of the sample size. 

  If the amount of gross errors (contamination) is unknown, try to find a 
  rough estimate for the amount of gross errors, such that it lies 
  between \code{eps.lower} and \code{eps.upper}.
  
  If neither \code{eps} nor \code{eps.upper} is provided, \code{eps.upper}
  will be estimated by applying function \code{outlier} to the RMX estimator
  with \code{eps.lower = 0} and \code{eps.upper = 0.5}.
  
  As models we have implemented so far:
  \enumerate{
    \item \code{"norm"}: normal location (mean) and scale (sd).
  }
}
\value{
  An object of class \code{"rmx"} is returned. It contails at least the 
  following arguments:
  \item{rmxEst}{ estimates }
  \item{rmxIF}{ object of class \code{optIF}; see \code{\link{optIF}}. }
  \item{initial.est}{ initial estimates. }
  \item{Infos}{ matrix with information about the estimator }
  \item{x}{ data used for the estimation. }
  \item{n}{ sample size }
  \item{eps.lower}{ lower bound for the amount of gross errors, if provided 
                    otherwise \code{NA}. }
  \item{eps.upper}{ upper bound for the amount of gross errors, if provided 
                    otherwise \code{NA}. }
  \item{eps}{ amount of gross errors, if provided otherwise \code{NA}. }
  \item{fsCor}{ finite-sample correction }
  \item{k}{ k-step construction }
  \item{call}{ matched call }
}
\references{
  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Rieder, H., Kohl, M. and Ruckdeschel, P. (2008) The Costs of not Knowing
  the Radius. Statistical Methods and Applications \emph{17}(1) 13-40.
  Extended version: \url{http://r-kurs.de/RRlong.pdf}
  
  M. Kohl, P. Ruckdeschel, and H. Rieder (2010). Infinitesimally Robust Estimation 
  in General Smoothly Parametrized Models. \emph{Statistical Methods and Application}, 
  \bold{19}(3):333-354. 
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{optIF}}, \code{\link{rowRmx}}, \code{\link{fsRadius}}}
\examples{
###########################################################
## generate some contaminated random data
###########################################################
ind <- rbinom(100, size=1, prob=0.05) 
x <- rnorm(100, mean=ind*3, sd=(1-ind) + ind*9)

###########################################################
## amount of gross errors in [1, 10]%
###########################################################
res1 <- rmx(x, eps.lower = 0.01, eps.upper = 0.1)
res1
coef(res1)
vcov(res1)
sqrt(diag(vcov(res1)))

summary(res1)

confint(res1) # method = "as"
confint(res1, method = "as.bias")
confint(res1, method = "boot", R = 999)
confint(res1, method = "boot", R = 999, parallel = TRUE, ncores = 2)

(cni.res1 <- cniper(res1))
getCnipers(cni.res1)
## or
getCnipers(res1)
plot(cni.res1)

(out.res1 <- outlier(res1))
getOutliers(out.res1)
## or
getOutliers(res1)

ifPlot(res1)
plot(res1, which = 1)

## only influence function
plot(res1$rmxIF)

qqPlot(res1)
plot(res1, which = 2)

ppPlot(res1)
plot(res1, which = 3)

dPlot(res1)
plot(res1, which = 4)

aiPlot(res1)
plot(res1, which = 5)

riPlot(res1)
plot(res1, which = 6)

iiPlot(res1)
plot(res1, which = 7)

\donttest{
## to save check time
###########################################################
## amount of gross errors known
###########################################################
res2 <- rmx(x, eps = 0.05)
res2
coef(res2)
vcov(res2)
sqrt(diag(vcov(res2)))

summary(res2)

confint(res2) # method = "as"
confint(res2, method = "as.bias")
confint(res2, method = "boot", R = 999)
confint(res2, method = "boot", R = 999, parallel = TRUE, ncores = 2)

(cni.res2 <- cniper(res2))
getCnipers(cni.res2)
## or
getCnipers(res2)
plot(cni.res2)

(out.res2 <- outlier(res2))
getOutliers(out.res2)
## or
getOutliers(res2)

ifPlot(res2)
plot(res2, which = 1)

## only influence function
plot(res2$rmxIF)

qqPlot(res2)
plot(res2, which = 2)

ppPlot(res2)
plot(res2, which = 3)

dPlot(res2)
plot(res2, which = 4)

aiPlot(res1)
plot(res1, which = 5)

riPlot(res1)
plot(res1, which = 6)

iiPlot(res1)
plot(res1, which = 7)

###########################################################
## estimator comparison
###########################################################
## classical optimal (non-robust)
c(mean(x), sd(x))
## most robust
c(median(x), mad(x))
## amount of gross errors unknown
res1$rmxEst
## amount of gross errors known
res2$rmxEst
}
}
\keyword{robust}
