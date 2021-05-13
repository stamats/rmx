\name{rmx}
\alias{rmx}
\title{Optimally Robust RMX Estimator}
\description{
  The function \code{rmx} computes optimally robust rmx estimators. 
  The definition of these estimators can be found in Kohl (2005) and 
  Rieder et al. (2008), respectively.
}
\usage{
rmx(x, model = "norm", eps.lower=0, eps.upper=0.5, eps=NULL, k = 3L, 
    initial.est=NULL, fsCor = TRUE, na.rm = TRUE)
}
\arguments{
  \item{x}{ numeric vector \code{x} of data values. }
  \item{model}{ character: short name of the model/distribution (default = \code{"norm"});
  see also details. }
  \item{eps.lower}{ positive real (0 <= \code{eps.lower} <= \code{eps.upper}): 
        lower bound for the amount of gross errors. See details below. }
  \item{eps.upper}{ positive real (\code{eps.lower} <= \code{eps.upper} <= 0.5): 
        upper bound for the amount of gross errors. See details below. }
  \item{eps}{ positive real (0 < \code{eps} <= 0.5): amount of gross errors. 
        See details below. }
  \item{k}{ positive integer: k-step is used to compute the optimally robust estimator. }
  \item{initial.est}{ initial estimate for \code{mean} and \code{sd}. If missing 
        median and MAD are used. }
  \item{fsCor}{ logical: perform finite-sample correction; see function \code{\link{fsRadius}}. }
  \item{na.rm}{logical: if \code{TRUE}, \code{NA} values are removed before the estimator is evaluated.}
}
\details{
  If the amount of gross errors (contamination) is known, it can be 
  specified by \code{eps}. The radius of the corresponding infinitesimal 
  contamination neighborhood is obtained by multiplying \code{eps} 
  by the square root of the sample size. 

  If the amount of gross errors (contamination) is unknown, try to find a 
  rough estimate for the amount of gross errors, such that it lies 
  between \code{eps.lower} and \code{eps.upper}.
  
  As models we have implemented so far:
  \enumerate{
    \item \code{"norm"}: normal location and scale. Parameters can be set 
    via argument \code{mean} and \code{sd}; see examples.
  }
}
\value{
  An object of class \code{"optIF"} is returned. It contails at least the 
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
\seealso{\code{\link{optIF}}}
\examples{
ind <- rbinom(100, size=1, prob=0.05) 
x <- rnorm(100, mean=ind*3, sd=(1-ind) + ind*9)

## amount of gross errors unknown
res1 <- rmx(x, eps.lower = 0.01, eps.upper = 0.1)
res1
summary(res1)
confint(res1)
confint(res1, method = "as.bias")

## amount of gross errors known
res2 <- rmx(x, eps = 0.05)
res2
summary(res2)
confint(res2)
confint(res2, method = "as.bias")

## estimator comparison
## classical optimal (non-robust)
c(mean(x), sd(x))
## most robust
c(median(x), mad(x))
## amount of gross errors unknown
res1$rmxEst
## amount of gross errors known
res2$rmxEst
}
\keyword{robust}
