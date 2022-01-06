\name{rowRmx and colRmx}
\alias{rowRmx}
\alias{colRmx}
\alias{rowRmx.norm}
\alias{rowRmx.binom}
\title{Optimally robust estimation for location and/or scale}
\description{
  The functions \code{rowRmx} and \code{colRmx} compute 
  optimally robust RMX estimates for (convex) contamination neighborhoods. 
  The definition of these estimators can be found in Kohl (2005) and 
  Rieder et al. (2008), respectively.
}
\usage{
rowRmx(x, model = "norm", eps.lower=0, eps.upper=0.5, eps=NULL, k = 3L, 
       initial.est=NULL, fsCor = NULL, na.rm = TRUE, message = TRUE, 
       computeSE = NULL, \dots)

colRmx(x, model = "norm", eps.lower=0, eps.upper=0.5, eps=NULL, k = 3L, 
       initial.est=NULL, fsCor = NULL, na.rm = TRUE, message = TRUE, 
       computeSE = NULL, \dots)

rowRmx.norm(x, eps.lower=0, eps.upper, eps = NULL, initial.est = NULL, k = 3L, 
            fsCor = TRUE, na.rm = TRUE, computeSE = TRUE)
            
rowRmx.binom(x, eps.lower=0, eps.upper, eps = NULL, initial.est = NULL, k = 3L, 
             fsCor = FALSE, na.rm = TRUE, size, computeSE = FALSE, 
             parallel = FALSE, ncores = NULL, aUp = 100*size, 
             cUp = 1e4, delta = 1e-9)
}
\arguments{
  \item{x}{ matrix or data.frame of (numeric) data values. }
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
  \item{na.rm}{logical: if \code{TRUE}, \code{NA} values are removed before the estimator is evaluated.}
  \item{message}{logical: if \code{FALSE}, messages are suppressed.}
  \item{size}{ size parameter (known!); see \code{\link[stats]{dbinom}}. }
  \item{computeSE}{ logical: compute asymptotic standard errors. }
  \item{parallel}{ if \code{computeSE = TRUE}: logical: use package parallel for the computation. }
  \item{ncores}{ if \code{parallel = TRUE}: number of cores used for the computation.
                 If missing, the maximum number of cores - 1 is used. }
  \item{aUp}{ numeric: upper limit for centering constant a. }
  \item{cUp}{ postive real: upper limit for clipping constant c. }
  \item{delta}{ positive real: desired accuracy (convergence tolerance). }
  \item{\dots}{ further arguments passed through; e.g., known parameters such as 
                \code{size} in case of the binomial model. }
}
\details{
  These functions are optimized for the situation where one has a matrix 
  and wants to compute the optimally robust RMX estimator for every row, 
  respectively column of this matrix. In particular, the amount of cross
  errors is assumed to be constant for all rows, respectively columns.

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
  An object of class \code{"RMXlist"} is returned. It contails at least the 
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
\seealso{\code{\link{rmx}}, \code{\link{optIF}}, \code{\link{fsRadius}}}
\examples{
ind <- rbinom(500, size=1, prob=0.05)
X <- matrix(rnorm(500, mean=ind*3, sd=(1-ind) + ind*9), nrow = 5)
rowRmx(X)
rowRmx(X, eps.lower = 0.01, eps.upper = 0.1)
rowRmx(X, eps.lower = 0.01, eps.upper = 0.1, computeSE = TRUE)
rowRmx(X, eps = 0.05)

X1 <- t(X)
colRmx(X1)
colRmx(X1, eps.lower = 0.01, eps.upper = 0.1)
colRmx(X1, eps.lower = 0.01, eps.upper = 0.1, computeSE = TRUE)
colRmx(X1, eps = 0.05)

## adaptive determination of outlier contamination
RMXlist <- apply(X, 1, rmx, eps.upper = 0.5)
eps.up <- max(sapply(RMXlist, function(x) outlier(x)$prop.outlier))
rowRmx(X, eps.upper = eps.up)
rowRmx(X, eps.upper = eps.up, computeSE = TRUE)
}
\keyword{robust}
