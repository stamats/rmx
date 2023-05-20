\name{cvm}
\alias{cvm}
\alias{cvm.norm}
\alias{cvm.binom}
\alias{cvm.pois}
\alias{cvm.exp}
\alias{rowCVM}
\title{Calculate Cramer von Mises Minimum Distance Estimates}
\description{
  The function \code{cvm} calculates Cramer von Mises Minimum Distance 
  Estimates. It is used to compute valid starting estimates for simple
  parametric models inside of function \code{rmx}. The functions are 
  rarely called directly and are mainly for internal use.
}
\usage{
cvm(x, model = "norm", mu = "model", na.rm = TRUE, \dots)

cvm.norm(x, mu, startPar = NULL)

cvm.binom(x, mu, size)

cvm.pois(x, mu)

cvm.exp(x, mu)

rowCVM(x, model, mu = "model", na.rm = TRUE, parallel = FALSE, ncores = NULL, \dots)
}
\arguments{
  \item{x}{ numeric vector of data. }
  \item{model}{ character: short name of the model/distribution 
                (default = \code{"norm"}); see also details. }
  \item{mu}{ possible options are \code{"model"} and \code{"data"}; see also 
             details. }
  \item{na.rm}{ logical: if \code{TRUE}, \code{NA} values are removed before the 
                estimator is evaluated. }
  \item{startPar}{ parameter values used as starting point for the optimization.
  If missing or \code{NULL}, appropriate default values are computed. }
  \item{size}{ size parameter of binomial distribution (known!); see \code{\link[stats]{dbinom}}. }
  \item{parallel}{ logical: use package \pkg{parallel} for computations. }
  \item{ncores}{ if \code{parallel = TRUE}: number of cores used for computations.
                 If missing, the maximum number of cores - 1 is used. }
  \item{\dots}{ further arguments passed through; e.g., known parameters as 
                for instance \code{size} in case of the binomial model. }
}
\details{
  The function is inspired by function \code{CvMMDEstimator} of package 
  \pkg{distrMod}. We also recommend to use package \pkg{distrMod} for 
  statistical inference based on the Cramer von Mises estimator. In our package
  the function is used as a helper function to provide valid starting estimates
  for our RMX estimators. For more details about minimum distance estimation we 
  refer to Section 6.3.2 of Rieder (1994).
  
  By setting \code{mu = "model"} the model distribution is used as weighting
  measure in the Cramer von Mises distance. Using \code{mu = "data"} leads 
  to weighting by the empiricial distribution of the provided data.
  
  As models we have implemented so far:
  \enumerate{
    \item \code{"norm"}: normal location (\code{mean}) and scale (\code{sd}).
    \item \code{"binom"}: probability of success (\code{prob}).
    \item \code{"pois"}: mean (\code{lambda}).
    \item \code{"exp"}: scale (\code{1/rate}).
  }
}
\value{
  Cramer von Mises minimum distance estimator.
}
\references{
  Kohl, M. (2005). \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.

  M. Kohl, P. Ruckdeschel, and H. Rieder (2010). Infinitesimally Robust Estimation 
  in General Smoothly Parametrized Models. \emph{Statistical Methods and Application}, 
  \bold{19}(3):333-354. 
  
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{rmx}}, \code{\link{mse}}}
\examples{
x <- rnorm(100)
cvm(x, model = "norm")
cvm(x, model = "norm", mu = "data")

x <- rbinom(100, prob = 0.2, size = 3)
cvm(x, model = "binom", size = 3)
cvm(x, model = "binom", size = 3, mu = "data")

x <- rpois(100, lambda = 1.5)
cvm(x, model = "pois")
cvm(x, model = "pois", mu = "data")

x <- rexp(100, rate = 3)
cvm(x, model = "exp")
cvm(x, model = "exp", mu = "data")

###############################################################################
## row-wise computations for data matrices
###############################################################################
X <- matrix(rnorm(200), ncol = 10)
rowCVM(X, model = "norm")

X <- matrix(rbinom(200, prob = 0.2, size = 3), ncol = 10)
rowCVM(X, model = "binom", size = 3)

X <- matrix(rpois(200, lambda = 1.5), ncol = 10)
rowCVM(X, model = "pois")

X <- matrix(rexp(200, rate = 3), ncol = 10)
rowCVM(X, model = "exp")
}
\keyword{robust}
