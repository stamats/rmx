\name{cvm}
\alias{cvm}
\title{Calculate Cramer von Mises Minimum Distance Estimates}
\description{
  The function \code{cvm} calculates Cramer von Mises Minimum Distance 
  Estimates. It is used to compute valid starting estimates for simple
  parametric models inside of function \code{rmx}. 
}
\usage{
cvm(x, model = "norm", mu = "model", na.rm = TRUE, startPar = NULL, \dots)
}
\arguments{
  \item{x}{ numeric vector of data. }
  \item{model}{ character: short name of the model/distribution 
                (default = \code{"norm"}); see also details. }
  \item{mu}{ possible options are \code{"model"} and \code{"data"}; see also 
             details. }
  \item{na.rm}{ logical: if \code{TRUE}, \code{NA} values are removed before the 
                estimator is evaluated. }
  \item{startPar}{ optional numeric vector with parameter values that should be
                   used as starting values in the numerical optimization via 
                   \code{optim}. The parameter is only of relevance, if a model
                   includes two or more unknown parameters. If it is not 
                   specified, starting values are computed inside the function. }
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
}
\keyword{robust}
