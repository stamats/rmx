\name{bias}
\alias{bias}
\alias{bias.rmx}
\title{Calculate Bias for a Fitted Model Object}
\description{
  The function \code{bias} extracts or calculates bias of a fitted model 
  object.
}
\usage{
bias(object, \dots)

\method{bias}{rmx}(object, \dots)
}
\arguments{
  \item{object}{ a fitted model object, especially an object of S3 class \code{rmx}. }
  \item{\dots}{ further arguments passed through. }
}
\details{
  The function is inspired by function \code{vcov}.
  
  In case of optimally-robust RMX estimators the maximum asymptotic bias 
  is extracted and returned.
}
\value{
  Numeric vector with the bias.
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
ind <- rbinom(100, size=1, prob=0.05) 
x <- rnorm(100, mean=ind*3, sd=(1-ind) + ind*9)
res <- rmx(x, eps.lower = 0.01, eps.upper = 0.1)
bias(res)
}
\keyword{robust}
