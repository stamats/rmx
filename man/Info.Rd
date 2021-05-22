\name{absInfo and relInfo}
\alias{absInfo}
\alias{absInfo.rmx}
\alias{relInfo}
\alias{relInfo.rmx}
\title{Compute Absolute and Relative Information}
\description{
  Functions to compute absolute and relative information of asymptotically 
  linear estimators.
}
\usage{
absInfo(x, \dots)
relInfo(x, \dots)

\method{absInfo}{rmx}(x, \dots)
\method{relInfo}{rmx}(x, \dots)
}
\arguments{
  \item{x}{ object of S3 class \code{rmx}. }
  \item{\dots}{ further arguments passed through. }
}
\details{
  The function is inspired by the respective functions of the RobASt-family of 
  packages.
  
  In case of optimally-robust RMX estimators computed with function \code{rmx}
  (S3 class \code{rmx}), the absolute and relative information is computed. 
  For more details about the concept we refer to Kohl (2005).
}
\value{
  Matrix with the data and the computed absolute and relative information, 
  respectively.
}
\references{
  Kohl, M. (2005). \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.

  M. Kohl, P. Ruckdeschel, and H. Rieder (2010). Infinitesimally Robust Estimation 
  in General Smoothly Parametrized Models. \emph{Statistical Methods and Application}, 
  \bold{19}(3):333-354. 
  
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Rieder, H., Kohl, M. and Ruckdeschel, P. (2008) The Costs of not Knowing
  the Radius. Statistical Methods and Applications \emph{17}(1) 13-40.
  Extended version: \url{http://r-kurs.de/RRlong.pdf}.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{rmx}}, \code{\link{optIF}}}
\examples{
ind <- rbinom(20, size=1, prob=0.05) 
x <- rnorm(20, mean=ind*3, sd=(1-ind) + ind*9)
res <- rmx(x, eps.lower = 0.01, eps.upper = 0.1)
absInfo(res)
relInfo(res)
ML <- rmx(x, eps = 0)
absInfo(ML)
relInfo(ML)
}
\keyword{robust}
