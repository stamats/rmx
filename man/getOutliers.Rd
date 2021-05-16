\name{getOutliers}
\alias{getOutliers}
\title{Outlier Points}
\description{
  The function \code{getOutliers} finds the outlier points in the data.
}
\usage{
getOutliers(x, \dots)
}
\arguments{
  \item{x}{ object of S3 class \code{outlier} or \code{rmx}. }
  \item{\dots}{ further arguments passed through. }
}
\details{
  The function determines the outliers points in the data used to compute 
  the RMX estimator.
}
\value{
  A list with at least the following arguments:
  \item{values}{ object of class \code{rmx}. }
  \item{indices}{ lower boundary of cniper region. }
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
  
  Ruckdeschel, P. (2010). Consequences of Higher Order Asymptotics for the MSE 
  of M-estimators on Neighborhoods. arXiv, \url{https://arxiv.org/abs/1006.0123}.

}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{rmx}}, \code{\link{cniper}}, \code{\link{getOutliers}}}
\examples{
ind <- rbinom(100, size=1, prob=0.05) 
x <- rnorm(100, mean=ind*3, sd=(1-ind) + ind*9)
res <- rmx(x, eps.lower = 0.01, eps.upper = 0.1)
(out.res <- outlier(res))
getOutliers(out.res)
}
\keyword{robust}
