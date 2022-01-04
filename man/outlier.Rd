\name{outlier}
\alias{outlier}
\alias{outlier.rmx}
\alias{print.outlier}
\title{Outlier Region}
\description{
  The function \code{outlier} computes outlier regions.
}
\usage{
outlier(x, \dots)

\method{outlier}{rmx}(x, prob = 0.001)

\method{print}{outlier}(x, digits = 3)
}
\arguments{
  \item{x}{ object of S3 class \code{rmx}. }
  \item{prob}{ probability used to define outliers.}
  \item{digits}{ minimal number of significant digits. }
  \item{\dots}{ further arguments passed through. }
}
\details{
  The function is inspired by the outlier rejection rule: median +/- 3 MAD.
  Since \code{pnorm(3)} is about 0.001, we use it as default.
  
  In case of optimally-robust RMX estimators computed with function \code{rmx}
  (S3 class \code{rmx}), the outliers are defined using the respective 
  quantiles of the fitted model. That is, the respective \code{prob} and 
  \code{1-prob} quantiles define the boundaries of the outlier region. In case 
  of normal location and scale, this is equivalent to replacing median and MAD 
  by the respective RMX estimates.
}
\value{
  An object of class \code{"outlier"} is returned. It contails at least the 
  following arguments:
  \item{rmx}{ object of class \code{rmx}. }
  \item{lower}{ lower boundary of outlier region. }
  \item{upper}{ upper boundary of outlier region. }
  \item{prop.outlier}{ proportion of data in the outlier region. }
  \item{p.outlier}{ probability of the outlier region under the fitted model. }
  \item{prop.lower}{ proportion of data in the lower outlier region. }
  \item{prop.upper}{ proportion of data in the upper outlier region. }
  \item{p.lower}{ probability of the lower outlier region under the fitted model. }
  \item{p.upper}{ probability of the upper outlier region under the fitted model. }
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
\seealso{\code{\link{rmx}}, \code{\link{getOutliers}}, \code{\link{cniper}}}
\examples{
ind <- rbinom(100, size=1, prob=0.05) 
x <- rnorm(100, mean=ind*3, sd=(1-ind) + ind*9)
res <- rmx(x, eps.lower = 0.01, eps.upper = 0.1)
outlier(res)
}
\keyword{robust}
