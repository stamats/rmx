\name{confint}
\alias{confint.rmx}
\title{Confidence Intervals for RMX Estimates}
\description{
  The function \code{confint} computes confindence interval for rmx estimates.
}
\usage{
\method{confint}{rmx}(object, parm, level = 0.95, method = "as", R = 9999, type = "all", ...)
}
\arguments{
  \item{object}{ object of S3 class \code{rmx}. }
  \item{parm}{ see \code{\link[stats]{confint}}. A specification of which parameters 
               are to be given confidence intervals, either a vector of numbers 
               or a vector of names. If missing, all parameters are considered. }
  \item{level}{ see \code{\link[stats]{confint}}, the confidence level required. }
  \item{method}{ see details.}
  \item{R}{ number of bootstrap replicates.}
  \item{type}{ A vector of character strings representing the type of intervals required. 
  The value should be any subset of the values c("norm","basic", "stud", "perc", "bca") 
  or simply "all" which will compute all five types of intervals; see also 
  \code{\link[boot]{boot.ci}}.}
  \item{\dots}{ further arguments passed through especially to \code{rowRmx}. }
}
\details{
  The function is inspired by the respective function of the RobASt-family of 
  packages.
  
  In case of optimally-robust RMX estimators computed with function \code{rmx}
  (S3 class \code{rmx}), confint intervals are computed either only using 
  the asymptotic (co)variance (\code{method = "as"}, default) or using the
  asymptotic (co)variance and the maximum asymptotic bias (\code{method = "as.bias"})
  or using bootstrap (\code{method = "boot"}).
}
\value{
  An object of class \code{"rmxCI"} is returned. It contails at least the 
  following arguments:
  \item{method}{ method selected for computing the confidence interval. }
  \item{conf.int}{ numeric matrix with the confidence bounds. }
  \item{rmxEst}{ point estimates. }
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
\seealso{\code{\link{rmx}}, \code{\link[stats]{confint}}}
\examples{
ind <- rbinom(100, size=1, prob=0.05) 
x <- rnorm(100, mean=ind*3, sd=(1-ind) + ind*9)
res <- rmx(x, eps.lower = 0.01, eps.upper = 0.1)
confint(res)
confint(res, method = "as.bias")
#confint(res, method = "boot", R = 999)
}
\keyword{robust}
