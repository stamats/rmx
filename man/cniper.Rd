\name{cniper}
\alias{cniper}
\alias{cniper.rmx}
\title{Cniper Region}
\description{
  The function \code{cniper} computes cniper regions by comparing the 
  risk of two estimators under contaminations by Dirac point measures.
}
\usage{
cniper(x, \dots)

\method{cniper}{rmx}(x, range.alpha = 1e-6, \dots)
}
\arguments{
  \item{x}{ object of S3 class \code{rmx}. }
  \item{range.alpha}{ alpha-quantile used to determine the search region for
                      cniper points. }
  \item{\dots}{ further arguments passed through. }
}
\details{
  The function is inspired by the respective functions of the RobASt-family of 
  packages.
  
  In case of optimally-robust RMX estimators computed with function \code{rmx}
  (S3 class \code{rmx}), cniper regions are computed by comparing the 
  maximum asymptotic MSE of the RMX estimator with the maximum asyptotic MSE
  of the ML estimator. For more details about the cniper concept we refer 
  to the Introduction of Kohl (2005) and Section 5 of Ruckdeschel (2010).
}
\value{
  An object of class \code{"cniper"} is returned. It contails at least the 
  following arguments:
  \item{rmx}{ object of class \code{rmx}. }
  \item{lower}{ lower boundary of cniper region. }
  \item{upper}{ upper boundary of cniper region. }
  \item{prop.cniper}{ proportion of data in the cniper region. }
  \item{p.cniper}{ probability of the cniper region under the fitted model. }
  \item{prop.lower}{ proportion of data in the lower cniper region. }
  \item{prop.upper}{ proportion of data in the upper cniper region. }
  \item{p.lower}{ probability of the lower cniper region under the fitted model. }
  \item{p.upper}{ probability of the upper cniper region under the fitted model. }
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
\seealso{\code{\link{rmx}}, \code{\link{getCnipers}}, \code{\link{outlier}}}
\examples{
ind <- rbinom(100, size=1, prob=0.05) 
x <- rnorm(100, mean=ind*3, sd=(1-ind) + ind*9)
res <- rmx(x, eps.lower = 0.01, eps.upper = 0.1)
(cni.res <- cniper(res))
plot(cni.res)
}
\keyword{robust}
