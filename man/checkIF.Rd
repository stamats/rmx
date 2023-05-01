\name{checkIF}
\alias{checkIF}
\alias{checkIF.optIF}
\alias{checkIF.rmx}
\alias{print.checkIF}
\title{Check IF Properties}
\description{
  The function \code{checkIF} checkes whether the IF properties are fulfilled, 
  which are Fisher consistency, centring and for positive finite radius the 
  MSE equation. For more details see Rieder (1994) or Kohl (2005).
}
\usage{
checkIF(x, \dots)

\method{checkIF}{optIF}(x, rel.tol = .Machine$double.eps^0.5, \dots)

\method{checkIF}{rmx}(x, rel.tol = .Machine$double.eps^0.5, \dots)

\method{print}{checkIF}(x, digits = getOption("digits"), prefix = " ", \dots)
}
\arguments{
  \item{x}{ object of S3 class \code{optIF} or \code{rmx}. }
  \item{rel.tol}{ relative accuracy requested; see \code{\link[stats]{integrate}}. }
  \item{digits}{ minimal number of significant digits. }
  \item{prefix}{ string, passed to \code{\link[base]{strwrap}}. }
  \item{\dots}{ further arguments passed through. }
}
\details{
  The function is inspired by the function \code{checkIC} of the RobASt-family 
  of packages.
}
\value{
  An object of class \code{"checkIF"} is returned. It contails at least the 
  following arguments:
  \item{Fisher}{ precision of Fisher consistency condition. }
  \item{center}{ precision of centering. }
  \item{mse.eq}{ precision of MSE equation. \code{NULL}, if neighborhood 
  \code{radius = 0} (ML estimator) or \code{radius = Inf} (minimum bias solution). }
}
\references{
  Kohl, M. (2005). \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
  
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{optIF}}, \code{\link{rmx}}}
\examples{
IF <- optIF(radius = 0.1)
checkIF(IF)
RMX <- rmx(rnorm(20))
checkIF(RMX)
}
\keyword{robust}
