\name{fsRadius}
\alias{fsRadius}
\title{Finite-sample Corrected Radius}
\description{
Given some radius and some sample size the function computes
the corresponding finite-sample corrected radius.
}
\usage{
fsRadius(r, n, model = "norm")
}
\arguments{
  \item{r}{ asymptotic radius (non-negative numeric) }
  \item{n}{ sample size }
  \item{model}{ character: short name of the model/distribution (default = \code{"norm"});
  see also details. }
}
\details{
  The finite-sample correction is based on empirical results obtained via 
  Monte-Carlo simulations. 
  
  Given some radius of a shrinking contamination neighborhood which leads 
  to an asymptotically optimal robust estimator, the finite-sample empirical 
  MSE based on contaminated samples was minimized for the class of 
  optimally robust estimators and the corresponding finite-sample 
  radius determined and saved. 
  
  The computation is based on the saved results of the Monte-Carlo simulations.
  
  As models we have implemented so far:
  \enumerate{
    \item \code{"norm"}: normal location and scale. Parameters can be set 
    via argument \code{mean} and \code{sd}; see examples.
  }
}
\value{Finite-sample corrected radius.}
\references{
  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Rieder, H., Kohl, M. and Ruckdeschel, P. (2008) The Costs of not Knowing
  the Radius. Statistical Methods and Applications \emph{17}(1) 13-40.
  Extended version: \url{http://r-kurs.de/RRlong.pdf}
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{rmx}}, \code{\link{rowRmx}}, \code{\link{colRmx}} }
\examples{
## finite-sample radius is larger (more conservative)
fsRadius(r = 0.5, n = 3, model = "norm")
fsRadius(r = 0.5, n = 10, model = "norm")
fsRadius(r = 0.5, n = 25, model = "norm")
fsRadius(r = 0.5, n = 50, model = "norm")
fsRadius(r = 0.5, n = 100, model = "norm")
fsRadius(r = 0.5, n = 500, model = "norm")
fsRadius(r = 0.5, n = 1000, model = "norm")
}
\keyword{robust}
