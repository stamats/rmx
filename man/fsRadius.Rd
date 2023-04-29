\name{fsRadius}
\alias{fsRadius}
\alias{fsRadius.norm}
\alias{fsRadius.binom}
\alias{fsRadius.pois}
\title{Finite-sample Corrected Radius}
\description{
Given some radius and some sample size the function computes
the corresponding finite-sample corrected radius. The user should use
function \code{fsRadius}. The other functions are rarely called directly and 
are mainly for internal use.
}
\usage{
fsRadius(r, n, model = "norm", \dots)

fsRadius.norm(r, n)

fsRadius.binom(r, n, prob, size, M = 10000, parallel = FALSE, ncores = NULL)

fsRadius.pois(r, n, lambda, M = 10000, parallel = FALSE, ncores = NULL)
}
\arguments{
  \item{r}{ asymptotic radius (non-negative numeric) }
  \item{n}{ sample size }
  \item{model}{ character: short name of the model/distribution (default = \code{"norm"});
  see also details. }
  \item{prob}{ prob parameter; see \code{\link[stats]{dbinom}}. }
  \item{size}{ size parameter (known!); see \code{\link[stats]{dbinom}}. }
  \item{lambda}{ lambda parameter; see \code{\link[stats]{dpois}}. }
  \item{M}{ number of Monte-Carlo simulations; see details below.}
  \item{parallel}{ logical: use package \pkg{parallel} for computations. }
  \item{ncores}{ if \code{parallel = TRUE}: number of cores used for computations.
                 If missing, the maximum number of cores - 1 is used. }
  \item{\dots}{ further arguments passed through; e.g., known parameters such as 
                \code{size} in case of the binomial model. }
}
\details{
  The finite-sample correction is based on empirical results obtained via 
  Monte-Carlo simulations. 
  
  Given some radius of a shrinking contamination neighborhood which leads 
  to an asymptotically optimal robust estimator, the finite-sample empirical 
  MSE based on contaminated samples is minimized for the class of 
  optimally robust estimators and the corresponding finite-sample 
  radius determined. 
  
  For some models (\code{"norm"}) the computation is based on the saved results 
  of the Monte-Carlo simulations, whereas for other models (\code{"binom"}, 
  \code{"pois"}) the results are computed inside of the function based on 
  \code{M} Monte-Carlo simulations.
  
  As models we have implemented so far:
  \enumerate{
    \item \code{"norm"}: normal location (mean) and scale (sd)
    \item \code{"binom"}: probability of success (prob)
    \item \code{"pois"}: Poisson mean (lambda)
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
\seealso{\code{\link{rmx}}, \code{\link{rowRmx}} }
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
