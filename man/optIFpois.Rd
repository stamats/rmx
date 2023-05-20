\name{optIF.pois}
\alias{optIF.pois}
\title{Computation of Optimally Robust IFs for Poisson Mean}
\description{
  The function computes the optimally robust IF for Poisson mean. The function 
  is rarely called directly, but via function \code{\link{optIF}} and is mainly 
  for internal use.
}
\usage{
optIF.pois(radius, lambda = 1, aUp = 100*lambda, cUp = 1e4, 
           delta = 1e-9)
}
\arguments{
  \item{radius}{ non-negative real: neighborhood radius.}
  \item{lambda}{ lambda parameter; see \code{\link[stats]{dpois}}. }
  \item{aUp}{ numeric: upper limit for centering constant a. }
  \item{cUp}{ postive real: upper limit for clipping constant c. }
  \item{delta}{ positive real: desired accuracy (convergence tolerance). }
}
\details{The Lagrange multipliers contained in the expression
  of the optimally robust IF are computed; i.e., clipping, centering
  and standardising constant; see Chapter 4 of Kohl (2005).
}
\value{
  An object of class \code{"optIF"} is returned. It contains the 
  following arguments:
  \item{model}{ short name of the model / distribution }
  \item{modelName}{ full name of the model }
  \item{parameter}{ parameter values of the model }
  \item{A}{ standardizing matrix }
  \item{a}{ centering vector }
  \item{b}{ clipping constant }
  \item{IFun}{ influence function }
  \item{range}{ function to generate grid for evaluating \code{IFun} }
  \item{asMSE}{ maximum asymptotic mean squared error (MSE) }
  \item{asVar}{ asymptotic (co)variance }
  \item{asBias}{ maximum asymptotic bias }
  \item{radius}{ neighborhood radius }
}
\references{ 
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{optIF}}}
\examples{
## IF of ML estimator
optIF.pois(radius = 0)
## IF of MB estimator
optIF.pois(radius = Inf)
## IF of optimally robust AL estimator
optIF.pois(radius = 0.5)
}
\keyword{robust}
