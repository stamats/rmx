\name{optIF.exp}
\alias{optIF.exp}
\title{Computation of Optimally Robust IFs for Exponential Scale}
\description{
  The function computes the optimally robust IF for Exponential scale. The function 
  is rarely called directly, but via function \code{\link{optIF}} and is mainly 
  for internal use.
}
\usage{
optIF.exp(radius, scale = 1, aUp = 1, cUp = 0.97, delta = 1e-9)
}
\arguments{
  \item{radius}{ non-negative real: neighborhood radius.}
  \item{scale}{ scale parameter; \code{rate} of \code{\link[stats]{dexp}}
  corresponds to \code{1/scale}. }
  \item{aUp}{ numeric: upper limit for centering constant a. }
  \item{cUp}{ postive real: upper limit for clipping constant c. }
  \item{delta}{ positive real: desired accuracy (convergence tolerance). }
}
\details{The Lagrange multipliers contained in the expression
  of the optimally robust IF are computed; i.e., clipping, centering
  and standardising constant; see Chapter 5 of Kohl (2005).
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
optIF.exp(radius = 0)
## IF of MB estimator
optIF.exp(radius = Inf)
## IF of optimally robust AL estimator
optIF.exp(radius = 0.5)
}
\keyword{robust}
