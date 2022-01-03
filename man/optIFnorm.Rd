\name{optIF.norm}
\alias{optIF.norm}
\title{Computation of Optimally Robust IFs for Normal Location and Scale}
\description{
  The function computes the optimally robust IF for normal location (mean)
  and scale (sd). The function is rarely called directly, but via function
  \code{\link{optIF}} and is mainly for internal use.
}
\usage{
optIF.norm(radius, mean = 0, sd = 1, A.loc.start = 1, A.sc.start = 0.5, 
           a.sc.start = 0, bUp = 1000, delta = 1e-6, itmax = 100L, check = FALSE)
}
\arguments{
  \item{radius}{ non-negative real: neighborhood radius.}
  \item{mean}{ mean parameter; see \code{\link[stats]{dnorm}}. }
  \item{sd}{ sd parameter; see \code{\link[stats]{dnorm}}. }
  \item{A.loc.start}{ numeric: initial value for location component of 
  standardising matrix A. }
  \item{A.sc.start}{ numeric: initial value for scale component of 
  standardising matrix A. }
  \item{a.sc.start}{ numeric: initial value for scale component of 
  centering vector a. }
  \item{bUp}{ postive real: upper limit for clipping constant b. }
  \item{delta}{ positive real: desired accuracy (convergence tolerance). }
  \item{itmax}{ positive ingeger: maximum number of iterations. }
  \item{check}{ logical: check accuracy of the computation (only for
  0 < radius < Inf.}
}
\details{The Lagrange multipliers contained in the expression
  of the optimally robust IF are computed; i.e., clipping constant, centering
  vector and standardising matrix; see Chapter 8 of Kohl (2005).
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
optIF.norm(radius = 0)
## IF of MB estimator
optIF.norm(radius = Inf)
## IF of optimally robust AL estimator
optIF.norm(radius = 0.5, check = TRUE)
}
\keyword{robust}
