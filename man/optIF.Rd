\name{optIF}
\alias{optIF}
\title{Computation of the optimally robust IF for AL estimators}
\description{
  The function computes the optimally robust IF for AL estimators 
  in case of shrinking contamination neighborhoods / gross error models. 
  The definition of these estimators can be found in Rieder (1994) and 
  Kohl (2005).
}
\usage{
optIF(model = "norm", radius = NULL, check = FALSE, ...)
}
\arguments{
  \item{model}{ character: short name of the model/distribution (default = \code{"norm"});
  see also details. }
  \item{radius}{ non-negative real: neighborhood radius.}
  \item{check}{ logical: check accuracy of the computation.}
  \item{\dots}{ further model specific parameters such as parameter value(s); 
  see details. }
}
\details{The Lagrange multipliers contained in the expression
  of the optimally robust IF are computed; i.e., clipping constant, centering
  vector and standardising matrix. For more details we refer to Rieder (1994) 
  or Kohl (2005).
  
  As models we have implemented so far:
  \enumerate{
    \item \code{"norm"}: normal location (mean) and scale (sd); see also
    \code{\link{optIF.norm}}.
    \item \code{"binom"}: binomial probability (size known); see also
    \code{\link{optIF.binom}}.
  }
}
\value{
  An object of class \code{"optIF"} is returned. It contains at least the 
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
  \item{call}{ matched call }
}
\references{ 
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{rmx}}}
\examples{
## normal location and scale
IF <- optIF(model = "norm", radius = 0.5, mean = -2, sd = 3, check = TRUE)
IF
summary(IF)
plot(IF)
## binomial probability
IF <- optIF(model = "binom", radius = 0.5, size = 10, check = TRUE)
IF
summary(IF)
plot(IF)
}
\keyword{robust}
