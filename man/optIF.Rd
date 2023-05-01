\name{optIF}
\alias{optIF}
\alias{print.optIF}
\alias{summary.optIF}
\alias{plot.optIF}
\title{Computation of the optimally robust IF for AL estimators}
\description{
  The function computes the optimally robust IF for AL estimators 
  in case of shrinking contamination neighborhoods / gross error models. 
  The definition of these estimators can be found in Rieder (1994) and 
  Kohl (2005).
}
\usage{
optIF(model = "norm", radius = NULL, \dots)

\method{print}{optIF}(x, digits = getOption("digits"), prefix = " ", \dots)

\method{summary}{optIF}(object, digits = getOption("digits"), prefix = " ", \dots)

\method{plot}{optIF}(x, alpha = 1e-6, digits = 2, plot = TRUE, n = 501, \dots)
}
\arguments{
  \item{model}{ character: short name of the model/distribution (default = \code{"norm"});
  see also details. }
  \item{radius}{ non-negative real: neighborhood radius. }
  \item{x}{ an \code{optIF} object. }
  \item{object}{ an \code{optIF} object. }
  \item{alpha}{ alpha-quantile used to compute plot range, which is from
  \code{alpha}-quantile to 1-\code{alpha}-quantile of the fitted model.}
  \item{digits}{ minimal number of significant digits. }
  \item{prefix}{ string, passed to \code{\link[base]{strwrap}}. }
  \item{plot}{ logical: generate plot, sometimes usefull to suppress plotting
  and modify the returned \code{ggplot} object. }
  \item{n}{ number of points to evaluate the IF in case of continuous models. }
  \item{\dots}{ further arguments passed through sucha as model parameter(s); 
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
    \item \code{"pois"}: Poisson mean; see also \code{\link{optIF.pois}}.
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
IF1 <- optIF(model = "norm", radius = 0.5, mean = -2, sd = 3)
IF1
summary(IF1)
plot(IF1)

## binomial probability
IF2 <- optIF(model = "binom", radius = 0.5, size = 10)
IF2
summary(IF2)
plot(IF2)

## Poisson mean
IF3 <- optIF(model = "pois", radius = 0.5)
IF3
summary(IF3)
plot(IF3)
}
\keyword{robust}
