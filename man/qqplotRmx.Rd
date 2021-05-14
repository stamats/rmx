\name{qqplotRmx}
\alias{qqplotRmx}
\title{qq-Plot of Model Estimated by Optimally Robust RMX Estimator}
\description{
  The function \code{qqplotRmx} generates a qq-Plot of a model fitted 
  by an optimally robust RMX estimator (S3 class \code{rmx}) computed
  with function \code{rmx}. 
}
\usage{
qqplotRmx(x, param.digits = 3, ggplot.xlab = "Theoretical Quantiles", 
          ggplot.ylab = "Sample Quantiles", ggplot.ggtitle = NULL,
          point.col = "#0072B5", point.alpha = 1)
}
\arguments{
  \item{x}{ object of S3 class \code{rmx}. }
  \item{param.digits}{ number of digits used for the estimated parameter values, if
                       default \code{ggplot.ggtitle} is used. }
  \item{ggplot.xlab}{ label of x-axis. }
  \item{ggplot.ylab}{ label of y-axis. }
  \item{ggplot.ggtitle}{ plot title. If \code{NULL} default plot title is generated. }
  \item{point.col}{ single color used for colouring the data points. }
  \item{point.alpha}{ alpha used for color shading.}
}
\details{
  The function is based on the function provided by packages \pkg{ggplot2} 
  and \pkg{qqplotr}.
}
\value{
  Invisible object of class \code{ggplot}.
}
\references{
  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Rieder, H., Kohl, M. and Ruckdeschel, P. (2008) The Costs of not Knowing
  the Radius. Statistical Methods and Applications \emph{17}(1) 13-40.
  Extended version: \url{http://r-kurs.de/RRlong.pdf}
  
  M. Kohl, P. Ruckdeschel, and H. Rieder (2010). Infinitesimally Robust Estimation 
  in General Smoothly Parametrized Models. \emph{Statistical Methods and Application}, 
  \bold{19}(3):333-354. 
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{rmx}}, \code{\link{optIF}}}
\examples{
ind <- rbinom(100, size=1, prob=0.05) 
x <- rnorm(100, mean=ind*3, sd=(1-ind) + ind*9)
res <- rmx(x, eps.lower = 0.01, eps.upper = 0.1)
qqplotRmx(res)

## plot-method
plot(res, which = 2)

## setting and passing argument
qqplotRmx(res, ggplot.xlab = "Quantiles of estimated model", 
               point.col = "black")
plot(res, which = 2, 
     control = list(qqplotRmx = list(ggplot.xlab = "Quantiles of estimated model", 
                                     point.col = "black")))
}
\keyword{robust}
