\name{iiPlot}
\alias{iiPlot}
\alias{iiPlot.rmx}
\title{Compare Absolute Information of Estimators}
\description{
  The function \code{iiPlot} compares the absolute information (AI) of 
  RMX estimators with AI of ML estimators.
}
\usage{
iiPlot(x, \dots)

\method{iiPlot}{rmx}(x, param.digits = 2, ggplot.ylim = NULL,
        ggplot.xlab = "Absolute Information of RMX", 
        ggplot.ylab = "Absolute Information of ML",
        ggplot.ggtitle = NULL, color.line = "#E18727",
        point.col = "#0072B5", point.alpha = 0.4, \dots)
}
\arguments{
  \item{x}{ object of S3 class \code{rmx}. }
  \item{param.digits}{ number of digits used for the estimated parameter values, if
                       default \code{ggplot.ggtitle} is used. }
  \item{ggplot.ylim}{ numeric vector of length two defining the limits of the 
                      y-axis. }
  \item{ggplot.xlab}{ label(s) of x-axis, recycled if length is equal to 1 and more
                      than 1 parameter has been estimated. }
  \item{ggplot.ylab}{ label(s) of y-axis, recycled if length is equal to 1 and more
                      than 1 parameter has been estimated. }
  \item{ggplot.ggtitle}{ if \code{NULL} a default plot title is. If it is specified,
                         it should have length 1. }
  \item{color.line}{ single color used for colouring the line that represents the
                     maximum absolute information in case of RMX. }
  \item{point.col}{ single color used for colouring the data points. }
  \item{point.alpha}{ alpha used for color shading.}
  \item{\dots}{ further arguments passed through. }
}
\details{
  The function is inspired by the \code{plot}-methods and function \code{PlotIC} 
  of the RobASt-family of packages.
  
  In case of optimally-robust RMX estimators computed with function \code{rmx}
  (S3 class \code{rmx}), the absolute information is plotted versus the absolute
  information of the ML estimator. 
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
ind <- rbinom(50, size=1, prob=0.05) 
x <- rnorm(50, mean=0, sd=(1-ind) + ind*2)
res <- rmx(x, eps.lower = 0.01, eps.upper = 0.1)
gg <- iiPlot(res)
gg
gg$data

## plot-method
plot(res, which = 7)

## setting and passing argument
iiPlot(res, ggplot.ylim = c(0.5, 10), point.col = "darkred")
plot(res, which = 7, 
     control = list(iiPlot = list(ggplot.ylim = c(0.5, 10), 
                                  point.col = "#0072B5")))
}
\keyword{robust}
