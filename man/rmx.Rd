\name{rmx}
\alias{rmx}
\alias{rmx.norm}
\alias{rmx.binom}
\alias{print.rmx}
\alias{summary.rmx}
\alias{coef.rmx}
\alias{vcov.rmx}
\alias{plot.rmx}
\title{Optimally Robust RMX Estimator}
\description{
  The function \code{rmx} computes optimally robust rmx estimators. 
  The definition of these estimators can be found in Kohl (2005) and 
  Rieder et al. (2008), respectively. The other functions are rarely 
  called directly and are mainly for internal use.
}
\usage{
rmx(x, model = "norm", eps.lower=0, eps.upper=NULL, eps=NULL, k = 3L, 
    initial.est=NULL, fsCor = NULL, na.rm = TRUE, message = TRUE, \dots)
    
rmx.norm(x, eps.lower=0, eps.upper=NULL, eps=NULL, k = 3L, 
         initial.est=NULL, fsCor = TRUE, na.rm = TRUE, mad0 = 1e-4)

rmx.binom(x, eps.lower=0, eps.upper=NULL, eps=NULL, k = 3L, 
          initial.est=NULL, fsCor = FALSE, na.rm = TRUE, 
          size, M = 10000, parallel = FALSE, ncores = NULL, 
          aUp = 100*size, cUp = 1e4, delta = 1e-9)
          
\method{print}{rmx}(x, digits = getOption("digits"), prefix = " ", \dots)

\method{summary}{rmx}(object, digits = getOption("digits"), prefix = " ", \dots)

\method{coef}{rmx}(object, complete = TRUE, \dots)

\method{vcov}{rmx}(object, \dots)

\method{plot}{rmx}(x, which = 1, 
              control = list(ifPlot = NULL, qqPlot = NULL,
                             ppPlot = NULL, dPlot = NULL,
                             aiPlot = NULL, riPlot = NULL,
                             iiPlot = NULL, plot = TRUE), \dots)
}
\arguments{
  \item{x}{ numeric vector \code{x} of data values. Object of class \code{rmx}
  in case of \code{print} and \code{plot} method. }
  \item{object}{ object of class \code{rmx}. }
  \item{model}{ character: short name of the model/distribution 
                (default = \code{"norm"}); see also details. }
  \item{eps.lower}{ positive real (0 <= \code{eps.lower} <= \code{eps.upper}): 
        lower bound for the amount of gross errors; see details below. }
  \item{eps.upper}{ positive real (\code{eps.lower} <= \code{eps.upper} <= 0.5): 
        upper bound for the amount of gross errors; see details below. }
  \item{eps}{ positive real (0 < \code{eps} <= 0.5): amount of gross errors. 
        See details below. }
  \item{k}{ positive integer: k-step is used to compute the optimally robust estimator. }
  \item{initial.est}{ initial estimate for \code{mean} and \code{sd}. If missing 
        median and MAD are used. }
  \item{fsCor}{ \code{NULL} or logical: perform finite-sample correction; see 
  function \code{\link{fsRadius}}. }
  \item{na.rm}{ logical: if \code{TRUE}, \code{NA} values are removed before the estimator is evaluated.}
  \item{message}{ logical: if \code{FALSE}, messages are suppressed.}
  \item{mad0}{ if \code{mad(x)} equal to 0, it is replaced by \code{mad0}. }
  \item{size}{ size parameter (known!); see \code{\link[stats]{dbinom}}. }
  \item{M}{ number of Monte-Carlo simulations; see \code{\link{fsRadius}}.}
  \item{parallel}{ logical: use package \pkg{parallel} for computations. }
  \item{ncores}{ if \code{parallel = TRUE}: number of cores used for computations.
                 If missing, the maximum number of cores - 1 is used. }
  \item{aUp}{ numeric: upper limit for centering constant a. }
  \item{cUp}{ postive real: upper limit for clipping constant c. }
  \item{delta}{ positive real: desired accuracy (convergence tolerance). }
  \item{digits}{ minimal number of significant digits. }
  \item{prefix}{ string, passed to \code{\link[base]{strwrap}}. }
  \item{complete}{ logical indicating if the full coefficient vector should
  be returned; see \code{\link[stats]{coef}}. }
  \item{which}{ numeric number from 1 to 7 indicating the plot that shall be
  generated; see details below. }
  \item{control}{ list of parameters passed through to the respective plot. }
  \item{\dots}{ further arguments passed through; e.g., known parameters such as 
                \code{size} in case of the binomial model. }
}
\details{
  If the amount of gross errors (contamination) is known, it can be 
  specified by \code{eps}. The radius of the corresponding infinitesimal 
  contamination neighborhood is obtained by multiplying \code{eps} 
  by the square root of the sample size. 

  If the amount of gross errors (contamination) is unknown, try to find a 
  rough estimate for the amount of gross errors, such that it lies 
  between \code{eps.lower} and \code{eps.upper}.
  
  If neither \code{eps} nor \code{eps.upper} is provided, \code{eps.upper}
  will be estimated by applying function \code{outlier} to the RMX estimator
  with \code{eps.lower = 0} and \code{eps.upper = 0.5}.
  
  As models we have implemented so far:
  \enumerate{
    \item \code{"norm"}: normal location (mean) and scale (sd).
    \item \code{"binom"}: binomial probability (size known!).
  }
  
  As plots we have implemented:
  \enumerate{
    \item \code{\link{ifPlot}}: plot of IF.
    \item \code{\link{qqPlot}}: qq-plot for fitted model, only implemented
    for continuous models.
    \item \code{\link{ppPlot}}: pp-plot for fitted model, only implemented
    for continuous models.
    \item \code{\link{dPlot}}: density plot for fitted model.
    \item \code{\link{aiPlot}}: plot of absolute information.
    \item \code{\link{riPlot}}: plot of relative information, only implemented
    for models where at least two parameters have to be estimated.
    \item \code{\link{iiPlot}}: compare absolute information of RMX estimator
    with absolute information of ML estimator.
  }
}
\value{
  An object of class \code{"rmx"} is returned. It contails at least the 
  following arguments:
  \item{rmxEst}{ estimates }
  \item{rmxIF}{ object of class \code{optIF}; see \code{\link{optIF}}. }
  \item{initial.est}{ initial estimates. }
  \item{Infos}{ matrix with information about the estimator }
  \item{x}{ data used for the estimation. }
  \item{n}{ sample size }
  \item{eps.lower}{ lower bound for the amount of gross errors, if provided 
                    otherwise \code{NA}. }
  \item{eps.upper}{ upper bound for the amount of gross errors, if provided 
                    otherwise \code{NA}. }
  \item{eps}{ amount of gross errors, if provided otherwise \code{NA}. }
  \item{fsCor}{ finite-sample correction }
  \item{k}{ k-step construction }
  \item{call}{ matched call }
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
\seealso{\code{\link{optIF}}, \code{\link{rowRmx}}, \code{\link{fsRadius}}}
\examples{
###########################################################
## For more details see vignette (A Guide to rmx)
###########################################################

## normal location (mean) and scale (sd)
ind <- rbinom(100, size=1, prob=0.05) 
x <- rnorm(100, mean=ind*3, sd=(1-ind) + ind*9)
res1 <- rmx(x, eps.lower = 0.01, eps.upper = 0.1)
res1
summary(res1)
confint(res1) # method = "as"
confint(res1, method = "as.bias")
plot(res1, which = 1)
plot(res1, which = 2)

## binomial probability
ind <- rbinom(100, size=1, prob=0.05) 
x <- rbinom(100, size=5, prob=(1-ind)*0.1 + ind*0.5)
res2 <- rmx(x, model = "binom", size = 5, eps.lower = 0.01, eps.upper = 0.1)
res2
summary(res2)
confint(res2) # method = "as"
confint(res2, method = "as.bias")
plot(res2, which = 1)
plot(res2, which = 4)
}
\keyword{robust}
