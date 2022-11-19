# <img src="https://github.com/stamats/rmx/raw/main/hex-rmx.png" alt="rmx" width="120"/> &emsp; rmx: Radius-minimax Estimators

The repository includes the development version of R package rmx

[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)


## Description
Radius-minimax (rmx) estimators are a special class of optimally-robust statistical 
procedures (Kohl (2005), <https://epub.uni-bayreuth.de/839/2/DissMKohl.pdf>; 
Rieder et al. (2008), <doi:10.1007/s10260-007-0047-7>). It is a class of 
asymptotically linear estimators where infinitesimal (shrinking) neighborhoods 
around parametric models are assumed (Rieder (1994), ISBN:978-1-4684-0624-5; 
Ruckdeschel (2006), <doi:10.1007/s00184-005-0020-0>; 
Kohl et al. (2010), <doi:10.1007/s10260-010-0133-0>; 
Kohl (2012), <doi:10.1080/02331888.2010.540668>). In addition, various diagnostic 
plots and finite-sample corrections are included.


## Motivation
The rmx estimators are also implemented in the RobASt-family of packages (i.e.,
packages RobAStBase, ROptEst, RobLox, RobExtremes, ROptRegTS, RobRex). However, 
package rmx provides a much simpler interface and lower computation times at 
the cost of being less general and flexible.


## Installation

```{r, eval = FALSE}
## Development version from GitHub
# install.packages("remotes")
remotes::install_github("stamats/rmx")
```


## Getting started

```{r}
library(rmx)
```
