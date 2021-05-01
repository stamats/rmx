# rmx: Radius-minimax Estimators

The repository includes the development version of R package rmx

[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

## Description
Radius-minimax estimators are a special class of optimally-robust statistical 
procedures (Kohl (2005), <https://epub.uni-bayreuth.de/839/2/DissMKohl.pdf>; 
Rieder et al. (2008), <doi:10.1007/s10260-007-0047-7>). It is a class of 
asymptotically linear estimators where infinitesimal (shrinking) neighborhoods 
around parametric models are assumed (Rieder (1994), ISBN:978-1-4684-0624-5; 
Ruckdeschel (2006), <doi:10.1007/s00184-005-0020-0>; 
Kohl et al. (2010), <doi:10.1007/s10260-010-0133-0>; 
Kohl (2012), <doi:10.1080/02331888.2010.540668>).

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
