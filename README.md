
# localSCR

<!-- badges: start -->
<!-- badges: end -->

The goal of ‘localSCR’ is to provide user-friendly functions to
implement Bayesian spatial capture-recapture models (Royle et al. 2014)
using the ‘nimble’ package (de Valpine et a. 2022) in R. The package
currently has functions to 1) assist with defining the state-space grid
and extent for a given 2-dimensional or 3-dimensional trap array (i.e.,
when traps are clustered in space), 2) simulate data under different
encounter distributions and other parameters, 3) create habitat masks
from either raster data or spatial polygons, 4) provide template SCR
models that are easily customizable, 5) fit and summarize SCR models
using ‘nimble’ (de Valpine et a. 2022) with options for parallel
processing, and 6) create realized density surfaces from MCMC output.
Future functionality will include discrete state-space models and
implementing localized approaches as in Milleret et al. (2019) and
Woodruff et al. (2020).

## Installation

You can install the development version of ‘localSCR’ like so:

``` r
library(remotes)
install_github("sitkensis22/localSCR")
```

Be sure to see important information about using ‘nimble’ on your
computer (including installing rtools): <https://r-nimble.org/download>.

## Vignettes

-   <a href="https://sitkensis22.github.io/localSCR/articles/classic_scr.html">Classic
    SCR models: continuous state-space and marked individuals</a>
-   <a href="https://sitkensis22.github.io/localSCR/articles/unmarked_scr.html">Spatial
    count models: continuous state-space and unmarked individuals</a>
-   <a href="https://sitkensis22.github.io/localSCR/articles/mark_resight_scr.html">Mark-resight
    models: continuous state-space with both marked and unmarked
    individuals</a>
-   <a href="https://sitkensis22.github.io/localSCR/articles/discrete_scr.html">Discrete
    SCR models: discrete state-space and marked individuals</a>
-   <a href="https://sitkensis22.github.io/localSCR/articles/local_classic_scr.html">Local
    classic SCR models: continuous state-space and marked
    individuals</a>

## References

de Valpine P, C. Paciorek, D. Turek, N. Michaud, C. Anderson-Bergman, F.
Obermeyer, C. C. Wehrhahn, A. Rodrìguez, L. D. Temple, and S. Paganin.
2022. *NIMBLE: MCMC, Particle Filtering, and Programmable Hierarchical
Modeling*. doi: 10.5281/zenodo.1211190 (URL:
<https://doi.org/10.5281/zenodo.1211190>), R package version 0.12.2,
URL:<https://cran.r-project.org/package=nimble>.

Milleret, C., P. Dupont, C. Bonenfant, H. Henrik Brøseth, Ø. Flagstad,
C. Sutherland, and R. Bischof. 2019. A local evaluation of the
individual state‐space to scale up Bayesian spatial capture‐recapture.
Ecology and Evolution 9:352–363.

Royle, J. A., R. B. Chandler, R. Sollmann, and B. Gardner. 2014. Spatial
capture‐recapture. Academic Press, Waltham, Massachusetts, USA.

Woodruff, S., D. R. Eacker, and L. Waits. 2020. Estimating coyote
density in local, discrete Bayesian capture-recapture models. Journal of
Wildlife Management 10.1002/jwmg.21967.
