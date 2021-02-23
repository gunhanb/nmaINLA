nmaINLA R package
===============
[![CRAN status](https://www.r-pkg.org/badges/version/nmaINLA)](https://cran.r-project.org/package=nmaINLA)

[![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/nmaINLA)](https://cranlogs.r-pkg.org/badges/grand-total/nmaINLA)


`nmaINLA` is an R package for Network Meta-Analysis using INLA. It is the accompanying R package for our paper: Günhan BK, Friede T, Held L. A design‐by‐treatment interaction model for network meta‐analysis and meta‐regression with integrated nested Laplace approximations. Reserach Synthesis Methods. 2018;9:179–194. https://doi.org/10.1002/jrsm.1285.


Building
--------
The installation of R package 'INLA' is compulsory for successful usage. The 'INLA' package is not available on CRAN, and it can be obtained from <http://www.r-inla.org>. We recommend the testing version, which can be downloaded by running:

```{r, echo=TRUE, eval=FALSE}
install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
```


Stable version of `nmaINLA` package is available on CRAN. It can be installed as follows:

```{r, echo=TRUE, eval=FALSE}
install.packages("nmaINLA")
```

For how to use the package, type

```{r, echo=TRUE, eval=FALSE}
vignette("nmaINLA")
```


Testing
-------

The `testthat` package is used for testing. Tests reside in the
`tests/testthat` directory. 


License
-------

    nmaINLA is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    nmaINLA is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with nmaINLA  If not, see <http://www.gnu.org/licenses/>.

