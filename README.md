EpiLPS: a fast and flexible Bayesian tool for estimation of the
time-varying reproduction number
================
Oswaldo Gressani

<!-- Introduce badges -->

![Version](https://img.shields.io/badge/Version-1.0.8-lightgrey)
![Languages](https://img.shields.io/badge/Languages-R%2C%20C%2B%2B-informational)
![Lifecycle](https://img.shields.io/badge/lifecycle-postexperimental-yellow)
![CodeSize](https://img.shields.io/github/languages/code-size/oswaldogressani/EpiLPS?color=orange&label=Code%20size&style=plastic)
![UserStars](https://img.shields.io/github/stars/oswaldogressani/EpiLPS?style=social)
![MyTwitter](https://img.shields.io/twitter/follow/OswaldoGressani?style=social)

<img src="man/figures/gplv3-or-later.png" width="10%" style="display: block; margin: auto auto auto 0;" />

## Aim and scope

EpiLPS <span style="color: blue;"> (Gressani et al. 2022)</span>, is the
acronym for **Epi**demiological modeling (tool) with
**L**aplacian-**P**-**S**plines. It proposes a new Bayesian methodology
for estimating the instantaneous reproduction number $R_t$, i.e. the
average number of secondary cases generated by an infectious agent at
time $t$ <span style="color: blue;"> (White et al., 2020) </span>; a key
metric for assessing the transmission dynamics of an infectious disease
and a useful reference for guiding interventions of governmental
institutions in a public health crisis. The EpiLPS project builds upon
two strong pillars in the statistical literature, namely Bayesian
P-splines and Laplace approximations, to deliver a fast and flexible
methodology for inference on $R_t$. EpiLPS requires two (external)
inputs: (1) a time series of incidence counts and (2) a serial or
generation interval distribution.

The underlying model for smoothing incidence counts is based on the
Negative Binomial distribution to account for possible overdispersion in
the data. EpiLPS has a *two-phase engine* for estimating $R_t$. First,
Bayesian P-splines are used to smooth the epidemic curve and to compute
estimates of the mean incidence count of the susceptible population at
each day $t$ of the outbreak. Second, in the philosophy of <span
style="color: blue;"> Fraser (2007)</span>, the renewal equation is used
as a bridge to link the estimated mean incidence and the reproduction
number. As such, the second phase delivers a closed-form expression of
$R_t$ as a function of the B-spline coefficients and the serial interval
distribution.

Another key strength of EpiLPS is that two different strategies can be
used to estimate $R_t$. The first approach called LPSMAP is completely
sampling-free and fully relies on Laplace approximations and *maximum a
posteriori* computation of model hyperparameters for estimation.
Routines for Laplace approximations and B-splines evaluations are
typically the ones that are computationally most intensive and are
therefore coded in C++ and integrated in R via the [Rcpp
package](https://www.rcpp.org/), making them time irrelevant. The second
approach is called LPSMALA (Laplacian-P-splines with a
Metropolis-adjusted Langevin algorithm) and is fully stochastic. It
samples the posterior of the model by using Langevin diffusions in a
Metropolis-within-Gibbs framework. Of course, LPSMAP has a computational
advantage over LPSMALA. Thanks to the lightning fast implementation of
Laplace approximations, LPSMAP typically delivers estimate of $R_t$ in a
fraction of a second. The chart below, summarizes the skeleton and
mechanisms behind EpiLPS for estimating $R_t$.

<br> <!-- Include a white space -->

<img src="man/figures/EpiLPSchart.jpg" alt="The EpiLPS tool: path from incidence data and serial interval to the estimated reproduction number." width="80%" style="display: block; margin: auto;" />

## Getting started

As the EpiLPS package includes C++ code, Windows users will need to
install Rtools to include the required compilers for a smooth
experience. Rtools is free and can be downloaded from
<https://cran.r-project.org/bin/windows/Rtools/>. To install the Github
version of EpiLPS (with
[devtools](https://cran.r-project.org/package=devtools)) type the
following lines in the R console:

``` r
install.packages("devtools")
devtools::install_github("oswaldogressani/EpiLPS")
```

The package can then be loaded as follows:

``` r
library("EpiLPS")
```

The EpiLPS package structure is fairly simple. The most important
routines are:

- `estimR()` The main routine to estimate the reproduction number.
- `plot.epilps()` S3 method to plot an object of class `epilps`.
- `episim()` A routine to simulate epidemic data.

## A simulated example

To simulate data with `episim()`, a serial interval distribution and a
pattern for the true reproduction number curve has to be specified. Six
patterns are available for the moment. The data generating process is
based on Poisson counts or negative binomial counts and the epidemic
renewal equation for establishing the link between the mean number of
infections and the reproduction number. The default duration of the
simulated outbreak is 50 days but other choices are possible. The code
below simulates an epidemic according to pattern 4 and gives summarizing
plots by setting the option `plotsim = TRUE`:

``` r
si <- c(0.12, 0.28, 0.30, 0.25, 0.05) # Specify a serial interval distribution
simepi <- episim(serial_interval = si, Rpattern = 4, plotsim = TRUE)
```

<img src="README_files/figure-gfm/simex1-1.png" style="display: block; margin: auto;" />

<br>

The simulated incidence count data can be accessed by typing:

``` r
simepi$y
```

    ##  [1]  10   4  13  14  18  34  39  50  75 114 122 164 192 234 297 324 341 413 451
    ## [20] 476 465 484 472 441 431 386 338 330 248 198 197 151 106 106  74  60  44  35
    ## [39]  20  17  10  10  10   5   2   2   0   1   0   1

The `estimR()` routine can be used to fit the epidemic data. By default,
the LPSMAP approach is used with 30 B-splines in the interval $[1;50]$
and a second order penalty. The `plot()` routine on the `epifit_LPSMAP`
object can be used to plot the estimated reproduction number.

``` r
epifit_LPSMAP <- estimR(incidence = simepi$y, si = si)
plot(epifit_LPSMAP)
```

<img src="README_files/figure-gfm/lpsmap-1.png" style="display: block; margin: auto;" />

<br>

Several options can be specified in the `plot()` routine. For instance,
graphical parameters such as `theme` and `col` can be used to control
the theme and color of the fitted $R_t$ curve. .

``` r
plot(epifit_LPSMAP, theme = "light", col = "steelblue")
```

<img src="README_files/figure-gfm/plotepi-1.png" style="display: block; margin: auto;" />

<br>

Estimates related to $R(t)$ at days $t=8,\dots,14$ can be obtained by
typing:

``` r
knitr::kable(epifit_LPSMAP$RLPS[8:14,])
```

|     | Time |        R |      Rsd |  Rq0.025 |   Rq0.05 |   Rq0.25 |   Rq0.50 |   Rq0.75 |   Rq0.95 |  Rq0.975 |
|:----|-----:|---------:|---------:|---------:|---------:|---------:|---------:|---------:|---------:|---------:|
| 8   |    8 | 2.347971 | 1.429948 | 2.071469 | 2.113619 | 2.248883 | 2.347971 | 2.451425 | 2.608308 | 2.661382 |
| 9   |    9 | 2.316242 | 1.408964 | 2.083892 | 2.119611 | 2.233496 | 2.316242 | 2.402053 | 2.531113 | 2.574498 |
| 10  |   10 | 2.239280 | 1.361276 | 2.039719 | 2.070560 | 2.168492 | 2.239280 | 2.312378 | 2.421748 | 2.458366 |
| 11  |   11 | 2.118863 | 1.287335 | 1.954625 | 1.980144 | 2.060842 | 2.118863 | 2.178519 | 2.267301 | 2.296902 |
| 12  |   12 | 1.983053 | 1.204475 | 1.842537 | 1.864437 | 1.933527 | 1.983053 | 2.033848 | 2.109216 | 2.134286 |
| 13  |   13 | 1.856523 | 1.127282 | 1.739468 | 1.757777 | 1.815377 | 1.856523 | 1.898602 | 1.960817 | 1.981456 |
| 14  |   14 | 1.746807 | 1.060498 | 1.644445 | 1.660488 | 1.710881 | 1.746807 | 1.783488 | 1.837614 | 1.855541 |

## Real data examples

To illustrate EpiLPS on real data, we work with the Covid19 R Interface
Data Hub <https://covid19datahub.io/>. Four countries are considered
(Luxembourg, Italy, Canada and Japan) and the reproduction number is
estimated with LPSMAP over the period April 2020 - October 2021 with a
uniform serial interval over 5 days.

``` r
library("COVID19")

# Uniform serial interval over 5 days
si <- c(0.2, 0.2, 0.2, 0.2, 0.2)

# Luxembourg
Luxembourg <- covid19(country = "LUX", level = 1, verbose = FALSE)
dateLUX <- Luxembourg$date[75:649]
inciLUX <- Luxembourg$hosp[75:649]

# Italy
Italy <- covid19(country = "ITA", level = 1, verbose = FALSE)
dateITA <- Italy$date[42:616]
inciITA <- Italy$hosp[42:616]

# Canada
Canada <- covid19(country = "CAN", level = 1, verbose = FALSE)
dateCAN <- Canada$date[92:645]
inciCAN <- Canada$hosp[92:645]

# Japan
Japan<- covid19(country = "JPN", level = 1, verbose = FALSE)
dateJPN <- Japan$date[95:649]
inciJPN <- Japan$hosp[95:649]

# Fit with EpiLPS
epiLUX <- estimR(incidence = inciLUX, si = si)
epiITA <- estimR(incidence = inciITA, si = si)
epiCAN <- estimR(incidence = inciCAN, si = si)
epiJPN <- estimR(incidence = inciJPN, si = si)

gridExtra::grid.arrange(
plot(epiLUX, dates = dateLUX, datelab = "3m", tcol = "steelblue",
     title = "Estimated R Luxembourg"),
plot(epiITA, dates = dateITA, datelab = "3m", rtcol = "chartreuse4",
     title = "Estimated R Italy"),
plot(epiCAN, dates = dateCAN, datelab = "3m", rtcol = "brown2",
     title = "Estimated R Canada"),
plot(epiJPN, dates = dateJPN, datelab = "3m", rtcol = "darkorchid1",
     title = "Estimated R Japan"),
nrow = 2, ncol = 2)
```

<img src="README_files/figure-gfm/realdata-1.png" style="display: block; margin: auto;" />

## Package version

This is version 1.0.8 (2023-02-12) - “EpiLPS Kernels”.

## Acknowledgments

This project is funded by the European Union’s Research and Innovation
Action under the H2020 work programme, EpiPose (grant number 101003688).

## License

EpiLPS: a fast and flexible Bayesian tool for estimation of the
time-varying reproduction number. Copyright (C) 2021-2023 Oswaldo
Gressani.

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along
with this program. If not, see <https://www.gnu.org/licenses/>.

## References

Gressani, O., Wallinga, J., Althaus, C. L., Hens, N. and Faes, C.
(2022). EpiLPS: A fast and flexible Bayesian tool for estimation of the
time-varying reproduction number. *PLoS Comput Biol* **18**(10):
e1010618. <https://doi.org/10.1371/journal.pcbi.1010618>

White, L. F., Moser, C. B., Thompson, R. N., & Pagano, M. (2021).
Statistical estimation of the reproductive number from case notification
data. *American Journal of Epidemiology*, **190**(4), 611-620.

Fraser C (2007) Estimating Individual and Household Reproduction Numbers
in an Emerging Epidemic. *PLoS ONE* **2**(8): e758.
<https://doi.org/10.1371/journal.pone.0000758>

Cori, A., Ferguson, N.M., Fraser, C., Cauchemez, S. (2013) A new
framework and software to estimate time-varying reproduction numbers
during epidemics, *American Journal of Epidemiology*, **178**(9),
1505–1512. <https://doi.org/10.1093/aje/kwt133>

<hr>
