
<!-- README.md is generated from README.Rmd. Please edit that file -->

# didgformula

<!-- badges: start -->

<!-- badges: end -->

The R package `didgformula` implements inverse-probability weighted,
iterated conditional, and outcome regression estimators for the
difference-in-differences g-formula.

The difference-in-differences (DID) g-formula extends conventional DID methods, which target average treatment effects in the treated (ATTs), to more general causal effects defined as average counterfactual outcomes under a (possibly time-varying) exposure regime. DID g-formula also extends conventional g-formula methods, which base identification on an assumption called *exchangeability*, which requires that a sufficient set of covariates to control confounding is measured. Instead, DID g-formula relies on an assumption called *parallel trends*, standard in DID approaches, which requires that under the specified treatment regime, conditional on covariates, counterfactual outcome *trends* would have been parallel (on an appropriate scale) regardless of actual treatment received. 

## Installation

Only the development version is available so far. You can install it
from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("audreyrenson/didgformula")
```

## Examples

### Continuous outcomes

First we simulate some data under the stated assumptions:

``` r
library(didgformula)

set.seed(10)

time_periods = 5
N_obs        = 1e4
parameters   = generate_parameters(Tt=time_periods)
df           = generate_data(N=N_obs, Tt=time_periods, Beta=parameters, ylink = 'rnorm_identity')

head(df)
#>   uid U0 L0 A0 L1 A1 L2 A2 L3 A3 L4 A4 L5 A5        Y0         Y1        Y2
#> 1   1  0  0  0  0  0  0  0  0  0  1  0  0  0 -3.719838 -0.4517163 -1.284486
#> 2   2  1  0  0  0  0  1  0  0  0  0  1  1  1 -3.864279 -1.2299652 -1.214372
#> 3   3  1  1  0  1  0  0  0  1  0  0  0  0  0 -1.755012 -2.2640601 -1.644146
#> 4   4  0  0  0  1  0  1  0  0  0  0  1  0  1 -4.721246 -3.6573464 -3.067212
#> 5   5  0  1  0  0  1  0  1  0  1  1  1  0  1 -5.334336 -2.0876061 -2.514040
#> 6   6  1  0  0  0  0  0  0  1  0  0  0  0  0 -4.848701 -0.9340174 -1.354781
#>          Y3         Y4        Y5
#> 1  0.675204 -1.5812527 -1.543083
#> 2  1.049619 -3.6433462 -3.606828
#> 3  1.047379 -0.3792723 -3.277435
#> 4 -1.239129 -2.8298270 -1.840702
#> 5  2.821422 -3.1627642 -2.315439
#> 6 -1.014045 -3.8651784 -4.365492
```

We can calculate the true parameters by generating a large number of
potential outcomes under the same data-generating mechanism:

``` r
df_po = generate_data(N=N_obs*10, Tt=time_periods, Beta=parameters, ylink='rnorm_identity', potential_outcomes = TRUE)

truth = colMeans(calc_ydiffs(df_po, Tt=time_periods)) #calc_ydiffs simply takes Y_t-Y_{t-1} for t=1,...,T
truth
#> [1]  2.141214558 -0.023444726  2.386090954 -3.893338022  0.001502573
```

We can estimate this using IPTW:

``` r
estimates_iptw = iptw_pipeline(data = df, rhs_formula = '~L{t}', Tt=time_periods)
estimates_iptw
#> # A tibble: 5 x 2
#>       t estimate
#>   <int>    <dbl>
#> 1     1  2.12   
#> 2     2  0.00755
#> 3     3  2.39   
#> 4     4 -3.91   
#> 5     5  0.0128
```

ICE:

``` r
estimates_ice = ice_pipeline(data = df, inside_formula = '~L{t}', outside_formula = '~L{k}', Tt=time_periods)
estimates_ice
#> # A tibble: 5 x 2
#>       t estimate
#>   <int>    <dbl>
#> 1     1  2.12   
#> 2     2  0.00940
#> 3     3  2.39   
#> 4     4 -3.91   
#> 5     5  0.0188
```

Outcome regression:

``` r
estimates_or = or_pipeline(data = df, y_formula = '~L{t}', l_formula = '~1', Tt=time_periods, nreps=N_obs) #usually nreps should be much larger but in this example it appears fine
estimates_or
#> # A tibble: 5 x 2
#>       t estimate
#>   <int>    <dbl>
#> 1     1  2.12   
#> 2     2  0.00936
#> 3     3  2.39   
#> 4     4 -3.91   
#> 5     5  0.0184
```
