
<!-- README.md is generated from README.Rmd. Please edit that file -->

# didgformula

<!-- badges: start -->
<!-- badges: end -->

The R package `didgformula` implements inverse-probability weighted,
iterated conditional g-computation, and doubly robust targeted maximum
likelihood estimators for sustained intervention effects under parallel
trends assumptions.

## Installation

Only the development version is available so far. You can install it
from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("audreyrenson/didgformula")
```

## Examples

Here is a basic example using simulated data:

``` r
library(didgformula)

set.seed(10)

time_periods = 5
N_obs        = 1e4
parameters   = generate_parameters(Tt=time_periods)
df           = generate_data(N=N_obs, Tt=time_periods, Beta=parameters, ylink = 'rnorm_identity')

head(df)
#>   uid U0 L0         W0 A0 L1         W1 A1 L2          W2 A2 L3         W3 A3
#> 1   1  1  1 0.42889213  0  0 -2.7811919  0  0  0.18963549  0  0 -2.4548948  1
#> 2   2  1  0 1.53171732  0  1  1.0080488  0  1  0.06726136  0  0  1.4825164  0
#> 3   3  0  1 0.65063303  0  1 -0.3276539  1  0 -1.68579853  1  1  1.5884367  1
#> 4   4  1  0 0.92350947  0  0 -1.3764935  0  0 -0.19326923  0  0 -1.0013806  1
#> 5   5  1  1 1.40689657  0  0 -0.5103987  0  0  0.57436194  0  0 -1.0926517  0
#> 6   6  1  0 0.07771018  0  1 -0.3072492  0  0 -0.35724493  0  0  0.4300903  0
#>   L4         W4 A4 L5          W5 A5         Y0        Y1          Y2
#> 1  0 -0.2439675  1  1  2.09577923  1 -2.1574128 -1.766963  0.11996427
#> 2  0 -0.6986895  0  0 -0.08951730  0 -0.7277109 -3.525773 -1.24276996
#> 3  0 -0.7389423  1  0  0.03679411  1 -2.9713932 -4.313408 -0.04549021
#> 4  0  1.0629917  1  0 -0.85402635  1 -1.8168327 -3.201612 -0.58378946
#> 5  0 -0.3215103  0  0 -2.07722722  0 -1.0881852 -3.542649  0.45214506
#> 6  1 -0.7852357  0  0 -0.07780930  0 -1.5003276 -2.923101 -1.09694604
#>           Y3         Y4          Y5
#> 1 -2.4577156 -0.4874623  3.61649314
#> 2  0.4199241  0.2020285  1.22008428
#> 3 -0.7302781  1.1281042  1.50973592
#> 4 -0.0680779  0.8496794 -0.73340196
#> 5 -1.9242997  0.1638116  0.04507994
#> 6 -0.5373790 -0.5839260  3.11250228
```

We can calculate the true parameters by generating a large number of
potential outcomes under the same data-generating mechanism:

``` r
df_po = generate_data(N=N_obs*10, Tt=time_periods, Beta=parameters, ylink='rnorm_identity', potential_outcomes = TRUE)

truth = colMeans(calc_ydiffs(df_po, Tt=time_periods)) #calc_ydiffs simply takes Y_t-Y_{t-1} for t=1,...,T
truth
#> [1] -1.84302146  2.95303511 -0.01427988  0.92387834  0.74953343
```

We can estimate this using IPTW:

``` r
estimates_iptw = iptw_pipeline(data = df, den_formula = '~L{t}', Tt=time_periods)
estimates_iptw
#> # A tibble: 5 x 2
#>       t estimate
#>   <int>    <dbl>
#> 1     1  -1.91  
#> 2     2   3.00  
#> 3     3   0.0259
#> 4     4   0.907 
#> 5     5   0.681
```

ICE:

``` r
estimates_ice = ice_pipeline(data = df, inside_formula_t = '~L{t}', inside_formula_tmin1='~L{t-1}', outside_formula = '~L{k}', Tt=time_periods)
estimates_ice
#> # A tibble: 5 x 2
#>       t estimate
#>   <int>    <dbl>
#> 1     1  -1.91  
#> 2     2   3.00  
#> 3     3   0.0251
#> 4     4   0.908 
#> 5     5   0.680
```
