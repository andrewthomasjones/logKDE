# Kernel Density Estimates of strictly positive distributions using FFT.

The function `logdensity_fft` computes kernel density estimates (KDE) of
strictly positive distributions by performing the KDE via fast fourier
transform utilizing the `fft` function. The syntax and function
structure is largely borrowed from the function `density` in package
stats.

## Usage

``` r
logdensity_fft(
  x,
  bw = "nrd0",
  adjust = 1,
  kernel = "gaussian",
  weights = NULL,
  n = 512,
  from,
  to,
  cut = log(3),
  na.rm = FALSE
)
```

## Arguments

- x:

  the data from which the estimate is to be computed.

- bw:

  the smoothing bandwidth to be used. Can also be can also be a
  character string giving a rule to choose the bandwidth. Like `density`
  defaults to "nrd0". All options in
  [`help(bw.nrd)`](https://rdrr.io/r/stats/bandwidth.html) are available
  as well as `"bw.logCV"` and `"bw.logG"`.

- adjust:

  the bandwidth used is actually `adjust*bw`.

- kernel:

  a character string giving the smoothing kernel to be used. Choose from
  "gaussian", "epanechnikov", "triangular", "uniform", "laplace" and
  "logistic". Default value is "gaussian".

- weights:

  numeric vector of non-negative observation weights of the same length
  as `x`.

- n:

  the number of equally spaced points at which the density is to be
  estimated. Note that these are equally spaced in the log domain for
  `logdensity_fft`, and thus on a log scale when transformed back to the
  original domain.

- from, to:

  the left and right-most points of the grid at which the density is to
  be estimated; the defaults are cut \* bw outside of range(x).

- cut:

  by default, the values of from and to are cut bandwidths beyond the
  extremes of the data

- na.rm:

  logical; if TRUE, missing values are removed from x. If FALSE any
  missing values cause an error.

## Value

An object with class "density". See
[`help(density)`](https://rdrr.io/r/stats/density.html) for details.

## References

Charpentier, A., & Flachaire, E. (2015). Log-transform kernel density
estimation of income distribution. L'Actualite economique, 91(1-2),
141-159.

Cooley, J. W., & Tukey, J. W. (1965). An algorithm for the machine
calculation of complex Fourier series. Mathematics of computation,
19(90), 297-301.

Wand, M. P., Marron, J. S., & Ruppert, D. (1991). Transformations in
density estimation. Journal of the American Statistical Association,
86(414), 343-353.

## See also

[`density`](https://rdrr.io/r/stats/density.html),
[`plot.density`](https://rdrr.io/r/stats/plot.density.html),
[`logdensity`](http://andrewthomasjones.com.au/logKDE/reference/logdensity.md),
[`bw.nrd`](https://rdrr.io/r/stats/bandwidth.html),
[`bw.logCV`](http://andrewthomasjones.com.au/logKDE/reference/bw.logCV.md),
[`bw.logG`](http://andrewthomasjones.com.au/logKDE/reference/bw.logG.md).

## Examples

``` r
logdensity_fft(abs(rnorm(100)), from =0.01, to= 2.5, kernel = 'logistic')
#> 
#> Call:
#>  logdensity_fft(x = abs(rnorm(100)), kernel = "logistic", from = 0.01,     to = 2.5)
#> 
#> Data: abs(rnorm(100)) (100 obs.);    Bandwidth 'bw' = 0.3077
#> 
#>        x                 y          
#>  Min.   :0.01000   Min.   :0.04666  
#>  1st Qu.:0.03976   1st Qu.:0.45039  
#>  Median :0.15812   Median :0.61538  
#>  Mean   :0.45254   Mean   :0.67347  
#>  3rd Qu.:0.62872   3rd Qu.:0.74895  
#>  Max.   :2.50000   Max.   :3.53755  
```
