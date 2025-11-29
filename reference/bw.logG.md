# Bandwidth estimation for strictly positive distributions.

Computes bandwidth for log domain KDE using the Silverman rule.

## Usage

``` r
bw.logG(x)
```

## Arguments

- x:

  numeric vector of the data. Must be strictly positive, will be log
  transformed during estimation.

## Value

bw the optimal bandwidth.

## References

Silverman, B. W. (1986). Density estimation for statistics and data
analysis. Monographs on Statistics and Applied Probability. 26.

Wand, M. P., Marron, J. S., & Ruppert, D. (1991). Transformations in
density estimation. Journal of the American Statistical Association,
86(414), 343-353.

## Examples

``` r
bw.logG(rchisq(100,10))
#> [1] 0.1836065
```
