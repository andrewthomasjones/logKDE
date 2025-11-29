# Optimal CV BW estimation for strictly positive distributions.

Computes least squares cross-validation (CV) bandwidth (BW) for log
domain KDE.

## Usage

``` r
bw.logCV(x, grid = 21, NB = 512)
```

## Arguments

- x:

  numeric vector of the data. Must be strictly positive, will be log
  transformed during estimation.

- grid:

  number of points used for BW selection CV grid.

- NB:

  number of points at which to estimate the KDE at during the CV loop.

## Value

bw the optimal least squares CV bandwidth.

## References

Silverman, B. W. (1986). Density estimation for statistics and data
analysis. Monographs on Statistics and Applied Probability. 26.

Stone, C. J. (1984). An asymptotically optimal window selection rule for
kernel density estimates. The Annals of Statistics, 12(4), 1285-1297.

## Examples

``` r
bw.logCV(rchisq(100,10), grid=21, NB=512)
#> [1] 0.2684358
```
