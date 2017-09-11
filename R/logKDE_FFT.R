#' Computes Kernel Density Estimates in Log Domain using FFT
#'
#' @param x the data from which the estimate is to be computed.
#' @param bw bandwidth.
#' @param adjust the bandwidth used is actually adjust*bw.
#' @param kernel choose from "gaussian", "epanechnikov", "triangular", "uniform", "laplace" and "logistic".
#' @param weights numeric vector of non-negative observation weights.
#' @param n the number of equally spaced points at which the density is to be estimated.
#' @param from,to the left and right-most points of the grid at which the density is to be estimated; the defaults are cut * bw outside of range(x).
#' @param cut by default, the values of from and to are cut bandwidths beyond the extremes of the data
#' @return Density object. See help(density).
#' @examples
#' logdensity_fft(abs(rnorm(100)), kernel = 'logistic')
#' logdensity_fft(abs(rnorm(100)), from =0.001, to= 2.5, kernel = 'logistic')
#'
#'@export
logdensity_fft <-
  function(x, bw = "nrd0", adjust = 1,
           kernel = "gaussian",
           weights = NULL, give.Rkern = FALSE,
           n = 512, from, to, cut = log(3), na.rm = FALSE, ...)
  {


    chkDots(...)

    if(give.Rkern)
      ##-- sigma(K) * R(K), the scale invariant canonical bandwidth:
      return(switch(kernel,
                    gaussian = 1/(2*sqrt(pi)),
                    uniform = sqrt(3)/6,
                    triangular = sqrt(6)/9,
                    epanechnikov = 3/(5*sqrt(5)),
                    laplace=2^(-0.5),
                    logistic=pi/(4*sqrt(3))*cosh(pi/(2*sqrt(3)))^-2 #????

      ))


    #x<-log(x)
    if (!is.numeric(x))
      stop("argument 'x' must be numeric")
    name <- deparse(substitute(x))
    x <- as.vector(x)
    x.na <- is.na(x)
    if (any(x.na)) {
      if (na.rm) x <- x[!x.na]
      else stop("'x' contains missing values")
    }
    N <- nx <- as.integer(length(x))
    if(is.na(N)) stop("invalid value of length(x)")
    x.finite <- is.finite(x)
    if(any(!x.finite)) {
      x <- x[x.finite]
      nx <- length(x) # == sum(x.finite)
    }

    if(any(x<=0)) {
      x <- x[x>0]
      nx <- length(x) # == sum(x.finite)
      warning("Non-positive values of x have been removed!")
    }

    if(!missing(from)){
      if (from<0) stop("negative 'from'")
      from<-log(from)
    }
    if(!missing(to)){
      if (to<0) stop("negative  'to'")
      to<-log(to)
    }
    ## Handle 'weights'
    if(is.null(weights))  {
      weights <- rep.int(1/nx, nx)
      totMass <- nx/N
    }
    else {
      if(length(weights) != N)
        stop("'x' and 'weights' have unequal length")
      if(!all(is.finite(weights)))
        stop("'weights' must all be finite")
      if(any(weights < 0))
        stop("'weights' must not be negative")
      wsum <- sum(weights)
      if(any(!x.finite)) {
        weights <- weights[x.finite]
        totMass <- sum(weights) / wsum
      } else totMass <- 1

      ## No error, since user may have wanted "sub-density"
      if (!isTRUE(all.equal(1, wsum)))
        warning("sum(weights) != 1  -- will not get true density")
    }

    n.user <- n
    n <- max(n, 512)
    if (n > 512) n <- 2^ceiling(log2(n)) #- to be fast with FFT


    x<-log(x)
    if (is.character(bw)) {
      if(nx < 2)
        stop("need at least 2 points to select a bandwidth automatically")
      bw <- switch(tolower(bw),
                   nrd0 = bw.nrd0(x),
                   nrd = bw.nrd(x),
                   ucv = bw.ucv(x),
                   bcv = bw.bcv(x),
                   sj = , "sj-ste" = bw.SJ(x, method="ste"),
                   "sj-dpi" = bw.SJ(x, method="dpi"),
                   stop("unknown bandwidth rule"))
    }
    if (!is.finite(bw)) stop("non-finite 'bw'")
    bw <- adjust * bw
    if (bw <= 0) stop("'bw' is not positive.")

    if (missing(from)){
      from <- max(min(x) - cut * bw, log(0.01))
      if(min(x) - cut * bw<log(0.01)){
        warning("Auto-range choice cut-off at 0.")
      }
    }
    if (missing(to))
      to   <- max(x) + cut * bw
    if (!is.finite(from)) stop("non-finite 'from'")
    if (!is.finite(to)) stop("non-finite 'to'")
    lo <- from - 4 * bw
    up <- to + 4 * bw
    ## This bins weighted distances
    y <- .Call(stats:::C_BinDist, x, weights, lo, up, n) * totMass

    kords <- seq.int(0, 2*(up-lo), length.out = 2L * n)
    kords[(n + 2):(2 * n)] <- -kords[n:2]
    kords <- switch(kernel,
                    gaussian = dnorm(kords, sd = bw),
                    ## In the following, a := bw / sigma(K0), where
                    ##	K0() is the unscaled kernel below
                    uniform = {
                      a <- bw*sqrt(3)
                      ifelse(abs(kords) < a, .5/a, 0) },
                    triangular = {
                      a <- bw*sqrt(6) ; ax <- abs(kords)
                      ifelse(ax < a, (1 - ax/a)/a, 0) },
                    epanechnikov = {
                      a <- bw*sqrt(5) ; ax <- abs(kords)
                      ifelse(ax < a, 3/4*(1 - (ax/a)^2)/a, 0) },
                    laplace = 2^(-0.5)*exp(-1*2^(-0.5)*abs(kords)/bw),
                    logistic = {
                      a <- 2*bw*sqrt(3); ax <- abs(kords)
                      (pi/2*a)*(cosh(pi*ax/a)^(-2))}

    )





    kords <- fft( fft(y)* Conj(fft(kords)), inverse=TRUE)
    kords <- pmax.int(0, Re(kords)[1L:n]/length(y))
    xords <- seq.int(lo, up, length.out = n)
    x <- (seq.int(from, to, length.out = n.user))
    structure(list(x = exp(x), y = approx(xords, kords, x)$y*(1/exp(x)), bw = bw, n = N,
                   call=match.call(), data.name=name, has.na = FALSE),
              class="density")
  }
