#' Computes Kernel Density Estimates in Log Domain.
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
#' logdensity(abs(rnorm(100)), from =.1, to=2, kernel='triangular')
#'
#'@export
logdensity <- function(x, bw = "nrd0", adjust = 1,
           kernel = "gaussian",
           weights = NULL, n = 512, from, to, cut = 3, na.rm = FALSE, ...)
  {

   #c("gaussian", "epanechnikov", "triangular", "uniform", "laplace", "logistic")

    #Taken from standard density function for consistency
    chkDots(...)
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

    if(any(!x>0)) {
      x <- x[x>0]
      nx <- length(x) # == sum(x.finite)
      warning("Non-positive values of x have been removed!")
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
      if(any(!x.finite|x<=0)) {
        weights <- weights[x.finite|x<=0]
        totMass <- sum(weights) / wsum
      } else totMass <- 1

      ## No error, since user may have wanted "sub-density"
      if (!isTRUE(all.equal(1, wsum)))
        warning("sum(weights) != 1  -- will not get true density")
    }

    n.user <-n
    dat<-x

    #need to actually do something with weights

    if (is.character(bw)) {
      if(nx < 2)
        stop("need at least 2 points to select a bandwidth automatically")
      bw <- switch(tolower(bw),
                   nrd0 = bw.nrd0(log(x)),
                   nrd = bw.nrd(log(x)),
                   ucv = bw.ucv(log(x)),
                   bcv = bw.bcv(log(x)),
                   sj = , "sj-ste" = bw.SJ(log(x), method="ste"),
                   "sj-dpi" = bw.SJ(log(x), method="dpi"),
                   stop("unknown bandwidth rule"))
    }
    if (!is.finite(bw)) stop("non-finite 'bw'")
    bw <- adjust * bw
    if (bw <= 0) stop("'bw' is not positive.")

    if (missing(from)){
      from <- max(min(x) - cut * bw, 0.0001)
      if(min(x) - cut * bw<0.0001){
        warning("Auto-range choice cut-off at 0.")
      }
    }
    if (missing(to))
      to   <- max(x) + cut * bw
    if (!is.finite(from)) stop("non-finite 'from'")
    if (!is.finite(to)) stop("non-finite 'to'")
    if (from<0) stop("negative 'from'")
    if (to<0) stop("negative  'to'")
    lo <- from - 4 * bw
    up <- to + 4 * bw

    x <- seq.int(from, to, length.out = n.user)

    #kernel = "gaussian"
    #need to do something when pass kernel name
    y<-logKDE(dat,x,bw,kernel)

    structure(list(x = x, y = y, bw = bw, n = N,
                   call=match.call(), data.name=name, has.na = FALSE),
              class="density")
  }



