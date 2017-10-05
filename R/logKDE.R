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
#' @param na.rm standard deal
#' @return Density object. See help(density).
#' @examples
#' logdensity(abs(rnorm(100)), from =.1, to=2, kernel='triangular')
#'
#'@export
logdensity <- function(x, bw = "nrd0", adjust = 1,
           kernel = "gaussian",
           weights = NULL, n = 512, from, to, cut = 3, na.rm = FALSE)
  {

   #c("gaussian", "epanechnikov", "triangular", "uniform", "laplace", "logistic")

    #Taken from standard density function for consistency
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

    if(!missing(from)){
      if (from<0) stop("negative 'from'")
     }
    if(!missing(to)){
      if (to<0) stop("negative  'to'")
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
                   nrd0 = stats::bw.nrd0(log(x)),
                   nrd = stats::bw.nrd(log(x)),
                   ucv = stats::bw.ucv(log(x)),
                   bcv = stats::bw.bcv(log(x)),
                   logg = bw.logG(x),
                   logcv = bw.logCV(x),
                    "sj-ste" = stats::bw.SJ(log(x), method="ste"),
                   "sj-dpi" = stats::bw.SJ(log(x), method="dpi"),
                   stop("unknown bandwidth rule"))
    }
    if (!is.finite(bw)) stop("non-finite 'bw'")
    bw <- adjust * bw
    if (bw <= 0) stop("'bw' is not positive.")

    if (missing(from)){
      from <- max(min(x) - cut * bw, (0.01))
      if(min(x) - cut * bw<(0.01)){
        warning("Auto-range choice cut-off at 0.")
      }
    }
    if (missing(to))
      to   <- max(x) + cut * bw
    if (!is.finite(from)) stop("non-finite 'from'")
    if (!is.finite(to)) stop("non-finite 'to'")
    x <- seq.int(from, to, length.out = n.user)

    y<-logKDE(dat,x,bw,kernel)

    structure(list(x = x, y = y, bw = bw, n = N,
                   call=match.call(), data.name=name, has.na = FALSE),
              class="density")
  }


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
#' @param na.rm standard dea
#' @return Density object. See help(density).
#' @examples
#' logdensity_fft(abs(rnorm(100)), from =0.01, to= 2.5, kernel = 'logistic')
#'
#'@export
logdensity_fft <-
  function(x, bw = "nrd0", adjust = 1,
           kernel = "gaussian",
           weights = NULL, n = 512, from, to, cut = log(3), na.rm = FALSE)
  {

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
                   nrd0 = stats::bw.nrd0(x),
                   nrd = stats::bw.nrd(x),
                   ucv = stats::bw.ucv(x),
                   bcv = stats::bw.bcv(x),
                   logg = bw.logG(exp(x)),
                   logcv = bw.logCV(exp(x)),
                   sj = , "sj-ste" = stats::bw.SJ(x, method="ste"),
                   "sj-dpi" = stats::bw.SJ(x, method="dpi"),
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
    y <- BinDist( x, weights, lo, up, n) * totMass

    kords <- seq.int(0, 2*(up-lo), length.out = 2L * n)
    kords[(n + 2):(2 * n)] <- -kords[n:2]
    kords <- switch(kernel,
                    gaussian = stats::dnorm(kords, sd = bw),
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
                    laplace = (1/bw)*(sqrt(2)/2)*exp(-1*2^(0.5)*abs(kords)/bw),
                    logistic = {
                      a <- 2*bw*sqrt(3); ax <- abs(kords)
                      (pi/2*a)*(cosh(pi*ax/a)^(-2))}

    )





    kords <- stats::fft( stats::fft(y)* Conj(stats::fft(kords)), inverse=TRUE)
    kords <- pmax.int(0, Re(kords)[1L:n]/length(y))
    xords <- seq.int(lo, up, length.out = n)
    x <- (seq.int(from, to, length.out = n.user))
    structure(list(x = exp(x), y = exp(log(stats::approx(xords, kords, x)$y)- x), bw = bw, n = N,
                   call=match.call(), data.name=name, has.na = FALSE),
              class="density")
  }


#replaces y <- .Call(stats:::C_BinDist, x, weights, lo, up, n)
BinDist <- function(x, w, lo, hi, n){
  y<-array(0,2*n)
  ixmin = 0
  ixmax = n - 2
  xdelta = (hi - lo) / (n - 1)
  for(i in 0:(length(x)-1)){
    if(is.finite(x[i+1])){
      xpos <- (x[i+1] - lo) / xdelta;
      ix <- floor(xpos);
      fx <- xpos - ix;
      wi <- w[i+1];

      if(ixmin <= ix && ix <= ixmax) {
        y[ix+1] <- y[ix+1] + (1 - fx) * wi;
        y[ix + 2] <- y[ix + 2] + fx * wi;
      }else if(ix == -1) {y[1] <- y[1] + fx * wi}
       else if(ix == ixmax+1){ y[ix+1]  <- y[ix+1] + (1 - fx) * wi}
      }
    }
  return(y)
}




#' Computes least squares bandwidth for log domain KDE using modified silverman rule.
#'
#' @param x the data from which the estimate is to be computed.
#' @return h the optimal bandwidth
#' @examples
#' bw.logG(rchisq(100,10))
#'
#'@export
bw.logG<-function(x){
 s<-stats::sd(x)
 n<-length(x)
 h = ((16*exp(0.25*s^2))/(s^4 + 4*s^2 + 12))^(0.2)*(s/(n^0.2))
 if(!is.finite(h)){
   h<-stats::bw.nrd0(log(x))
 }
 return(h)
}



#' Computes least squares CV bandwidth for log domain KDE.
#'
#' @param y the data from which the estimate is to be computed.
#' @param grid number of points for bw grid
#' @param NB number of points to estimate the KDE at during the CV loop
#' @return bw the optimal CV bandwidth
#' @examples
#' bw.logCV(rchisq(100,10), grid=21, NB=512)
#'
#'@export
bw.logCV<-function(y,  grid=21, NB=512){
  fit<-logdensity_fft(y)
  n<-length(y)
  ### Bandwidth
  HH <- 2^seq(-5,2,length.out = grid)*fit$bw
  ## Storage for CV
  CVSTORE <- c()

  for (hh in 1:grid) {
    LD <- logdensity_fft(y,bw=HH[hh],n=NB)
    FF <- stats::approxfun(LD$x,LD$y^2)

    #CV <- integrate(FF,min(LD$x),max(LD$x))$value

    CV = tryCatch({
      x_seq<-seq(min(LD$x),max(LD$x),length.out=NB)
      pracma::trapz(x_seq, FF(x_seq))
      #integrate(FF,min(LD$x),max(LD$x), subdivisions = 100)$value
    }, warning = function(w) {
      NA
    }, error = function(e) {
      NA
    }, finally = {
      NA
    })

    ### Compute CV for given bandwidth and NB
    CVlist<-array(0,n)
    for (x in 1:n) {
      FF <- stats::approxfun(LD$x,logdensity_fft(y[-x],bw=HH[hh],from = min(LD$x),to = max(LD$x),n=NB)$y)
      temp<-FF(y[x])
      #print(temp)
      CVlist[x] <- temp
    }
    CV <- CV - 2*mean(CVlist, na.rm=T)
    CVSTORE[hh] <- CV
  }

  ### Get the optimal BW
  return(HH[which.min(CVSTORE)])
}



#' @import stats
