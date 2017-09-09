logKDEdensity <-function(x, bw = "nrd0", adjust = 1,
           kernel = "gaussian",
           weights = NULL, window = kernel, n = 512, from, to, cut = 3, na.rm = FALSE, ...)
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

    if (missing(from))
      from <- min(x) - cut * bw
    if (missing(to))
      to   <- max(x) + cut * bw
    if (!is.finite(from)) stop("non-finite 'from'")
    if (!is.finite(to)) stop("non-finite 'to'")
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



