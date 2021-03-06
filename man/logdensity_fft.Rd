% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logKDE.R
\name{logdensity_fft}
\alias{logdensity_fft}
\title{Kernel Density Estimates of strictly positive distributions using FFT.}
\usage{
logdensity_fft(x, bw = "nrd0", adjust = 1, kernel = "gaussian",
  weights = NULL, n = 512, from, to, cut = log(3), na.rm = FALSE)
}
\arguments{
\item{x}{the data from which the estimate is to be computed.}

\item{bw}{the smoothing bandwidth to be used. Can also be can also be a character string giving a rule to choose the bandwidth. Like \code{density} defaults to "nrd0". All options in  \code{help(bw.nrd)} are available as well as \code{"bw.logCV"} and \code{"bw.logG"}.}

\item{adjust}{the bandwidth used is actually \code{adjust*bw}.}

\item{kernel}{a character string giving the smoothing kernel to be used. Choose from "gaussian", "epanechnikov", "triangular", "uniform", "laplace" and "logistic". Default value is "gaussian".}

\item{weights}{numeric vector of non-negative observation weights of the same length as \code{x}.}

\item{n}{the number of equally spaced points at which the density is to be estimated. Note that these are equally spaced in the log domain for \code{logdensity_fft}, and thus on a log scale when transformed back to the original domain.}

\item{from, to}{the left and right-most points of the grid at which the density is to be estimated; the defaults are cut * bw outside of range(x).}

\item{cut}{by default, the values of from and to are cut bandwidths beyond the extremes of the data}

\item{na.rm}{logical; if TRUE, missing values are removed from x. If FALSE any missing values cause an error.}
}
\value{
An object with class "density". See \code{help(density)} for details.
}
\description{
The function \code{logdensity_fft} computes kernel density estimates (KDE) of strictly positive distributions by performing the KDE via fast fourier transform utilizing the \code{fft} function. The syntax and function structure is largely borrowed from the function \code{density} in package {stats}.
}
\examples{
logdensity_fft(abs(rnorm(100)), from =0.01, to= 2.5, kernel = 'logistic')

}
\references{
Charpentier, A., & Flachaire, E. (2015). Log-transform kernel density estimation of income distribution. L'Actualite economique, 91(1-2), 141-159.

Cooley, J. W., & Tukey, J. W. (1965). An algorithm for the machine calculation of complex Fourier series. Mathematics of computation, 19(90), 297-301.

Wand, M. P., Marron, J. S., & Ruppert, D. (1991). Transformations in density estimation. Journal of the American Statistical Association, 86(414), 343-353.
}
\seealso{
\code{\link{density}}, \code{\link{plot.density}}, \code{\link{logdensity}}, \code{\link{bw.nrd}}, \code{\link{bw.logCV}},  \code{\link{bw.logG}}.
}
