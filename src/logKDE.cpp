//[[Rcpp::plugins(cpp11)]]
//[[Rcpp::depends(Rcpp)]]
//[[Rcpp::depends(BH)]]


//'@importFrom Rcpp sourceCpp
//'@useDynLib logKDE
//'
#include <Rcpp.h>
#include <vector>
#include <string>
#include <algorithm>
#include <math.h>


#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/bind.hpp>

using namespace Rcpp;
using namespace boost::accumulators;


// [[Rcpp::export]]
double logEpanechnikov(double x){
  double bound = sqrt(5);
  double y = 3*(5-(x*x))/(20*sqrt(5)*x)*(x>(-1*bound))*(x<(1*bound));
  return y;
}

// [[Rcpp::export]]
double logGaussian(double x){
  double y=std::pow(2.0*M_PI,-0.5)*(1/x)*std::exp(-0.5*x*x);
  return y;
}

// [[Rcpp::export]]
double logLaplace(double x){
  double val = std::pow(2.0,-0.5);
  double y=(1/x)*val*std::exp(-1*val*std::abs(x));
  return y;
}

// [[Rcpp::export]]
double logLogistic(double x){
  double bound = sqrt(3);
  double tmp = std::pow(std::cosh(M_PI*x/(2*bound)),2.0);
  double y=(M_PI/(4*bound))*(1/x)*tmp;
  return y;
}

// [[Rcpp::export]]
double logTriangular(double x){
  double bound = sqrt(6);
  double y = (1/x)*(std::pow(6, -0.5)-std::abs(x))*(x>(-1*bound))*(x<(1*bound));
  return y;
}


// [[Rcpp::export]]
double uniform(double x){
  double bound = sqrt(3);
  double y=(1/(2*bound))*(x>(-1*bound))*(x<(1*bound));
  return y;
}


// [[Rcpp::export]]
double KDE(double x, std::vector<double> xi, double h){
  int n = xi.size();
  std::vector<double> tmp(n);
  std::transform(xi.begin(), xi.end(), tmp.begin(), [x,h](double xi){return uniform((x-xi)/h);});
  double dens_x = (1/(n*h))*std::accumulate(tmp.begin(), tmp.end(), 0.0);
  return(dens_x);
}


//silvermans rule of thumb for bandwidth
// [[Rcpp::export]]
double silverman(std::vector<double> x){
  accumulator_set<double, features<tag::variance(lazy)>> acc;
  for_each(x.begin(), x.end(), boost::bind<void> (boost::ref(acc), _1));
  return 1.06*sqrt(variance(acc))*std::pow((x.size()),-0.2);

}


//'@export
// [[Rcpp::export]]
std::vector<double> logKDE(const std::vector<double>& input, const std::vector<double>& support,  double h, std::string method="logUniform"){
  int n = input.size();
  int N = support.size();

  std::vector<double> xi(n);

  std::vector<double> ret(N);
  std::vector<double> ret2(N);
  std::vector<double> X(N);

  //take log of input
  std::transform(input.begin(), input.end(), xi.begin(), [](double x){return log(x);});

  //take log of support
  std::transform(support.begin(), support.end(), X.begin(), [](double x){return log(x);});

  //Rcout<< "h = " << h << std::endl;
  //do kernels
  std::transform(X.begin(), X.end(), ret.begin(), [h, xi](double X){return KDE(X, xi, h);});
  std::transform( ret.begin(), ret.end(),support.begin(), ret2.begin(), [](double d, double x) {return (1/x)*d;} );
  return ret2;
}

