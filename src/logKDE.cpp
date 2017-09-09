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
double epanechnikov(double x){
  double bound = sqrt(5);
  double y = 3*(5-(x*x))/(20*sqrt(5))*(x>(-1*bound))*(x<(1*bound));
  return y;
}

// [[Rcpp::export]]
double gaussian(double x){
  double y=std::pow(2.0*M_PI,-0.5)*std::exp(-0.5*x*x);
  return y;
}

// [[Rcpp::export]]
double laplace(double x){
  double val = std::pow(2.0,-0.5);
  double y=val*std::exp(-1*val*std::abs(x));
  return y;
}

// [[Rcpp::export]]
double logistic(double x){
  double bound = sqrt(3);
  double tmp = std::pow(std::cosh(M_PI*x/(2*bound)),-2.0);
  double y=(M_PI/(4*bound))*tmp;
  return y;
}

// [[Rcpp::export]]
double triangular(double x){
  double bound = sqrt(6);
  double y = (std::pow(6, 0.5)-std::abs(x))*(x>(-1*bound))*(x<(1*bound));
  return y;
}


// [[Rcpp::export]]
double uniform(double x){
  double bound = sqrt(3);
  double y=(1/(2*bound))*(x>(-1*bound))*(x<(1*bound));
  return y;
}


// [[Rcpp::export]]
double KDE_e(double x, std::vector<double> xi, double h){
  int n = xi.size();
  std::vector<double> tmp(n);
  std::transform(xi.begin(), xi.end(), tmp.begin(), [x,h](double xi){return epanechnikov((x-xi)/h);});
  double dens_x = (1/(n*h))*std::accumulate(tmp.begin(), tmp.end(), 0.0);
  return(dens_x);
}

// [[Rcpp::export]]
double KDE_g(double x, std::vector<double> xi, double h){
  int n = xi.size();
  std::vector<double> tmp(n);
  std::transform(xi.begin(), xi.end(), tmp.begin(), [x,h](double xi){return gaussian((x-xi)/h);});
  double dens_x = (1/(n*h))*std::accumulate(tmp.begin(), tmp.end(), 0.0);
  return(dens_x);
}
// [[Rcpp::export]]
double KDE_lo(double x, std::vector<double> xi, double h){
  int n = xi.size();
  std::vector<double> tmp(n);
  std::transform(xi.begin(), xi.end(), tmp.begin(), [x,h](double xi){return logistic((x-xi)/h);});
  double dens_x = (1/(n*h))*std::accumulate(tmp.begin(), tmp.end(), 0.0);
  return(dens_x);
}
// [[Rcpp::export]]
double KDE_t(double x, std::vector<double> xi, double h){
  int n = xi.size();
  std::vector<double> tmp(n);
  std::transform(xi.begin(), xi.end(), tmp.begin(), [x,h](double xi){return triangular((x-xi)/h);});
  double dens_x = (1/(n*h))*std::accumulate(tmp.begin(), tmp.end(), 0.0);
  return(dens_x);
}
// [[Rcpp::export]]
double KDE_la(double x, std::vector<double> xi, double h){
  int n = xi.size();
  std::vector<double> tmp(n);
  std::transform(xi.begin(), xi.end(), tmp.begin(), [x,h](double xi){return laplace((x-xi)/h);});
  double dens_x = (1/(n*h))*std::accumulate(tmp.begin(), tmp.end(), 0.0);
  return(dens_x);
}
// [[Rcpp::export]]
double KDE_u(double x, std::vector<double> xi, double h){
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
std::vector<double> logKDE(const std::vector<double>& input, const std::vector<double>& support,  double h, std::string method="uniform"){
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

  //do kernels, this bit needs to be generalised

  if(method=="gaussian"){
    std::transform(X.begin(), X.end(), ret.begin(), [h, xi](double X){return KDE_g(X, xi, h);});
  }
  if(method=="epanechnikov"){
    std::transform(X.begin(), X.end(), ret.begin(), [h, xi](double X){return KDE_e(X, xi, h);});
  }
  if(method=="triangular"){
    std::transform(X.begin(), X.end(), ret.begin(), [h, xi](double X){return KDE_t(X, xi, h);});
  }
  if(method=="uniform"){
    std::transform(X.begin(), X.end(), ret.begin(), [h, xi](double X){return KDE_u(X, xi, h);});
  }
  if(method=="laplace"){
    std::transform(X.begin(), X.end(), ret.begin(), [h, xi](double X){return KDE_la(X, xi, h);});
  }
  if(method=="logistic"){
    std::transform(X.begin(), X.end(), ret.begin(), [h, xi](double X){return KDE_lo(X, xi, h);});
  }



  std::transform( ret.begin(), ret.end(),support.begin(), ret2.begin(), [](double d, double x) {return (1/x)*d;} );
  return ret2;
}

