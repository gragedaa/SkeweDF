#include <Rcpp.h>
#include <vector>
#include <math.h>
#include <boost/math/special_functions.hpp>

// [[Rcpp::depends(BH)]]
using namespace Rcpp;

double hypergeometric( double a, double b, double c, double x )
{
  const double TOLERANCE = 1.0e-10;
  double term = a * b * x / c;
  double value = 1.0 + term;
  int n = 1;

  while (abs( term ) > TOLERANCE )
  {
    a++, b++, c++, n++;
    term *= a * b * x / c / n;
    value += term;
  }

  return value;
}

//' @export
// [[Rcpp::export]]
double boost_hypergeometric_2F1(double a1, double a2, double b, double theta){
  return(boost::math::hypergeometric_pFq({a1,a2}, {b}, theta));
}

//' @export
// [[Rcpp::export]]
double Kolmogorov_Waring_P0(double a1, double a2, double b, double theta){
  double out;
  
  if(theta > 1){
    out = hypergeometric(a1,a2,b+1.0, theta);
  }else{
     out = boost::math::hypergeometric_pFq({a1,a2}, {b+1.0}, theta);
  }
  
  return(1/out);
}


//' Kolmogorov Waring Factorial Moment
//'
//' Calculates P0 of Kolmogorov Waring distribution function given parameters
//' @param a1 Parameter of the Kolmogorov Waring distribution function
//' @param b1 Parameter of the Kolmogorov Waring distribution function
//' @param theta Parameter of the Kolmogorov Waring distribution function
//' @param r Index of factorial moment to calculate
//' @export
// [[Rcpp::export]]
double Kolmogorov_Waring_Moment(double a1, double b1, double theta, double r){
  double fact = 1.0;
  for(int i = 1; i<= r; i++){
    fact = fact * i;
  }
  
  double mid_theta = pow(theta,r);
  double mid_a = pow(a1,r);
  double mid_b = pow(b1+1.0,r);
  double mid = mid_a/mid_b;
  mid = mid * mid_theta *  Kolmogorov_Waring_P0(a1, 1, b1, theta);
  
  double hyp;
  if(theta > 1){
    hyp = hypergeometric(a1+r, 1.0+r, b1+1.0+r, theta);
  }else{
    hyp = boost::math::hypergeometric_pFq({a1+r, 1.0+r}, {b1+1.0+r}, theta);
  }
  
  return(fact * mid * hyp);
}


//' Kolmogorov Waring
//'
//' Calculates vector of n length of Kolmogorov distribution function given parameters
//' @param n Length of vector to be generated
//' @param a Vector of parameters of the Kolmogorov Waring distribution function
//' @param b Vector of parameters of the Kolmogorov Waring distribution function
//' @param theta Parameter of the Kolmogorov Waring distribution function
//' @export
// [[Rcpp::export]]
std::vector<double> Kolmogorov_Waring(int n, Rcpp::NumericVector a, Rcpp::NumericVector b, double theta){

  std::vector<double> p_values(n+1);


  if(a.length() > b.length()){
    for(int i = 0; i < (a.length() - b.length()); i++){
      b.push_back(0.0);
    }
  }else if(b.length() > a.length()){
    for(int i = 0; i < (b.length() - a.length()); i++){
      a.push_back(0.0);
    }
  }

  if((a.length() == 1) & (b.length() == 1)){
    p_values[0] = Kolmogorov_Waring_P0(a[0], 1, b[0], theta);
  }else{
    
    p_values[0] = Kolmogorov_Waring_P0(a[0],a[1], b[0], theta);
  }
  
  double sum_p = p_values[0];

  for(int i = 1; i < n+1; i++){
    p_values[i] = p_values[i-1] * theta;

    for(int j = 0; j < a.length(); j++){
      p_values[i] =  p_values[i] * (a[j] + i);
    }
    for(int j = 0; j < b.length();j++){
      p_values[i] =  p_values[i] / (b[j] + i + 1);
    }

    sum_p += p_values[i];
  }

  for(int i = 0; i < n+1; i++){
    p_values[i] = p_values[i] / sum_p;
  }

  return(p_values);

}

//' Kolmogorov Waring no trunc
//'
//' Calculates vector of n length of Kolmogorov distribution function given parameters
//' @param n Length of vector to be generated
//' @param a Vector of parameters of the Kolmogorov Waring distribution function
//' @param b Vector of parameters of the Kolmogorov Waring distribution function
//' @param theta Parameter of the Kolmogorov Waring distribution function
//' @export
// [[Rcpp::export]]
std::vector<double> Kolmogorov_Waring_no_trunc(int n, Rcpp::NumericVector a, Rcpp::NumericVector b, double theta){
  
  std::vector<double> p_values(n+1);
  
  
  if(a.length() > b.length()){
    for(int i = 0; i < (a.length() - b.length()); i++){
      b.push_back(0.0);
    }
  }else if(b.length() > a.length()){
    for(int i = 0; i < (b.length() - a.length()); i++){
      a.push_back(0.0);
    }
  }
  
  if((a.length() == 1) & (b.length() == 1)){
    double ratio = a[0] / b[0];
    p_values[0] = 1 - ratio;
  }else{
    
    p_values[0] = Kolmogorov_Waring_P0(a[0],a[1], b[0], theta);
  }
  
  double sum_p = p_values[0];
  
  for(int i = 1; i < n+1; i++){
    p_values[i] = p_values[i-1] * theta;
    
    for(int j = 0; j < a.length(); j++){
      p_values[i] =  p_values[i] * (a[j] + i);
    }
    for(int j = 0; j < b.length();j++){
      p_values[i] =  p_values[i] / (b[j] + i + 1);
    }
    
    sum_p += p_values[i];
  }
  
  
  return(p_values);
  
}
