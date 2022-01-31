#include <Rcpp.h>
#include <vector>
#include <iostream>
//#include <boost/math/distributions/hypergeometric.hpp>
//#include <boost/math/distributions/negative_binomial.hpp>
//#include <gsl/specfunc/gsl_sf_hyperg.h>
#include <4PKolmogorovWaring.cpp>


using namespace Rcpp;

/*
double hypergeometric( double a, double b, double c, double x )
{
  const double TOLERANCE = 1.0e-10;
  double term = a * b * x / c;
  double value = 1.0 + term;
  int n = 1;
  
  while ( abs( term ) > TOLERANCE )
  {
    a++, b++, c++, n++;
    term *= a * b * x / c / n;
    value += term;
  }
  
  return value;
}*/

//' Kolmogorov Waring P0 calculation
//' 
//' Calculates P0 of Kolmogorov Waring distribution function given parameters
//' @param a Parameter of the Kolmogorov Waring distribution function
//' @param b Parameter of the Kolmogorov Waring distribution function
//' @param theta Parameter of the Kolmogorov Waring distribution function
//' @export
// [[Rcpp::export]]
double Generalized_Kolmogorov_Waring_P0_calc(double a1, double a2, double b, double theta){
	int i = 0, m = 0;
	double p0 = -1, p1 = 0, p = 1 , sp = 0;
	while (fabs(p1-p0)> 0.000001)
	{
		p0 = p1;
		m++;
		i++;
		p = p * (a1-1+i) * (a2-1+i) * theta / (i*(b + i));
		sp = p + sp;
		p1 = 1 / (1 + sp);
	}
	return(p1);
}


//' Kolmogorov Waring P0
//' 
//' Calculates P0 of Kolmogorov Waring distribution function given parameters. Approximation is used if parameters meet a specific criteria.
//' @param a Parameter of the Kolmogorov Waring distribution function
//' @param b Parameter of the Kolmogorov Waring distribution function
//' @param theta Parameter of the Kolmogorov Waring distribution function
//' @export
// [[Rcpp::export]]
double Generalized_Kolmogorov_Waring_P0(double a1, double a2, double b, double theta){
  return(P4Kolmogorov_Waring_P0_calc(a1,a2,b,theta));
}

// [[Rcpp::export]]
double Generalized_Kolmogorov_Waring_P0_hyp(double a1, double a2, double b, double theta){
  double out = hypergeometric(a1,a2,b+1.0, theta);
  return(1/out);
}

//' Kolmogorov Waring
//' 
//' Calculates vector of n length of Kolmogorov distribution function given parameters
//' @param n Length of vector to be generated
//' @param a Parameter of the Kolmogorov Waring distribution function
//' @param b Parameter of the Kolmogorov Waring distribution function
//' @param theta Parameter of the Kolmogorov Waring distribution function
//' @export
// [[Rcpp::export]]
std::vector<double> Generalized_Kolmogorov_Waring(int n, Rcpp::NumericVector a, Rcpp::NumericVector b, double theta){
  
  std::vector<double> p_values(n+1);

  
  if(a.length() > 1){
    p_values[0] = P4Kolmogorov_Waring_P0_calc(a[0],a[1], b[0], theta);
  }else{
    p_values[0] = Kolmogorov_Waring_P0_hyp(a[0], b[0], theta);
  }
    
  if(a.length() > b.length()){
    for(int i = 0; i < (a.length() - b.length()); i++){
      b.push_back(0.0);
    }
  }else if(b.length() > a.length()){
    for(int i = 0; i < (b.length() - a.length()); i++){
      a.push_back(1.0);
    }
  }
  
  double sum_p = p_values[0];
 
  for(int i = 1; i < n+1; i++){
    p_values[i] = p_values[i-1] * theta;
    
    for(int j = 0; j < a.length(); j++){
      p_values[i] =  p_values[i] * (a[j] + i - 1);
    }
    for(int j = 0; j < b.length();j++){
      p_values[i] =  p_values[i] / (b[j] + i);
    }
    
    sum_p += p_values[i];
  }
 
  for(int i = 0; i < n+1; i++){
    p_values[i] = p_values[i] / sum_p;
  }
  
  return(p_values);
  
}
