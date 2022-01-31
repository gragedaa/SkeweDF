#include <Rcpp.h>
#include <vector>

using namespace Rcpp;

//[[Rcpp::export]]
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
}

//' Kolmogorov Waring P0 calculation
//' 
//' Calculates P0 of Kolmogorov Waring distribution function given parameters
//' @param a Parameter of the Kolmogorov Waring distribution function
//' @param b Parameter of the Kolmogorov Waring distribution function
//' @param theta Parameter of the Kolmogorov Waring distribution function
//' @export
// [[Rcpp::export]]
double Kolmogorov_Waring_P0_calc(double a, double b, double theta){
	int i = 0, m = 0;
	double p0 = -1, p1 = 0, p = 1 , sp = 0;
	while (fabs(p1-p0)> 0.000001)
	{
		p0 = p1;
		m++;
		i++;
		p = p * ( a-1+i) * theta / (b + i);
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
double Kolmogorov_Waring_P0(double a, double b, double theta){
  if ((b > a) & (theta > 0) & (theta < 1)){
    return(1 - (a/b));
  }
  else{
   return(Kolmogorov_Waring_P0_calc(a,b,theta));
  }
}

// [[Rcpp::export]]
double Kolmogorov_Waring_P0_hyp(double a, double b, double theta){
  double out = hypergeometric(a,1,b+1.0, theta);
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
std::vector<double> Kolmogorov_Waring(int n, double a, double b, double theta, bool use_hypergeometric = true){
  
  std::vector<double> p_values(n+1);
  
  if(use_hypergeometric){
    p_values[0] = Kolmogorov_Waring_P0_hyp(a,b,theta);
  }else{
    p_values[0] = Kolmogorov_Waring_P0(a,b,theta);
  }
    
  double sum_p = p_values[0];
  
  for(int i = 1; i < n+1; i++){
    p_values[i] = p_values[i-1] * theta * (a + i -1)/(b+i);
    sum_p += p_values[i];
  }
  
  for(int i = 0; i < n; i++){
    p_values[i] = p_values[i] / sum_p;
  }

  return(p_values);
  
}
