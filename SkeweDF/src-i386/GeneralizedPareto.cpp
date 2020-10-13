#include <Rcpp.h>
#include <vector>
#include <iostream>

//using namespace Rcpp;


// [[Rcpp::export]]
double Generalized_Pareto_calc_function(int m, double c, double b, double rho){
  return(1 + ((c-1) / pow((m+b), rho)));
}

// [[Rcpp::export]]
double Generalized_Pareto_calc_pi(double k, double theta, double c, double b, double rho){
  double pi = 1;
  
  for(int m = 0; m < k; m++){
    pi = pi * Generalized_Pareto_calc_function(m,c,b,rho);
  }
  
  double output = pow(theta,k) / pow((k+b), rho);
  output = output * pi;
  
  return(output);
}

//' Generalized Pareto Distribution Function P0 with defined delta
//' 
//' Returns P0 of Generalized Pareto given a parameters theta, c, b , and rho
//' @param theta Parameter of the Generalized Pareto function
//' @param c Parameter of the Generalized Pareto function
//' @param b  Parameter of the Generalized Pareto function
//' @param rho Parameter of the Generalized Pareto function
//' @param delta Value of difference between iterations in order to output a result. Decreasing this parameter will increase accuracy of P0. Delta > 0
//' @export
// [[Rcpp::export]]
double Generalized_Pareto_calc_P0_delta(double theta, double c, double b, double rho, double delta){
  
  double prev_sum = Generalized_Pareto_calc_pi(1,theta,c,b,rho);
  double cur_sum = prev_sum + Generalized_Pareto_calc_pi(2,theta,c,b,rho);
  
  int n = 3;
  
  while(cur_sum - prev_sum > delta){
    prev_sum = cur_sum;
    cur_sum += Generalized_Pareto_calc_pi(n,theta,c,b,rho);
    n++;
  }
  
  return(cur_sum + 1);
}

//' Generalized Pareto Distribution Function P0 with defined number of iterations
//' 
//' Returns P0 of Generalized Pareto given a parameters theta, c, b , and rho
//' @param theta Parameter of the Generalized Pareto function
//' @param c Parameter of the Generalized Pareto function
//' @param b  Parameter of the Generalized Pareto function
//' @param rho Parameter of the Generalized Pareto function
//' @param iter Number of iterations to be performed for summation calcuation. Increasing this parameter will increase accuracy of P0
//' @export
// [[Rcpp::export]]
double Generalized_Pareto_calc_P0_iter(double theta, double c, double b, double rho, int iter){
  
  double cur_sum = Generalized_Pareto_calc_pi(1,theta,c,b,rho);
  
  for(int i = 2; i <= iter; i++){
    cur_sum += Generalized_Pareto_calc_pi(i,theta,c,b,rho);
  }
  
  return(cur_sum + 1);
}

//' Generalized Pareto Distribution Function
//' 
//' Returns vector of length k of Generalized Pareto given a parameters theta, c, b , and rho
//' @param k Length of vector to be generated
//' @param theta Parameter of the Generalized Pareto function
//' @param c Parameter of the Generalized Pareto function
//' @param b  Parameter of the Generalized Pareto function
//' @param rho Parameter of the Generalized Pareto function
//' @export
// [[Rcpp::export]]
std::vector<double> Generalized_Pareto(int k, double theta, double c, double b, double rho){
 
  double p0 = Generalized_Pareto_calc_P0_iter(theta,c,b,rho,100);
  std::vector<double> p(k);

  double sum_p = 0;
  
  for(int i = 1; i <= k; i++){
    p[i-1] = (Generalized_Pareto_calc_pi(i, theta, c, b, rho) / p0);
    sum_p += p[i-1];
  }
  
  for(int i = 0; i < k; i++){
    p[i] /= sum_p;
  }
  
  return(p);
}
