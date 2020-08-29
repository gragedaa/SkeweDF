#include <Rcpp.h>
#include <vector>
#include <iostream>
            
//using namespace Rcpp;


// [[Rcpp::export]]
double RGHD_P0_calc_function(int k, int m, Rcpp::NumericVector r, Rcpp::NumericVector q){

  //Rcpp::NumericVector rc = clone(r);
  //Rcpp::NumericVector qc = clone(q);
  
  std::vector<double> rc(m);
  std::vector<double> qc(m);
  //std::vector<double> rc = as<std::vector<double>>(r);
  //std::vector<double> qc = as<std::vector<double>>(q);
	
	/**
  for(int i = 0; i < m; i++){
	  rc[i] += k;
	  qc[i] += k;
  }
	*/
	
	for(int i = 0; i < m; i++){
	  rc[i] = r[i] + (k * (i+1));
	  qc[i] = q[i] + (k * (i+1));
	}
	
  double r_product = 1;
  double q_product = 1;
	
	for(int i = 0; i < m; i++){
	  	r_product *= rc[i];
	  	q_product *= qc[i];
  	}
	
	return(r_product / q_product);
}

/*
RGHD_P0_calc_pi - Calculates Pi with k from 0 to y-1 for RGHD function given a set of parameters

int y - index of Sigma indicating upper limit of Pi
int m - number indicating (number of parameters / 2)
Rcpp::NumericVector r - vector of parameters with length m in the numerator
Rcpp::NumericVector q - vector of parameters with length m in the denominator

*/
            
// [[Rcpp::export]]
double RGHD_P0_calc_pi(int y, int m, Rcpp::NumericVector r, Rcpp::NumericVector q){
  double cumulative_product = 1;
  
  for(int i = 0; i < y; i++){
    cumulative_product *= RGHD_P0_calc_function(i,m,r,q);
  }
  return(cumulative_product);
}

/*
 RGHD_P0_calc - Calculates P0 with Sigma lower limit of 1 and upper limit of argument sigma_upper
 
 int sigma_upper - upper limit of Sigma
 int m - number indicating (number of parameters / 2)
 Rcpp::NumericVector r - vector of parameters with length m in the numerator
 Rcpp::NumericVector q - vector of parameters with length m in the denominator
 
 */

// [[Rcpp::export]]
double RGHD_P0_calc(int sigma_upper, int m, Rcpp::NumericVector r, Rcpp::NumericVector q){
  double sum = 0;
  
  for(int i = 0; i < sigma_upper; i++){// this loop with calculate and store Pi where y = (1,2,...,sigma_upper)
    sum +=  RGHD_P0_calc_pi(i+1,m,r,q);
  }
  
  return(1 / (1 + sum)); // return P0
}

//' 2m-RGHD Distribution Function P0
//' 
//' Calculates P0 given a set of parameters
//' @param m Parameter of the 2m-RGHD function, this defines number of r and q parameters of the function
//' @param r R vector containing r parameters from 1:m
//' @param q R vector containing q parameters from 1:m
//' @export
// [[Rcpp::export]]
double RGHD_P0(int m, Rcpp::NumericVector r, Rcpp::NumericVector q){
  double delta = 0.0001; // this is the minimum required difference from one P0 calculation to the next
  double diff = 1;
  double p0 = 0.0;
  double prev_p0 = 1.0;
  int y = 1;
  
  while(diff > delta){
    p0 = RGHD_P0_calc(y,m,r,q);
    diff = prev_p0 - p0;
    diff = (diff * diff) / diff;
    prev_p0 = p0;
    y++;
  }

  return(p0);
}


//' 2m-RGHD Distribution Function
//' 
//' Returns doubly truncated vector of 2m-RGHD function values where input is 1-J
//' @param J Length of vector to be generated
//' @param m Parameter of the 2m-RGHD function, this defines number of r and q parameters of the function
//' @param r R vector containing r parameters from 1:m
//' @param q R vector containing q parameters from 1:m
//' @param P0_iter Integer indicating number of iterations to use for calculation of P0, increasing this parameter will increase accuracy of P0
//' @param P0_included Boolean used to include P0 in vector or not
//' @export
// [[Rcpp::export]]
std::vector<double> RGHD(int J, int m, Rcpp::NumericVector r, Rcpp::NumericVector q, int P0_iter=100, bool P0_included = false){
	std::vector<double> p(J+1);
  
	double sum_p = 0;
	//p[0] = RGHD_P0(m, r, q); // calculate P0
	p[0] = RGHD_P0_calc(P0_iter,m, r, q); // calculate P0
	
	Rcpp::NumericVector rc = clone(r);
	Rcpp::NumericVector qc = clone(q);
	//std::vector<double> rc(m);
	//std::vector<double> qc(m);
	
	for(int i = 1; i < J+1; i++){// loop will calculate P1 - PJ
		
		double r_product = 1;
	  double q_product = 1;
	  
	  for(int j = 0; j < m; j++){ // increment all parameters by 1 then get the product 
	    rc[j]++;
	    r_product *= rc[j];
	  }
	  for(int j = 0; j < m; j++){ // increment all parameters by 1 then get the product 
	    qc[j]++;
	    q_product *= qc[j];
	  }
	  
	  p[i] = p[i-1] * (r_product / q_product); // recursive calculation
		
		sum_p += p[i]; // calculate sum from P1 - PJ, needed to right side truncation
	}
	
	sum_p += p[0];
	
	if(P0_included){
	  std::vector<double> output(J+1);
	  for(int i = 0; i < J+1; i++){// this loop will interate through P1 - PJ normalizing each by the previously calculated sum, starting at P1 is the left side truncation and the normalization is the right side truncation
	    output[i] =  p[i] / sum_p;
	  }
	  return(output);
	}else{
	  std::vector<double> output(J);
	  sum_p -= p[0];
	  for(int i = 0; i < J; i++){// this loop will interate through P1 - PJ normalizing each by the previously calculated sum, starting at P1 is the left side truncation and the normalization is the right side truncation
	    output[i] =  p[i+1] / sum_p;
	  }
	  return(output);
	}
}
            
