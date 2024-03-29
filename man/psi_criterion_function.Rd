% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/psi_criterion_function.R
\name{psi_criterion_function}
\alias{psi_criterion_function}
\title{Psi Criterion given a function}
\usage{
psi_criterion_function(
  params,
  data,
  model_fn,
  pmf_weight = 0,
  weighted_rt = FALSE,
  left_trunc = 1,
  right_trunc = left_trunc + length(data) - 1
)
}
\arguments{
\item{params}{Vector of parameters for model_fn, not including n. For example, for Generalized_Yule(n, rho, alpha), params will be c(rho, alpha)}

\item{data}{Vector of observed values}

\item{model_fn}{Function of theoretical model to be used. For example, for Generalized_Yule(n, rho, alpha), model_fn <- Generalied_Yule}

\item{pmf_weight}{Numeric of weight given to probability mass function for generation of Psi Criterion. For example, if pmf_weight <- 0.5, 50 percent of the Psi Criterion value will be attributed to the probability mass function while the other 50 percent will be attributed to the right-tail cumulative distribution function.}

\item{weighted_rt}{Boolean used to determine if the weighted right-tail cumulative distribution function should be used or not.}

\item{left_trunc}{Int used to determine starting index of model to use for optimization}

\item{right_trunc}{Int used to determine ending index of model to use for optimization}
}
\description{
This function generates the Psi Criterion goodness of fit value given an empirical distribution. The function and parameters are given, as well as desired weight of pmf and use of the weighted right-tail cumulative distribution function.
}
\examples{
obs_data <- c(100,75,20,1)
parameters <- c(1,2,0.8)
psi <- psi_criterion_function(parameters, obs_data, Kolmogorov_Waring)
}
