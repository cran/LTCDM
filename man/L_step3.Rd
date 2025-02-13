% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/step3_est.R
\name{L_step3}
\alias{L_step3}
\title{Step 3 estimation for latent logistic regression coefficients}
\usage{
L_step3(par, z_t1, z_t2, weight, k)
}
\arguments{
\item{par}{Coefficients of latent logistic regression to be estimated.}

\item{z_t1}{covariates at Time 1, which has already had the intercept column (1s).}

\item{z_t2}{covariates at Time 2, which has already had the intercept column (1s).}

\item{weight}{Correction weight obtained from CEP.}

\item{k}{The k-th attribute.}
}
\value{
log-likelihood value.
}
\description{
Function to estimate the latent logistic regression models at the initial state and transition
}
\references{
Liang, Q., de la Torre, J., & Law, N. (2023). Latent transition cognitive diagnosis model with covariates: A three-step approach. \emph{Journal of Educational and Behavioral Statistics}.\doi{10.3102/10769986231163320}
}
