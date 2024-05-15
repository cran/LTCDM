#' @title Step 3 estimation for latent logistic regression coefficients
#'
#' @description
#' Function to estimate the latent logistic regression models at the initial state and transition
#'
#' @references
#' Liang, Q., de la Torre, J., & Law, N. (2023). Latent transition cognitive diagnosis model with covariates: A three-step approach. \emph{Journal of Educational and Behavioral Statistics}.\doi{10.3102/10769986231163320}
#' @importFrom stats optim
#' @importFrom stats pnorm
#'
#' @param z_t1 covariates at Time 1, which has already had the intercept column (1s).
#' @param z_t2 covariates at Time 2, which has already had the intercept column (1s).
#' @param par Coefficients of latent logistic regression to be estimated.
#' @param weight Correction weight obtained from CEP.
#' @param k The k-th attribute.
#'
#' @return log-likelihood value.
#'@import GDINA
#' @export


###############################################################
#                    Step 3 functions
#           The latent logistic regression model
#                and the objective function
#                   for two time points
###############################################################

### The latent logistic regression model and the objection function
### For attribute-level associations
L_step3 <- function(par, z_t1, z_t2, weight, k){

  N = nrow(z_t1) # sample size
  C1 = ncol(z_t1) #number of covariates at time 1; already has the intercept column
  C2 = ncol(z_t2) #number of covariates at time 2; already has the intercept column

  beta = as.matrix(par[1:C1])
  ga01 = as.matrix(par[(C1+1) : (C1+C2)])
  ga10 = as.matrix(par[(C1+C2+1) : (C1+C2+C2)])

  # initial state
  mp1 = exp(z_t1 %*% beta)/(1+exp(z_t1 %*% beta))

  # transition probability matrix
  tran_pro = data.frame(p00 = rep(0,N), p01 = rep(0,N), p10 = rep(0,N), p11 = rep(0,N))

  tran_pro$p01 = exp(z_t2 %*% ga01)/(1+exp(z_t2 %*% ga01))
  tran_pro$p00 = 1-tran_pro$p01

  tran_pro$p10 = exp(z_t2 %*% ga10)/(1+exp(z_t2 %*% ga10))
  tran_pro$p11 = 1-tran_pro$p10

  mp_t1.k = cbind(1-mp1, mp1)
  mp_t2.k = matrix(0, N, 2)
  for(i in 1:N){
    # i=1
    tran.matrix = matrix(c(tran_pro[i,1],tran_pro[i,2], tran_pro[i,3],tran_pro[i,4]), 2, 2, T)
    mp_t2.k[i,] = t(tran.matrix) %*% mp_t1.k[i,]
  }

  # The objective function
  ll = sum( log( (weight$w10*mp_t1.k[,1] + weight$w11*mp_t1.k[,2]) * (weight$w20*mp_t2.k[,1] + weight$w21*mp_t2.k[,2]) ))

  ll = -1*ll # optimization will minimize function

  return(ll)
}

###  The optimization function
out.ltm <- function(cep, z_t1, z_t2, K, t, beta_in, ga01_in, ga10_in){

  C1 = ncol(z_t1) # number of covariates at pre-test (has the intercept column)
  C2 = ncol(z_t2) # number of covariates at post-test (has the intercept column)
  N = nrow(z_t1)
  eap = cep$eap

  param = list(beta = beta_in, ga01 = ga01_in, ga10 = ga10_in) # regression parameters to be estimated

  res = NULL # estimated parameters for all attributes using correction weights

  for (k in 1:K) {

    w = NULL
    for (tt in 1:t) {
      w = cbind(w, cep$w[[k]][[tt]])
    }
    w = as.data.frame(w)
    colnames(w) <- c("w10", "w11","w20","w21")

    par = as.matrix(c(param$beta[,k], param$ga01[,k], param$ga10[,k]))

    #================= Optimize the modified objective function =================
    out = optim(par = par, fn = L_step3, weight = w, z_t1 = z_t1, z_t2 = z_t2, k = k,
                method="L-BFGS-B", lower= -3, upper = 3, hessian = T)

    res = cbind(res, out$par)

  }

  return(list(res = res))
}

#### Compute the standard errors (SE) of the regression coefficients ####
SE <- function(beta, gamma_01, gamma_10, z_t1, z_t2, N, K){

  SE_in = NULL # SEs for initial state coefficients
  SE_t01 = NULL # SEs for coefficients of transition from 0 to 1
  SE_t10 = NULL # SEs for coefficients of transition from 1 to 0

  design_X1 = as.matrix(z_t1)
  design_X2 = as.matrix(z_t2)

  for(k in 1:K){

    #########Diagonal matrix##########
    # initial state regression
    p = exp(z_t1%*%beta[,k])/(1+exp(z_t1%*%beta[,k]))
    V_in = matrix(0, N, N)
    diag(V_in) <- p*(1-p)

    # transition probability regression: 0->1
    p_t01 = exp(z_t2%*%gamma_01[,k])/(1+exp(z_t2%*%gamma_01[,k]))
    V_t01 = matrix(0, N, N)
    diag(V_t01) <- p_t01*(1-p_t01)

    # transition probability regression: 1->0
    p_t10 = exp(z_t2%*%gamma_10[,k])/(1+exp(z_t2%*%gamma_10[,k]))
    V_t10 = matrix(0, N, N)
    diag(V_t10) <- p_t10*(1-p_t10)

    ##########compute SE###########
    # initial state regression
    cov_in = solve(t(design_X1)%*%V_in%*%design_X1)
    SE_in <- cbind(SE_in, sqrt(diag(cov_in)))

    # transition probability regression: 0->1
    cov_t01 = solve(t(design_X2)%*%V_t01%*%design_X2)
    SE_t01 <- cbind(SE_t01, sqrt(diag(cov_t01)))

    # transition probability regression: 1->1
    cov_t10 = solve(t(design_X2)%*%V_t10%*%design_X2)
    SE_t10 <- cbind(SE_t10, sqrt(diag(cov_t10)))

  }

  return(SE = list(SE_in = SE_in, SE_t01 = SE_t01, SE_t10 = SE_t10))

}


#### Compute the p-values and 95% CI of the regression coefficients ####
sig <- function(beta, gamma_01, gamma_10, SE){

  SE_in <- SE$SE_in
  SE_t01 <- SE$SE_t01
  SE_t10 <- SE$SE_t10
  K = ncol(gamma_01)

  wald_in = beta/SE_in
  sig_in.one = pnorm(-abs(wald_in)) #one-tail test
  sig_in = 2*pnorm(abs(wald_in),lower.tail = F) #two-tail test
  round(sig_in,4)
  CI_in_l = beta - 1.96*SE_in
  CI_in_r = beta + 1.96*SE_in

  wald_t01 = gamma_01/SE_t01
  sig_t01.one = pnorm(-abs(wald_t01)) #one-tail test
  sig_t01 = 2*pnorm(abs(wald_t01),lower.tail = F) #two-tail test
  round(sig_t01,4)
  CI_t01_l = gamma_01 - 1.96*SE_t01
  CI_t01_r = gamma_01 + 1.96*SE_t01

  wald_t10 = gamma_10/SE_t10
  sig_t10.one = pnorm(-abs(wald_t10)) #one-tail test
  sig_t10 = 2*pnorm(abs(wald_t10),lower.tail = F) #two-tail test
  round(sig_t10,4)
  CI_t10_l = gamma_10 - 1.96*SE_t10
  CI_t10_r = gamma_10 + 1.96*SE_t10

  res = NULL
  for(k in 1:K){

    npar = length(beta[,k])+length(gamma_01[,k])+length(gamma_10[,k])
    beta.name = paste0("beta_",c(0,seq_len(length(beta[,k])-1)))
    gamma01.name = paste0("gamma01_",c(0,seq_len(length(gamma_01[,k])-1)))
    gamma10.name = paste0("gamma10_",c(0,seq_len(length(gamma_10[,k])-1)))

    CI_l = round(c(CI_in_l[,k], CI_t01_l[,k], CI_t10_l[,k]),4)
    CI_r = round(c(CI_in_r[,k], CI_t01_r[,k], CI_t01_r[,k]),4)

    coef.est = c(beta[,k],gamma_01[,k],gamma_10[,k])

    res1 <- data.frame(Attribute = rep(k, length(coef.est)), Estimate = round(coef.est,4),
                       odds.ratio = round(exp(coef.est),4),
                       d = round(coef.est * sqrt(3)/pi, 4),
                       SE = round(c(SE_in[,k],SE_t01[,k],SE_t10[,k]),4),
                       CI.95 = paste0("(", CI_l, ",", CI_r, ")"), z.value = round(c(wald_in[,k],wald_t01[,k],wald_t10[,k]),4),
                       p.2tailed = round(c(sig_in[,k],sig_t01[,k],sig_t10[,k]),4),
                       p.1tailed = round(c(sig_in.one[,k],sig_t01.one[,k], sig_t10.one[,k]),4),
                       row.names = c(beta.name,gamma01.name,gamma10.name))
    res = rbind(res, res1)

  }

  return(res)

}

#' Data Set Q
#'
#' The Q-matrix empirically validated by Tan et al.(2023).
#'
#' @format A data frame with 40 rows and 4 columns.
"Q"

#' Data Set cep
#'
#' The classification error probabilities (CEP) can be obtained in this data example.
#'
#' @format A list containing:
#' \describe{
#'   \item{\code{cep.matrix}}{Each of the 4 lists includes two 2 by 2 matrices.}
#'   \item{\code{w}}{Each of the 4 lists includes two 2005 by 2 matrices.}
#'   \item{\code{mp}}{Each of the 4 lists includes two 2005 by 4 matrices.}
#'   \item{\code{eap}}{Each of the 4 lists includes two 2005 by 4 matrices.}
#' }
"cep"

#' Data Set dat0
#'
#' The dataset of time point 1 used in (Tan et al., 2023).
#'
#' @format A data frame with 719 rows and 40 columns.
"dat0"

#' Data Set dat1
#'
#' The longitudinal dataset used in this example. Items with a prefix "a" are for the pre-test, and items with a prefix "b" are for the post-test.
#'
#' @format A data frame with 2005 rows and 82 columns.
"dat1"

#' Data Set step3.output
#'
#' The output of step 3 estimation can be obtained in this data example.
#'
#' @format A list containing:
#' \describe{
#'   \item{\code{beta}}{A data frame with 2 rows and 4 columns.}
#'   \item{\code{gamma_01}}{A data frame with 4 rows and 4 columns.}
#'   \item{\code{gamma_10}}{A data frame with 4 rows and 4 columns.}
#'   \item{\code{result}}{A list containing the results of the estimation, with dimensions 40 by 9.}
#' }
"step3.output"

#' @title Step 3 estimation for latent logistic regression coefficients
#'
#' @param cep estimated classification error probabilities returned from \code{\link{CEP_t}}. The uncorrected attribute profile (EAP) is also stored in this object.
#' @param z_t1 covariates at Time 1, which has already had the intercept column (1s).
#' @param z_t2 covariates at Time 2, which has already had the intercept column (1s).
#' @param K the number of attributes.
#' @param t the number of time points. This package can only handle two time points can for the time being.
#' @param beta_in the initial values for the regression coefficients at Time 1 (initial state). Default are 0s.
#' @param ga01_in the initial values for the regression coefficients of transition from absence (0) to presence (1) at Time 2. Default are 0s.
#' @param ga10_in the initial values for the regression coefficients of transition from presence (1) to absence (0) at Time 2. Default are 0s.
#' @return a list with elements
#' \describe{
#' \item{beta}{A data frame with 2 rows and 4 columns, representing the estimated regression coefficients at Time 1 (initial state)}
#' \item{gamma_01}{A data frame with 4 rows and 4 columns, representing the estimated regression coefficients of transition from absence (0) to presence (1) at Time 2}
#' \item{gamma_10}{A data frame with 4 rows and 4 columns, representing the estimated regression coefficients of transition from absence (0) to presence (1) at Time 2}
#' \item{result}{A data frame with dimensions 40 by 9, containing the results of the estimation, including all regression coefficients and the corresponding odds ratios, Cohen's d, standard errors (SE), 95% confidence intervals, and p-values.}
#' }
#' @examples
#' t = 2 # the number of time points
#' K = ncol(Q) # the number of attributes
#' z_t1_test = matrix(sample(c(0, 1), size = 20, replace = TRUE), nrow = 10)
#' z_t2_test = matrix(sample(c(0, 1), size = 40, replace = TRUE), nrow = 10)
#' # Set appropriate initial values of the coefficients
#' # Initial values of initial state's regression coefficients
#' beta_in = matrix(0, ncol(z_t1_test), K)
#'
#' # Initial values of transition probability's regression coefficients
#' # These were computed using the raw data.
#' # When Gender coding is 1 = male, 0 = female:
#' ga01_in = cbind(c(-2.15, 0.56, 0.09, -0.79),
#'                 c(-1.6, 0.05, -0.01, -0.38),
#'                 c(-1.25, 0.06, -0.25, 0.14),
#'                 c(-1.18, -0.26, 0.04, 0.37))
#'                 #initial values of regression coefficients (for transition from 0 to 1)
#' ga10_in = cbind(c(-0.84, -0.18, -0.14, 0.23),
#'                 c(-0.18, 0.49, 0.44, -0.35),
#'                 c(-0.22, 0.18, 0.37, -0.45),
#'                 c(-0.49, 0.10, 0.43, 0.20))
#' cep_test = list()
#' cep_test[["mp"]][[1]] = matrix(runif(40,min = 0,max=1),nrow = 10)
#' cep_test[["mp"]][[2]] = matrix(runif(40,min = 0,max=1),nrow = 10)
#' cep_test[["eap"]][[1]] = matrix(runif(40,min = 0,max=1),nrow = 10)
#' cep_test[["eap"]][[2]] = matrix(runif(40,min = 0,max=1),nrow = 10)
#' for (i in 1:4){
#' cep_test[["cep_matrix"]][[i]]=list()
#' cep_test[["w"]][[i]]=list()
#' for (k in 1:2) {
#'   cep_test[["cep_matrix"]][[i]][[k]]=matrix(c(1,0.02,0.06,1),nrow = 2)
#'       cep_test[["w"]][[i]][[k]] = matrix(runif(20,min = 0,max=1),nrow = 10)
#'         }
#'        }
#' step3.output_test <- step3.est(cep = cep_test, z_t1 = z_t1_test, z_t2 = z_t2_test,
#'  K = K, t = t, beta_in, ga01_in, ga10_in)
#' \dontrun{
#' # The run is dependent on the output of the CEP_t() function
#' # And the process time takes more than 5s.
#' # It is not recommended to run it.
#' # Covariates
#' Z = dat1[, c(1,2)] # use intervention and gender as covariates
#' z_t1 = cbind(1, Z$gender)  # Covariate at time 1
#' z_t2 = cbind(1, Z$gender, Z$intervention, apply(Z,1,prod)) # Covariates at time 2
#' colnames(z_t1) <- c("intercept", "gender")
#' colnames(z_t2) <- c("intercept", "gender", "intervention", "intervention_gender")
#'
#' t = 2 # the number of time points
#' K = ncol(Q) # the number of attributes
#'
#' # Set appropriate initial values of the coefficients
#' # Initial values of initial state's regression coefficients
#' beta_in = matrix(0, ncol(z_t1), K)
#'
#' # Initial values of transition probability's regression coefficients
#' # These were computed using the raw data.
#' # When Gender coding is 1 = male, 0 = female:
#' ga01_in = cbind(c(-2.15, 0.56, 0.09, -0.79),
#'                 c(-1.6, 0.05, -0.01, -0.38),
#'                 c(-1.25, 0.06, -0.25, 0.14),
#'                 c(-1.18, -0.26, 0.04, 0.37))
#'                 #initial values of regression coefficients (for transition from 0 to 1)
#' ga10_in = cbind(c(-0.84, -0.18, -0.14, 0.23),
#'                 c(-0.18, 0.49, 0.44, -0.35),
#'                 c(-0.22, 0.18, 0.37, -0.45),
#'                 c(-0.49, 0.10, 0.43, 0.20))
#'                 #initial values of regression coefficients (for transition from 1 to 0)
#' # Step 3 estimation (This will take a few minutes)
#' step3.output <- step3.est(cep = cep, z_t1 = z_t1, z_t2 = z_t2, K = K,
#'                           t = t, beta_in = beta_in, ga01_in = ga01_in, ga10_in = ga10_in)
#'
#' # Obtain estimation results
#' step3.output$result
#'
#' # Latent logistic regression coefficients
#' beta = step3.output$beta
#' gamma_01 = step3.output$gamma_01
#' gamma_10 = step3.output$gamma_10
#' }
#' @export
step3.est <- function(cep, z_t1, z_t2, K, t, beta_in = matrix(0, ncol(z_t1), K), ga01_in = matrix(0, ncol(z_t2), K), ga10_in = matrix(0, ncol(z_t2))){

  N = nrow(z_t1) # Sample size

  if(t != 2)
    stop("This package can only handle two time points can for the time being.", call. = FALSE)

  #### Estimation ####
  res.obj <- out.ltm(cep = cep, z_t1 = z_t1, z_t2 = z_t2, K = K, t = t, beta_in = beta_in, ga01_in = ga01_in, ga10_in = ga10_in)

  #### Regression coefficients ####
  C1 = ncol(z_t1) #number of covariates at time 1; already has the intercept column
  C2 = ncol(z_t2) #number of covariates at time 2; already has the intercept column

  beta = res.obj$res[1:C1,]
  gamma_01 = res.obj$res[(C1+1) : (C1+C2),]
  gamma_10 = res.obj$res[(C1+C2+1) : (C1+C2+C2),]
  se <- SE(beta = beta, gamma_01 = gamma_01, gamma_10 = gamma_10, z_t1 = z_t1, z_t2 = z_t2, N = N, K = K)
  res <- sig(beta = beta, gamma_01 = gamma_01, gamma_10 = gamma_10, SE = se)

  return(list(beta = beta, gamma_01 = gamma_01, gamma_10 = gamma_10, result = res))
}
