#' @title Compute classification error probabilities for attributes at different time points
#'
#' @description
#' Function to compute classification error probabilities (CEP) for attributes at different time points. Only attribute-level CEP is available for the time being.
#'
#' @references
#' \itemize{
#'   \item{Liang, Q., de la Torre, J., & Law, N. (2023).Latent transition cognitive diagnosis model with covariates: A three-step approach. \emph{Journal of Educational and Behavioral Statistics}. \doi{10.3102/10769986231163320}}
#'   \item{Huebner, A., & Wang, C. (2011).A note on comparing examinee classification methods for cognitive diagnosis models. \emph{Educational and Psychological Measurement, 71}, 407-419. \doi{10.1177/0013164410388832}}
#' }
#'
#' @param fit.object a list of the G-DINA model objects return from \code{GDINA} R package at pre-and post-tests.
#' @param K the number of attributes.
#' @param t the number of time points. This package can only handle two time points can for the time being.
#' @param N the number of examinees (observations).
#'
#' @return a list with elements
#' \describe{
#' \item{cep.matrix}{the CEP matrix}
#' \item{w}{the correction weights}
#' \item{mp}{the estimated marginal posterior probabilities obtained from \code{GDINA} R package}
#' \item{eap}{the estimated EAP of attribute profiles obtained from \code{GDINA} R package}
#' }
#' @import ggplot2 ggpubr ggsignif utils
#'
#' @examples
#' if(requireNamespace("GDINA")){
#' library(GDINA)
#' # Assuming dat0, dat1, Q, and other necessary data and objects are predefined.
#' rdmodel <- c("GDINA","GDINA","GDINA","GDINA","GDINA","GDINA","GDINA","GDINA",
#' "GDINA","GDINA","GDINA","GDINA","GDINA","GDINA","GDINA","GDINA",
#' "GDINA","GDINA","GDINA","GDINA","RRUM","GDINA","GDINA","GDINA",
#' "GDINA","LLM","LLM","RRUM","ACDM","GDINA","GDINA","GDINA",
#' "GDINA","GDINA","GDINA","GDINA","GDINA","GDINA","GDINA","GDINA")
#' fitrd <- GDINA(dat = dat0, Q = Q, model= rdmodel, mono.constraint = TRUE, verbose=0)
#'
#' # Obtained the item parameters from Tan et al. (2022)
#' itemparm.rd = GDINA::extract(fitrd,"catprob.parm")
#'
#' # Fit the response data at pre-test to the selected models
#' fit.t1 = GDINA(dat = dat1[,3:42], Q = Q, mono.constraint = TRUE, model = rdmodel,
#' catprob.parm = itemparm.rd, att.dist = "independent", control=list(maxitr = 0), verbose=0)
#'
#' # Fit the response data at post-test to the selected models
#' fit.t2 = GDINA(dat = dat1[,43:82], Q = Q, mono.constraint = TRUE, model = rdmodel,
#' catprob.parm = itemparm.rd, att.dist = "independent", control=list(maxitr = 0), verbose=0)
#' fit.object = list()
#' fit.object[[1]] <- fit.t1
#' fit.object[[2]] <- fit.t2
#' t = 2 # the number of time points
#' K = ncol(Q) # the number of attributes
#' N = nrow(dat1) # the number of observations
#'
#' cep = CEP_t(fit.object = fit.object, t = t, K = K, N = N)
#'
#' # The CEP matrices of the attributes
#' cep$cep.matrix
#'
#' }
#' @export

CEP_t <- function(fit.object, t, K, N){

  # fit.object: The G-DINA objects, a list containing all measurement models for all time points
  # t: The number of time points
  # K: The number of attributes
  # N: The sample size

  mp = list()
  eap = list()
  prob = list()
  w.total = list()
  cep.total = list()

  for(tt in 1:t){
    GDINA.obj = fit.object[[tt]]
    mp[[tt]] <- personparm(GDINA.obj, what = "mp") # Person Marginal posterior
    eap[[tt]] <- personparm(GDINA.obj, what = "EAP")
    prob[[tt]] = GDINA::extract(GDINA.obj, what = "prevalence")
  }

  #================== CEP ================
  for(k in 1:K){

    indicator = list()
    Den = list()
    P_ik = list()

    for(tt in 1:t){
      indicator[[tt]] = cbind(1-(eap[[tt]][,k]),eap[[tt]][,k]) # the indicator function/w_irt
      Den[[tt]] = colSums(cbind(1-(mp[[tt]][,k]),mp[[tt]][,k])) # the denominator
      P_ik[[tt]] = cbind(1-(mp[[tt]][,k]),mp[[tt]][,k])
    }

    num_00 = num_01= num_10 = num_11 = 0 # the numerators
    for(tt in 1:t){

      # true=0, observed=0
      num_00 = num_00 + P_ik[[tt]][,1]*indicator[[tt]][,1]
      # true=1, observed=0
      num_10 = num_10 + P_ik[[tt]][,2]*indicator[[tt]][,1]
      # true=0, observed=1
      num_01 = num_01 + P_ik[[tt]][,1]*indicator[[tt]][,2]
      # true=1, observed=1
      num_11 = num_11 + P_ik[[tt]][,2]*indicator[[tt]][,2]
    }

    cep.matrix = list()
    for(tt in 1:t){
      # tt=1
      cep.matrix[[tt]] = matrix(0,2,2)

      # true=0, observed=0
      cep.matrix[[tt]][1,1] = sum(num_00)/(N*t*prob[[tt]]$all[k,1])
      # true=1, observed=0
      cep.matrix[[tt]][2,1] = sum(num_10)/(N*t*prob[[tt]]$all[k,2])
      # true=0, observed=1
      cep.matrix[[tt]][1,2] = sum(num_01)/(N*t*prob[[tt]]$all[k,1])
      # true=1, observed=1
      cep.matrix[[tt]][2,2] = sum(num_11)/(N*t*prob[[tt]]$all[k,2])
      # cep.matrix

    }
    cep.matrix
    cep.total[[k]] = cep.matrix

    #============Correction weights================
    w.k = list()
    for(tt in 1:t){
      w.k[[tt]] = matrix(0,N,2)
      for(i in 1:N){
        if(eap[[tt]][i,k]){
          w.k[[tt]][i,] = t(cep.matrix[[tt]][,2])
        }else{w.k[[tt]][i,] = t(cep.matrix[[tt]][,1])}
      }
    }

    w.total[[k]] = w.k

  }

  return(list(cep.matrix = cep.total, w = w.total, mp = mp, eap = eap))

}
#.onAttach <- function(libname, pkgname) {
#  data(list = "Data_example", package = pkgname, envir = .GlobalEnv)
#  packageStartupMessage("Data_example is now loaded in the global environment.")
#}

