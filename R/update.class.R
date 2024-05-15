#' @title Classification update using the Bayes' Theorem
#'
#' @description
#' Function to update classifications (attribute profiles) using the Bayes' Theorem
#'
#' @references
#' Liang, Q., de la Torre, J., & Law, N. (2023). Latent transition cognitive diagnosis model with covariates: A three-step approach. \emph{Journal of Educational and Behavioral Statistics}.\doi{10.3102/10769986231163320}
#'
#' @param cep estimated classification error probabilities returned from \code{\link{CEP_t}}. The uncorrected attribute profile (EAP) is also stored in this object.
#' @param K the number of attributes.
#' @param t the number of time points. This package can only handle two time points can for the time being.
#' @param z_t1 covariates at Time 1, which has already had the intercept column (1s).
#' @param z_t2 covariates at Time 2, which has already had the intercept column (1s).
#' @param beta the estimated regression coefficients at Time 1 (initial state)
#' @param gamma_01 the estimated regression coefficients of transition from absence (0) to presence (1) at Time 2
#' @param gamma_10 the estimated regression coefficients of transition from absence (0) to presence (1) at Time 2
#'
#' @return a list with elements
#' \describe{
#' \item{post.prob}{the corrected and updated posterior probability}
#' \item{att.prevalance}{the corrected and updated attribute prevalance}
#' \item{cor.profile}{the corrected and updated attribute profiles for different time points}
#' }
#'
#' @examples
#' \dontrun{
#' #The run is dependent on the output of the step3.est() function and CEP_t() function
#' #It is not recommended for run it.
#' t = 2 # the number of time points
#' K = ncol(Q) # the number of attributes
#' Z = dat1[, c(1,2)]
#' z_t1 = cbind(1, Z$gender)  # Covariate at time 1
#' z_t2 = cbind(1, Z$gender, Z$intervention, apply(Z,1,prod)) # Covariates at time 2
#' beta = step3.output$beta
#' gamma_01 = step3.output$gamma_01
#' gamma_10 = step3.output$gamma_10
#'
#' # Update classifications using the Bayes' Theorem
#' updated.class <- update.class(cep = cep, K = K, t = t, z_t1 = z_t1, z_t2 = z_t2,
#'                               beta = beta, gamma_01 = gamma_01, gamma_10 = gamma_10)
#'
#' # The corrected and updated attribute prevalance
#' updated.class$att.prevalance
#'
#' # The corrected and updated posterior probability
#' updated.class$post.prob
#' }
#' @export

update_class <- function(cep, K, t, z_t1, z_t2, beta, gamma_01, gamma_10){

  N <- nrow(z_t1)  # 或者是 nrow(z_t2)，取决于N代表的具体含义

  if(t != 2)
    stop("This package can only handle two time points can for the time being.", call. = FALSE)

  #### Get the posterior probabilities and attribute profiles (EAP) estimated in Step 1
  # Posterior probabilities
  mp = NULL
  for(tt in 1:t){
    mp = cbind(mp, cep$mp[[tt]])
  }
  colnames(mp) <- c(paste0("A",seq_len(K),"_t1"),
                    paste0("A",seq_len(K),"_t2"))

  # Attribute profiles (EAP)
  eap.dat = NULL
  for (tt in 1:t) {
    eap.dat = cbind(eap.dat, cep$eap[[tt]])
  }
  eap.dat = cbind(eap.dat, z_t1, z_t2)
  colnames(eap.dat)[1:8] <- c(paste0("A",seq_len(K),"_t1"),
                              paste0("A",seq_len(K),"_t2"))

  ### corrected attribute profiles (EAP)
  cor.profile = list()

  #================================initial state=================================#
  # time 1
  mp_t1 = mp[,c(paste0("A",seq_len(K),"_t1"))]
  mp.cep.t1 = exp(z_t1%*%beta) / (1 + exp(z_t1%*%beta))

  # compute updated EAP based on Bayes theorem
  mp.cep.comb.t1 = (mp_t1[,1:K]*mp.cep.t1)/(mp_t1[,1:K]*mp.cep.t1 + (1-mp_t1[,1:K])*(1-mp.cep.t1))

  #========= Obtain EAP estimates based on updated marginal posterior =========
  # time 1
  C.eap.t1 <- 1*(mp.cep.comb.t1 > 0.5)
  cor.profile[[1]] = C.eap.t1

  #=============================Transition probability===========================#
  # time 2
  mp_t2 = mp[,c(paste0("A",seq_len(K),"_t2"))]
  tran01.cep.t2 = exp(z_t2%*%gamma_01) / (1 + exp(z_t2%*%gamma_01)) # the transition probability from 0 to 1
  tran10.cep.t2 = exp(z_t2%*%gamma_10) / (1 + exp(z_t2%*%gamma_10)) # the transition probability from 1 to 1

  mp.cep.t2 = NULL

  for (k in 1:K) {

    mp.cep.comb.t1.k = cbind(1-mp.cep.comb.t1[,k], mp.cep.comb.t1[,k])
    mp.cep.t2.k = matrix(0, N, 2) #N的具体含义

    for(i in 1:N){
      # update the transition matrix
      tran.matrix.cep = matrix(c(1-tran01.cep.t2[i,k], tran01.cep.t2[i,k],
                                 tran10.cep.t2[i,k], 1-tran10.cep.t2[i,k]), 2, 2, T)

      mp.cep.t2.k[i,] = t(tran.matrix.cep) %*% mp.cep.comb.t1.k[i,]
    }

    mp.cep.t2 = cbind(mp.cep.t2, mp.cep.t2.k[,2])

  }

  # compute updated EAP based on Bayes theorem
  mp.cep.comb.t2 = (mp.cep.t2*mp_t2)/(mp.cep.t2*mp_t2 + (1-mp.cep.t2)*(1-mp_t2))
  colMeans(mp.cep.comb.t2)

  C.eap.t2 <- 1*(mp.cep.comb.t2 > 0.5)
  cor.profile[[2]] = C.eap.t2

  ###########
  ### posterior probability (corrected and updated)
  post.prob = rbind(colMeans(mp.cep.comb.t1), colMeans(mp.cep.comb.t2))

  ###########
  ### attribute prevalence (corrected and updated)
  att.prevalance = rbind(colMeans(C.eap.t1), colMeans(C.eap.t2))


  return(list(post.prob = post.prob, att.prevalance = att.prevalance, cor.profile = cor.profile))
}



