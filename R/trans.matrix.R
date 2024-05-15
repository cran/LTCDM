#' @title Compute transition matrix
#'
#' @description
#' Function to compute transition matrix using classification results
#'
#' @param X a matrix containing the initial state (first column) and the transition state (second column).
#'
#' @return a 2 \eqn{\times} 2 matrix where rows represent initial states (0 and 1) and the columns represent transition states (0 and 1).
#'
#' @examples
#' initial_states <- c(1, 2, 1, 2)
#' final_states <- c(1, 1, 2, 2)
#' transition_matrix <- trans.matrix(data.frame(initial_states, final_states))
#' print(transition_matrix)
#' \dontrun{
#' # transition probabilities (corrected and updated)
#' t = 2 # the number of time points
#' K = ncol(Q) # the number of attributes
#' Z = dat1[, c(1,2)]
#' z_t1 = cbind(1, Z$gender)  # Covariate at time 1
#' z_t2 = cbind(1, Z$gender, Z$intervention, apply(Z,1,prod)) # Covariates at time 2
#' beta = step3.output$beta
#' gamma_01 = step3.output$gamma_01
#' gamma_10 = step3.output$gamma_10
#' updated.class <- update.class(cep = cep, K = K, t = t, z_t1 = z_t1,
#'                               z_t2 = z_t2, beta = beta, gamma_01 = gamma_01, gamma_10 = gamma_10)
#' C.eap.t1 = updated.class$cor.profile[[1]]
#' C.eap.t2 = updated.class$cor.profile[[2]]
#' C.eap.t1t2 <- cbind(z_t2, C.eap.t1, C.eap.t2)
#' t.A1.c = trans.matrix(as.matrix(C.eap.t1t2[,c("A1_t1","A1_t2")]))
#' t.A1.c
#' }
#' @export

#### Transition matrix (computed using classification results)
trans.matrix <- function(X)
{
  tt <- table( c(X[,-ncol(X)]), c(X[,-1]) )
  tt <- tt / rowSums(tt)
  return(tt)
}



