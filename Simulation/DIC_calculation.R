deviance_fun <- function(beta_mat, phi_vec, psi_vec, E_it, y_mat, x_cube){
  # dim: beta -> p * n, phi -> n, psi -> T, y -> 
  n <- length(phi_vec)
  T <- length(psi_vec)
  eta <- matrix(0, n, T)
  loglik <- 0
  for(i in 1:n){
    eta[i, ] = x_cube[, ,i]%*%beta_mat[,i,drop = F] + psi_vec + phi_vec[i]
    loglik <- loglik + sum(dpois(y_mat[i, ], E_it[, i]*exp(eta[i, ]), log = TRUE))
  }
  return(-2*loglik)
}



calculate_DIC <- function(beta_post, phi_post, psi_post, E_it, y_mat, x_cube){
  iter_num <- ncol(phi_post)
  beta_est <- apply(beta_post, c(1,2), mean)
  phi_est <- apply(phi_post, c(1), mean)
  psi_est <- apply(psi_post, c(1), mean)
  D_post <- 0
  for(iter in 1:iter_num){
    D_post <- D_post + deviance_fun(beta_post[, , iter], phi_post[, iter], psi_post[, iter], E_it, y_mat, x_cube)/iter_num
  }
  D_postMean <- deviance_fun(beta_est, phi_est, psi_est, E_it, y_mat, x_cube)
  pD <- D_post - D_postMean
  return(D_postMean + 2*pD)
}


