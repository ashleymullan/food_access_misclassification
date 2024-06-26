#' Observed-data log-likelihood for the sieve maximum likelihood estimator (SMLE)
#'
#' This function returns the value of the observed-data log-likelihood (equation (#) in Lotspeich et al. (2023+))
#' for a given dataset and parameter values `theta` and `p`.
#
#'
#' @param Y Column name with the outcome
#' @param offset (Optional) Column name with the offset for \code{Y}. Default is \code{offset = 1}, no offset
#' @param X_unval Column(s) with the unvalidated covariates 
#' @param X_val Column(s) with the validated covariates 
#' @param Z (Optional) Column(s) with additional error-free covariates 
#' @param Bspline Vector of columns containing the B-spline basis functions 
#' @param comp_dat_val Dataset containing rows for each validated subjects' data with values from Phase II (a matrix)
#' @param comp_dat_unval Augmented dataset containing rows for each combination of unvalidated subjects' data with values from Phase II (a matrix)
#' @param theta_pred Vector of columns in \code{data} that pertain to the covariates in the analysis model.
#' @param theta Parameters for the analysis model (a column vector)
#' @param p B-spline coefficients for the approximated covariate error model (a matrix)
#' @return Scalar value of the function

smle_observed_data_loglik = function(Y = NULL, offset = NULL, X_unval = NULL, X_val = NULL, Z = NULL, Bspline = NULL, comp_dat_val, comp_dat_unval, theta_pred, theta, p) {
  # Determine error setting -----------------------------------------
  ## If unvalidated variable was left blank, assume error-free ------
  errorsX = !is.null(X_unval)
  ## ------ If unvalidated variable was left blank, assume error-free
  # ----------------------------------------- Determine error setting

  if (errorsX) {
    m = nrow(p)
  }

  # For validated subjects --------------------------------------------------------
  #################################################################################
  ## Sum over log[P_theta(Yi|Xi)] -------------------------------------------------
  mu_theta = exp(as.numeric((cbind(int = 1, comp_dat_val[, theta_pred]) %*% theta)))
  pY_X = dbinom(x = comp_dat_val[, Y], 
                size = comp_dat_val[, offset], 
                prob = mu_theta)
  log_pY_X = log(pY_X)
  log_pY_X[log_pY_X == -Inf] = 0
  return_loglik = sum(log_pY_X)
  ## ------------------------------------------------- Sum over log[P_theta(Yi|Xi)]
  #################################################################################
  if (errorsX) {
    ## Sum over I(Xi=xk)Bj(Xi*)log p_kj ---------------------------------------------
    pX = p[comp_dat_val[, "k"], ]
    log_pX = log(pX)
    log_pX[log_pX == -Inf] = 0
    return_loglik = return_loglik + sum(comp_dat_val[, Bspline] * log_pX)
    ## --------------------------------------------- Sum over I(Xi=xk)Bj(Xi*)log q_kj
  }
  #################################################################################
  # -------------------------------------------------------- For validated subjects

  # For unvalidated subjects ------------------------------------------------------
  ## Calculate P_theta(y|x) for all (y,xk) ----------------------------------------
  mu_theta = exp(as.numeric((cbind(int = 1, comp_dat_unval[, theta_pred]) %*% theta)))
  pY_X = dbinom(x = comp_dat_unval[, Y], 
                size = comp_dat_unval[, offset], 
                prob = mu_theta)
  ## ---------------------------------------- Calculate P_theta(y|x) for all (y,xk)
  ################################################################################
  if (errorsX) {
    ## Calculate Bj(Xi*) p_kj for all (k,j) ----------------------------------------
    pX = p[comp_dat_unval[, "k"], ]
    ## ---------------------------------------- Calculate Bj(Xi*) p_kj for all (k,j)
  } else {
    pX = rep(1, nrow(comp_dat_unval))
  }
  ################################################################################
  ## Calculate sum of P(y|xk) x Bj(X*) x p_kj ------------------------------------
  if (errorsX) {
    person_sum = rowsum(x = pY_X * pX * comp_dat_unval[, Bspline], 
                        group = rep(seq(1, nrow(comp_dat_unval)), times = m), 
                        reorder = FALSE)
  }
  person_sum = rowSums(person_sum)
  log_person_sum = log(person_sum)
  log_person_sum[log_person_sum == -Inf] = 0
  ## And sum over them all -------------------------------------------------------
  return_loglik = return_loglik + sum(log_person_sum)
  ################################################################################
  # ----------------------------------------------------- For unvalidated subjects
  return(return_loglik)
}
