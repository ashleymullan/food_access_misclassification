#' @title
#' Observed-data log-likelihood
#'
#' @description
#' This function returns the value of the observed-data log-likelihood (equation (2) in Lotspeich et al. (2020))
#' for a given dataset and parameter values `theta`, `gamma`, and `p`.
#
#'
#' @param N Phase I sample size
#' @param n Phase II sample size
#' @param Y Column with the number of events (can be name or numeric index)
#' @param Size Column with the number of trials (can be name or numeric index)
#' @param X_unval Column(s) with the unvalidated predictors (can be name or numeric index)
#' @param X Column(s) with the validated predictors (can be name or numeric index)
#' @param Z (Optional) Column(s) with additional error-free covariates (can be name or numeric index)
#' @param Bspline Vector of columns containing the B-spline basis functions (can be name or numeric index)
#' @param comp_dat_all Augmented dataset containing rows for each combination of unvalidated subjects' data with values from Phase II (a matrix)
#' @param theta_pred Predictors in the analysis model
#' @param theta Parameters for the analysis model (a column vector)
#' @param p B-spline coefficients for the approximated covariate error model (a matrix)
#' @return Scalar value of the function
#' @noRd

logistic2ph_summ_loglik <- function(N, n, Y = NULL, Size = NULL, X_unval = NULL, X = NULL, Z = NULL, Bspline = NULL, comp_dat_all, theta_pred, theta, p) {
  #sn <- ncol(p)
  m <- nrow(p)
  comp_dat_all <- data.matrix(comp_dat_all)

  # For validated subjects --------------------------------------------------------
  #################################################################################
  ## Sum over log[P_theta(Yi|Xi)] -------------------------------------------------
  pY_X <- 1 / (1 + exp(-as.numeric((cbind(int = 1, comp_dat_all[c(1:n), theta_pred]) %*% theta))))
  pY_X <- ifelse(as.vector(comp_dat_all[c(1:n), c(Y)]) == 0, 1 - pY_X, pY_X)
  return_loglik <- sum(log(pY_X))
  ## ------------------------------------------------- Sum over log[P_theta(Yi|Xi)]
  #################################################################################
  ## Sum over I(Xi=xk)Bj(Xi*)log p_kj ---------------------------------------------
  pX <- p[comp_dat_all[c(1:n), "k"], ]
  log_pX <- log(pX)
  log_pX[log_pX == -Inf] <- 0
  return_loglik <- return_loglik + sum(comp_dat_all[c(1:n), Bspline] * log_pX)
  ## --------------------------------------------- Sum over I(Xi=xk)Bj(Xi*)log q_kj
  #################################################################################
  # -------------------------------------------------------- For validated subjects

  # For unvalidated subjects ------------------------------------------------------
  ## Calculate P_theta(y|x) for all (y,xk) ----------------------------------------
  pY_X <- 1 / (1 + exp(-as.numeric((cbind(int = 1, comp_dat_all[-c(1:n), theta_pred]) %*% theta))))
  pY_X[which(comp_dat_all[-c(1:n), Y] == 0)] <- 1 - pY_X[which(comp_dat_all[-c(1:n), Y] == 0)]
  ## ---------------------------------------- Calculate P_theta(y|x) for all (y,xk)
  ################################################################################
  ## Calculate Bj(Xi*) p_kj for all (k,j) ----------------------------------------
  pX <- p[comp_dat_all[-c(1:n), "k"], ]
  ## ---------------------------------------- Calculate Bj(Xi*) p_kj for all (k,j)
  ################################################################################
  ## Calculate sum of P(y|xk) x P(Y*|X*,y,xk) x Bj(X*) x p_kj --------------------
  person_sum <- rowsum(c(pY_X * pX) * comp_dat_all[-c(1:n), Bspline],
                       group = rep(seq(1, (N - n)),
                                   times = m))
  person_sum <- rowSums(person_sum)
  log_person_sum <- log(person_sum)
  log_person_sum[log_person_sum == -Inf] <- 0
  ## And sum over them all -------------------------------------------------------
  return_loglik <- return_loglik + sum(log_person_sum)
  ################################################################################
  # ----------------------------------------------------- For unvalidated subjects
  return(return_loglik)
}
