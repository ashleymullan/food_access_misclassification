# Load data from GitHub
food_access = read.csv(
  "https://raw.githubusercontent.com/sarahlotspeich/food/main/piedmont-triad-data/analysis_data.csv"
) |> 
  dplyr::select(-CountyName, -Y_BPHIGH, -Y_CHD, -Y_OBESITY) 

## Create binary indicators from X and Xstar
food_access = food_access |> 
  dplyr::mutate(
    Xstar_0.5 = as.numeric(Xstar <= 0.5), ### 0.5 miles
    Xstar_1 = as.numeric(Xstar <= 1), ### 1.0 miles
    Xstar_5 = as.numeric(Xstar <= 5), ### 5.0 miles
    Xstar_10 = as.numeric(Xstar <= 10) ### 10.0 miles
  ) |> 
  tidyr::gather(key = "THRESH", value = "binXstar", -c(1:6)) |> 
  dplyr::mutate(
    THRESH = as.numeric(sub(pattern = ".*_", 
                            replacement = "", 
                            x = THRESH)), 
    binX_full = as.numeric(X_full <= THRESH), 
    binX_partial = as.numeric(X_partial <= THRESH)
  ) |> 
  dplyr::filter(THRESH == 1)

logbinom_fit = glm(formula = cbind(Y_DIABETES, O_POP - Y_DIABETES) ~ binX_full, 
                   family = binomial(link = log),
                   data = food_access)
coefficients(summary(logbinom_fit))
cc_logbinom_fit = glm(formula = cbind(Y_DIABETES, O_POP - Y_DIABETES) ~ binX_partial, 
                      family = binomial(link = log),
                      data = food_access)
coefficients(summary(cc_logbinom_fit))
# Estimate  Std. Error    z value      Pr(>|z|)
# (Intercept) -2.2456410 0.002697088 -832.61696  0.000000e+00
# binX_full    0.1549768 0.004696674   32.99714 8.926034e-239
logistic_fit = glm(formula = cbind(Y_DIABETES, O_POP - Y_DIABETES) ~ binX_full, 
                   family = binomial(link = logit),
                   data = food_access)
coefficients(summary(logistic_fit))
# Estimate  Std. Error    z value     Pr(>|z|)
# (Intercept) -2.1337484 0.003016403 -707.38167  0.00000e+00
# binX_full    0.1750227 0.005324248   32.87275 5.38945e-237
poisson_fit = glm(formula = Y_DIABETES ~ binX_full, 
                  family = poisson(link = log),
                  offset = log(O_POP),
                  data = food_access)
coefficients(summary(poisson_fit))
# Estimate  Std. Error    z value      Pr(>|z|)
# (Intercept) -2.2456410 0.002852281 -787.31418  0.000000e+00
# binX_full    0.1549768 0.005000516   30.99217 6.873773e-211

glm_fit <- glm(cbind(Y_DIABETES, O_POP - Y_DIABETES) ~ binX_full, 
               family = binomial(link = log),
               data = food_access)
## (Intercept)   binX_full 
## -2.2456410   0.1549768 
library(logbin)
logbin_fit <- logbin(formula = cbind(Y_DIABETES, O_POP - Y_DIABETES) ~ binX_full, 
                     data = food_access)# , start = rep(-0.2, 9))

neg_log_lik = function(beta_ppv, Y, Offset, Xstar = NULL, X) {
  ## Define query indicator
  Q = as.numeric(!is.na(X))
  
  ## Separate parameters
  beta0 = beta_ppv[1]
  beta1 = beta_ppv[2]
  if (length(beta_ppv) > 2) {
    phi = beta_ppv[3] ### logit of the PPV (since it's unconstrained)
    ppv = exp(phi) / (1 + exp(phi))
  }
  
  ## Compute log-likelihood contribution for queried neighborhoods
  ### First, log{P(Y|X)}
  mu = exp(beta0 + beta1 * X)
  # if (any(mu > 1)) {
  #   return(9999)
  # }
  log_pYgivX = Y * log(mu) + (Offset - Y) * log(1 - mu)
  #log_pYgivX[which(log_pYgivX == -Inf)] = 0 #### Replace -Inf with 0
  ll = sum(log_pYgivX, na.rm = TRUE)
  
  if (mean(Q) < 1) {
    ### Then, log{P(X|X*)}
    log_pXgivXstar = X * log(ppv) + (1 - X) * log(1 - ppv)
    #log_pXgivXstar[which(log_pXgivXstar == -Inf)] = 0 #### Replace -Inf with 0
    ll = ll + sum(Q * log_pXgivXstar, na.rm = TRUE)
    
    ## Compute log-likelihood contribution for unqueried neighborhoods
    ### Loop over unqueried neighborhoods and x = 0, 1
    for (i in which(Q == 0)) {
      l_i = 0
      for (x in 0:1) {
        #### First, P(Y|X) 
        mu_x = exp(beta0 + beta1 * x)
        pYgivX = dbinom(x = Y[i], size = Offset[i], prob = mu_x)
        #### Then, P(X|X*)
        pXgivXstar = (ppv ^ x) * ((1 - ppv) ^ (1 - x))
        #### Save P(Y|X)P(X|X*)
        l_i = l_i + pYgivX * pXgivXstar
      }
      ll_i = log(l_i) 
      ll_i = ifelse(test = ll_i == -Inf, yes = 0, no = ll_i)
      ll = ll + ll_i
    }
    
  }
  
  ## Return log-likelihood (negated)
  return(-ll)
}

nlm_fit <- nlm(f = neg_log_lik, 
               p = c(log(0.5), 0, 1),
               Y = food_access$Y_DIABETES, 
               Offset = food_access$O_POP, 
               X = food_access$binX_partial, 
               hessian = TRUE)

log_likelihood <- function(theta, Y, Offset, X)
{
  p = exp(theta[1] + theta[2] * X)
  ll = sum(Y * log(p) + (Offset - Y) * log(1 - p))
  return(ll)
}

neg_log_likelihood <- function(theta, Y, Offset, X)
{
  ll = log_likelihood(theta = theta, Y = Y, Offset = Offset, X = X)
  return(-ll)
}
nlm_fit2 <- nlm(f = neg_log_likelihood, 
               p = c(log(0.5), 0),
               Y = food_access$Y_DIABETES, 
               Offset = food_access$O_POP, 
               X = food_access$binX_full)
