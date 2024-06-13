# Install libraries
## Run once
## if needed: install.packages("devtools") 
## then: 
## devtools::install_github("sarahlotspeich/possum", 
##                          ref = "main") 

# Load libraries
library(possum) ## for SMLE
library(tictoc) ## to calculate runtime

# Set working directory 
setwd("~/Documents/food_access_misclassification/Simulations/")

# Random seed to be used for each simulation setting
sim_seed = 11422

# Number of replicates per simulation setting
num_reps = 1000
## Note: You may want to run this code was run in parallel a cluster instead of locally, as it can be slow.

# Set parameters that won't be varied in the loop
## These values will be set as the defaults in the sim_data() function for convenience
fix_beta0 = -2.2456410 ## outcome model intercept (leads to ~ 11% prevalence)
fix_beta1 = 0.1549768 ## log prevalence ratio for X on Y
fix_q = 0.1 ## proportion of neighborhoods queried 

# ---------------------------------------------------------------------------------
# Function to simulate data (arguments defined as follows) ------------------------
## N = number of neighborhoods (sample size)
## beta0 = model intercept
## beta1 = log prevalence ratio 
# ---------------------------------------------------------------------------------
sim_data = function(N, beta0 = fix_beta0, beta1 = fix_beta1, beta2 = fix_beta2, ppv) {
  ## Simulate straight-line access
  binXstar = rbinom(n = N, 
                    size = 1, 
                    prob = 0.5)
  
  ## Simulate map-based access
  binX = binXstar
  binX[binXstar == 1] = rbinom(n = sum(binXstar == 1), 
                               size = 1, 
                               prob = ppv)
  
  ## Simulate population
  P = rpois(n = N, 
            lambda = 4165)
  
  ## Simulate cases of health outcome
  lambda = exp(beta0 + beta1 * binX)
  Cases = rpois(n = N, 
                lambda = P * lambda)
  
  ## Create dataset
  dat = data.frame(id = 1:N, binX, binXstar, P, Cases)
  
  # Return dataset
  return(dat)
}

# Loop over different sample sizes: N = 390 (Piedmont Triad), 2200 (all of NC)
for (N in c(390, 2200)) {
  tic(paste("Sims with N =", N)) ## Start counting runtime for sims with current sample size N
  # And proportion to be queried for complete case/imputation analyses: 0.1, 0.25, 0.5, 0.75
  for (ppv in c(0.5, 0.6, 0.8, 0.9)){
    # Be reproducible
    set.seed(sim_seed) ## set random seed
    
    # Create dataframe to save results for setting
    sett_res = data.frame(
      sim = paste(sim_seed, 1:num_reps, sep = "-"), 
      N, beta0 = fix_beta0, beta1 = fix_beta1, q = fix_q,  ## simulation setting
      ppv = ppv, avg_prev = NA, fpr = NA, ## simulation setting
      beta0_gs = NA, se_beta0_gs = NA, beta1_gs = NA, se_beta1_gs = NA, eta0_gs = NA, se_eta0_gs = NA, ## gold standard analysis
      beta0_n = NA, se_beta0_n = NA, beta1_n = NA, se_beta1_n = NA, ## naive analysis
      beta0_cc = NA, se_beta0_cc = NA, beta1_cc = NA, se_beta1_cc = NA, eta0_cc = NA, se_eta0_cc = NA, ## complete case analysis
      beta0_mle = NA, se_beta0_mle = NA, beta1_mle = NA, se_beta1_mle = NA, eta0_mle = NA, se_eta0_mle = NA  ## MLE analysis
    )
    
    # Loop over replicates 
    for (r in 1:num_reps) {
      # Generate data
      dat = sim_data(N  = N, ## sample size
                     ppv = ppv) ## positive predictive value
      
      # Save average neighborhood prevalence
      sett_res$avg_prev[r] = mean(dat$Cases / dat$P)
      
      # Save false positive rate
      sett_res$fpr[r] = sum(dat$binX == 0 & dat$binXstar == 1) /
        (sum(dat$binX == 0 & dat$binXstar == 1) + 
           sum(dat$binX == 0 & dat$binXstar == 0))
      
      # Fit the gold standard models
      ## Analysis model
      fit_gs = glm(formula = Cases ~ binX, 
                   family = poisson,
                   offset = log(P),
                   data = dat)
      sett_res[r, c("beta0_gs", "beta1_gs")] = coefficients(fit_gs) ## estimated log prevalence ratio
      sett_res[r, c("se_beta0_gs", "se_beta1_gs")] = sqrt(diag(vcov(fit_gs))) ## and its standard error
      ## Error model
      fit_gs = glm(formula = binX ~ 1, 
                   family = binomial,
                   data = dat, 
                   subset = binXstar == 1)
      sett_res[r, "eta0_gs"] = coefficients(fit_gs) ## estimated log prevalence ratio
      sett_res[r, "se_eta0_gs"] = sqrt(diag(vcov(fit_gs))) ## and its standard error
      
      # Fit the gold standard model
      fit_n = glm(formula = Cases ~ binXstar, 
                  family = poisson,
                  offset = log(P),
                  data = dat)
      sett_res[r, c("beta0_n", "beta1_n")] = coefficients(fit_n) ## estimated log odds ratio
      sett_res[r, c("se_beta0_n", "se_beta1_n")] = sqrt(diag(vcov(fit_n))) ## and its standard error
      
      # Select subset of neighborhoods/rows for map-based measures
      query_rows = sample(x = 1:N, 
                          size = ceiling(fix_q * N), 
                          replace = FALSE)
      
      # Make X NA/missing for rows not in selected subset (query_rows)
      dat[!(dat$id %in% query_rows), "binX"] = NA 
      
      # Fit the complete case models
      ## Analysis model
      fit_cc = glm(formula = Cases ~ binX, 
                   family = poisson,
                   offset = log(P),
                   data = dat)
      sett_res[r, c("beta0_cc", "beta1_cc")] = coefficients(fit_cc) ## estimated log prevalence ratio
      sett_res[r, c("se_beta0_cc", "se_beta1_cc")] = sqrt(diag(vcov(fit_cc))) ## and its standard error
      ## Error model
      fit_cc = glm(formula = binX ~ 1, 
                   family = binomial,
                   data = dat, 
                   subset = binXstar == 1)
      sett_res[r, "eta0_cc"] = coefficients(fit_cc) ## estimated odds ratio
      sett_res[r, "se_eta0_cc"] = sqrt(diag(vcov(fit_cc))) ## and its standard error
      
      # Fit the MLE for both models at once 
      fit_mle = suppressMessages(
        mlePossum2(analysis_formula = Cases ~ binX + offset(log(P)), 
                   error_formula = binX ~ binXstar, 
                   data = dat, 
                   noFN = TRUE, 
                   noSE = FALSE)
      )
      sett_res[r, c("beta0_mle", "beta1_mle")] = fit_mle$coefficients$coeff ## estimated log prevalence ratio
      sett_res[r, c("se_beta0_mle", "se_beta1_mle")] = fit_mle$coefficients$se ## and its standard error
      sett_res[r, "eta0_mle"] = fit_mle$misclass_coefficients$coeff ## estimated log odds ratio
      sett_res[r, "se_eta0_mle"] = fit_mle$misclass_coefficients$se ## and its standard error
      
      # Save results
      write.csv(x = sett_res,
                file = paste0("try_mlePossum2/proximity_N", N, "_ppv", 100 * ppv, "_seed", sim_seed, ".csv"), 
                row.names = F)
    }
  }
  toc() ## End runtime for sims with current sample size N
}

# Timing from tictoc:
# Sims with N = 390: 503.197 sec elapsed
# Sims with N = 2200: 1789.26 sec elapsed