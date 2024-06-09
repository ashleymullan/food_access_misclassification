# Install libraries
## Run once
## if needed: install.packages("devtools") 
## then: 
## devtools::install_github("sarahlotspeich/possum", 
##                          ref = "main") 

# Load libraries
library(possum) ## for SMLE
library(auditDesignR) ## for subsample designs
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
fix_beta0 = -2.2456410 ## outcome model intercept (leads to ~ 7% prevalence)
fix_beta1 = 0.1549768 ## log prevalence ratio for X on Y
fix_muU = -0.7 ## error mean
fix_sigmaU = 0.8 ## error standard deviation
fix_access = 1 ## mile threshold to discretize X into access/no access
fix_q = 0.1 ## proportion of neighborhoods queried

# ---------------------------------------------------------------------------------
# Function to simulate data (arguments defined as follows) ------------------------
## N = number of neighborhoods (sample size)
## beta0 = model intercept
## beta1 = log prevalence ratio 
## muU = mean of the measurement error distribution
## sigmaU = standard deviation of the measurement error distribution
## access = threshold (in miles) used to discretize proximity into access/no access
# ---------------------------------------------------------------------------------
sim_data = function(N, beta0 = fix_beta0, beta1 = fix_beta1, muU = fix_muU, sigmaU = fix_sigmaU, access = fix_access) {
  ## Simulate true (map-based) proximity to grocery store
  X = rgamma(n = N,
             shape = 1,
             scale = 2.5) 
  
  ## Construct discretived true (map-based) access
  binX = as.numeric(X <= access)
  
  ## Simulate random errors
  U = truncnorm::rtruncnorm(n = N, 
                            a = -Inf, 
                            b = 0, 
                            mean = muU, 
                            sd = sigmaU)
  
  ## Construct error-prone (straight-line) proximity to grocery store
  Xstar = X + U ### assuming additive measurement error model
  
  ## Construct discretived error-prone (straight-line) access
  binXstar = as.numeric(Xstar <= access)
  
  ## Simulate population
  P = rpois(n = N, 
            lambda = 4165)
  
  ## Simulate cases of health outcome
  lambda = exp(beta0 + beta1 * binX)
  Cases = rpois(n = N, 
                lambda = P * lambda)
  
  ## Create dataset
  dat = data.frame(id = 1:N, X, Xstar, binX, binXstar, P, Cases)
  
  # Return dataset
  return(dat)
}

# Loop over different sample sizes: N = 390 (Piedmont Triad), 2200 (all of NC)
for (N in c(390, 2200)) {
  tic(paste("Sims with N =", N)) ## 
  fix_n = ifelse(test = N == 390, 
                 yes = 40, 
                 no = ceiling(fix_q * N))
  # And proportion to be queried for complete case/imputation analyses: 0.1, 0.25, 0.5, 0.75
  for (q in c(0.1, 0.25, 0.5, 0.75)){
    # Be reproducible
    set.seed(sim_seed) ## set random seed
    
    # Create dataframe to save results for setting
    sett_res = data.frame(
      sim = paste(sim_seed, 1:num_reps, sep = "-"), 
      N, beta0 = fix_beta0, beta1 = fix_beta1, muU = fix_muU, ## simulation setting
      sigmaU = fix_sigmaU, q = q, avg_prev = NA, fpr = NA, ## simulation setting
      beta0_srs = NA, se_beta0_srs = NA, beta1_srs = NA, se_beta1_srs = NA, ## gold standard analysis
      beta0_bcc = NA, se_beta0_bcc = NA, beta1_bcc = NA, se_beta1_bcc = NA, ## naive analysis
      beta0_cc = NA, se_beta0_cc = NA, beta1_cc = NA, se_beta1_cc = NA, ## complete case analysis
      beta0_resid = NA, se_beta0_resid = NA, beta1_resid = NA, se_beta1_resid = NA ## SMLE analysis
    )
    
    # Loop over replicates 
    for (r in 1:num_reps) {
      # Generate data
      dat = sim_data(N  = N) ## sample size
      
      # Save average neighborhood prevalence
      sett_res$avg_prev[r] = mean(dat$Cases / dat$P)
      
      # Save false positive rate
      sett_res$fpr[r] = sum(dat$binX == 0 & dat$binXstar == 1) /
        (sum(dat$binX == 0 & dat$binXstar == 1) + 
           sum(dat$binX == 0 & dat$binXstar == 0))
      
      # Save positive predictive value
      sett_res$ppv[r] = with(dat, 
                             sum(binX == 1 & binXstar == 1) / 
                               sum(binXstar == 1))
      
      # Design (i): SRS
      ## Select random subset of neighborhoods/rows for map-based measures
      query_rows = sample(x = 1:N, 
                          size = fix_n, 
                          replace = FALSE)
      
      ## Make X NA/missing for rows not in selected subset (query_rows)
      dat[, "binX_partial"] = dat[, "binX"]
      dat[!(dat$id %in% query_rows), "binX_partial"] = NA 
      
      ## Fit the SMLE with the SRS queried subsample
      fit_smle = smlePossum(Y = "Cases", 
                            offset = "P", 
                            X_unval = "binXstar", 
                            X = "binX_partial", 
                            data = dat)
      sett_res[r, c("beta0_srs", "beta1_srs")] = fit_smle$coeff$coeff ## estimated log prevalence ratio
      sett_res[r, c("se_beta0_srs", "se_beta1_srs")] = fit_smle$coeff$se ## and its standard error
      
      # Design (ii): BCC* 
      ## Define strata of above/below state diabetes prevalence
      dat[, "DIABETES_STRAT"] = ifelse(test = dat[, "Cases"] < median(dat[, "Cases"]), 
                                       yes = "Below", 
                                       no = "Above")
      
      ## Select random subset of neighborhoods/rows for map-based measures
      query_indicators = sample_bcc(dat = dat, 
                                    phI = N, 
                                    phII = fix_n, 
                                    sample_on = c("DIABETES_STRAT", "binXstar"))
      
      ## Make X NA/missing for rows not in selected subset (query_rows)
      dat[, "binX_partial"] = dat[, "binX"]
      dat[query_indicators == 0, "binX_partial"] = NA 
      
      ## Fit the SMLE with the SRS queried subsample
      fit_smle = smlePossum(Y = "Cases", 
                            offset = "P", 
                            X_unval = "binXstar", 
                            X = "binX_partial", 
                            data = dat)
      sett_res[r, c("beta0_bcc", "beta1_bcc")] = fit_smle$coeff$coeff ## estimated log prevalence ratio
      sett_res[r, c("se_beta0_bcc", "se_beta1_bcc")] = fit_smle$coeff$se ## and its standard error
      
      # Design (iii): CC* on Xstar
      ## Select random subset of neighborhoods/rows for map-based measures
      query_indicators = sample_cc(dat = dat, 
                                   phI = N, 
                                   phII = fix_n, 
                                   sample_on = "binXstar")
      
      ## Make X NA/missing for rows not in selected subset (query_rows)
      dat[, "binX_partial"] = dat[, "binX"]
      dat[query_indicators == 0, "binX_partial"] = NA 
      
      ## Fit the SMLE with the SRS queried subsample
      fit_smle = smlePossum(Y = "Cases", 
                            offset = "P", 
                            X_unval = "binXstar", 
                            X = "binX_partial", 
                            data = dat)
      sett_res[r, c("beta0_cc", "beta1_cc")] = fit_smle$coeff$coeff ## estimated log prevalence ratio
      sett_res[r, c("se_beta0_cc", "se_beta1_cc")] = fit_smle$coeff$se ## and its standard error
      
      # Design (iv): Residual sampling based on the naive model
      ## Select random subset of neighborhoods/rows for map-based measures
      query_indicators = sample_resid(formula = Cases ~ binXstar + offset(log(P)), 
                                      family = "poisson", 
                                      dat = dat, 
                                      phI = N, 
                                      phII = fix_n)

      ## Make X NA/missing for rows not in selected subset (query_rows)
      dat[, "binX_partial"] = dat[, "binX"]
      dat[query_indicators == 0, "binX_partial"] = NA 
      
      ## Fit the SMLE with the SRS queried subsample
      fit_smle = smlePossum(Y = "Cases", 
                            offset = "P", 
                            X_unval = "binXstar", 
                            X = "binX_partial", 
                            data = dat)
      sett_res[r, c("beta0_resid", "beta1_resid")] = fit_smle$coeff$coeff ## estimated log prevalence ratio
      sett_res[r, c("se_beta0_resid", "se_beta1_resid")] = fit_smle$coeff$se ## and its standard error
      
      # Save results
      write.csv(x = sett_res,
                file = paste0("vary_design/proximity_N", N, "_q", 100 * q, "_seed", sim_seed, ".csv"), 
                row.names = F)
    }
  }
  toc() ## End runtime for sims with current sample size N
}

# Timing from tictoc:
## Sims with N = 390: 154.936 sec elapsed
## Sims with N = 2200: 528.144 sec elapsed