# Load libraries
library(possum) ## for SMLE
library(tictoc) ## to calculate runtime

# Set working directory
local_string <- "/Users/ashleymullan/Documents/Food-Access/"
local_string <- "/Users/lotspes/Documents/"
repo_string <- "food_access_misclassification/Simulations/One_Sided_Vary_Q"
setwd(paste0(local_string, repo_string))

# Random seed to be used for each simulation setting
sim_seed <- 1031

# Number of replicates per simulation setting
num_reps = 1000
## Note: You may want to run this code was run in parallel a cluster instead of locally, as it can be slow.

# Set parameters that won't be varied in the loop
## These values will be set as the defaults in the sim_data() function for convenience
fix_beta0 = -2.2778541  ## outcome model intercept (leads to ~ 11% prevalence)
fix_beta1 = 0.1764544 ## log prevalence ratio for X on Y (given Z)
fix_beta2 = 0.1413503 ## log prevalence ratio for Z on Y (given X)
fix_eta1 = 0.3897812 ## log odds ratio for Z on X (among X* = 1)
fix_ppv = 0.6 ## positive predictive value
fix_eta0 = log(fix_ppv/(1 - fix_ppv))
# ---------------------------------------------------------------------------------
# Function to simulate data (arguments defined as follows) ------------------------
## N = number of neighborhoods (sample size)
## beta0 = model intercept
## beta1 = log prevalence ratio of access (given land area) on the outcome
## beta2 = log prevalence ratio of land area (given access) on the outcome
## eta1 = log odds ratio of land area on access
# ---------------------------------------------------------------------------------
sim_data = function(N, beta0 = fix_beta0, beta1 = fix_beta1, beta2 = fix_beta2, eta0 = fix_eta0, eta1 = fix_eta1) {
  ## Simulate land area Z
  Z = rgamma(n = N,
             shape = 0.6,
             scale = 0.2)

  ## Simulate straight-line access X*|Z
  binXstar = rbinom(n = N,
                    size = 1,
                    prob = 1 / (1 + exp(- (1 - Z))))

  ## Simulate map-based access X|X*,Z
  binX = binXstar
  binX[binXstar == 1] = rbinom(n = sum(binXstar == 1),
                               size = 1,
                               prob = 1 / (1 + exp(- (eta0 + eta1 * Z[binXstar == 1]))))

  ## Simulate population O
  P = rpois(n = N,
            lambda = 4165)

  ## Simulate cases of health outcome Y|X,Z,O
  lambda = exp(beta0 + beta1 * binX + beta2 * Z)
  Cases = rpois(n = N,
                lambda = P * lambda)

  ## Create dataset
  dat = data.frame(id = 1:N, binX, binXstar, P, Cases, Z)

  # Return dataset
  return(dat)
}

# Loop over different sample sizes: N = 390 (Piedmont Triad), 2200 (all of NC)
for (N in c(390, 2200)) {
  print(paste("System time =", Sys.time()))
  tic(paste("Sims with N =", N)) ## Start counting runtime for sims with current sample size N
  # And proportion to be queried for complete case/imputation analyses: 0.1, 0.25, 0.5, 0.75
  for (q in c(0.05, 0.1, 0.25, 0.5)){
    # Be reproducible
    set.seed(sim_seed) ## set random seed

    # Create dataframe to save results for setting
    sett_res = data.frame(
      sim = paste(sim_seed, 1:num_reps, sep = "-"),
      N, beta0 = fix_beta0, beta1 = fix_beta1, beta2 = fix_beta2, q = q,  ## simulation setting
      eta0 = fix_eta0, eta1 = fix_eta1, ppv = fix_ppv, queried_ppv = NA, queried_strat_ppv = NA, avg_prev = NA, fpr = NA, ## simulation setting
      beta0_gs = NA, se_beta0_gs = NA, beta1_gs = NA, se_beta1_gs = NA, beta2_gs = NA, se_beta2_gs = NA, ## gold standard analysis (outcome)
      eta0_gs = NA, se_eta0_gs = NA, eta1_gs = NA, se_eta1_gs = NA, ## gold standard analysis (misclassification)
      beta0_n = NA, se_beta0_n = NA, beta1_n = NA, se_beta1_n = NA, beta2_n = NA, se_beta2_n = NA, ## naive analysis (outcome)
      beta0_cc = NA, se_beta0_cc = NA, beta1_cc = NA, se_beta1_cc = NA, beta2_cc = NA, se_beta2_cc = NA, ## complete case analysis (outcome)
      #eta0_cc = NA, se_eta0_cc = NA, eta1_cc = NA, se_eta1_cc = NA, ## complete case analysis (misclassification)
      beta0_mle = NA, se_beta0_mle = NA, beta1_mle = NA, se_beta2_mle = NA, beta2_mle = NA, se_beta1_mle = NA, ## MLE analysis (outcome)
      eta0_mle = NA, se_eta0_mle = NA, eta1_mle = NA, se_eta1_mle = NA, mle_msg = "", ## MLE analysis (misclassification)
      beta0_mle_strat = NA, se_beta0_mle_strat = NA, beta1_mle_strat = NA, se_beta2_mle_strat = NA, beta2_mle_strat = NA, se_beta1_mle_strat = NA, ## MLE (sample only Y*=1) analysis (outcome)
      eta0_mle_strat = NA, se_eta0_mle_strat = NA, eta1_mle_strat = NA, se_eta1_mle_strat = NA, mle_strat_msg = "", ## MLE analysis (sample only Y*=1) (misclassification)
      sim_msg = "", sim_strat_msg = ""
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

      # Fit the gold standard models
      ## Analysis model
      fit_gs = glm(formula = Cases ~ binX + Z,
                   family = poisson,
                   offset = log(P),
                   data = dat)
      sett_res[r, c("beta0_gs", "beta1_gs", "beta2_gs")] = coefficients(fit_gs) ## estimated log prevalence ratio
      sett_res[r, c("se_beta0_gs", "se_beta1_gs", "se_beta2_gs")] = sqrt(diag(vcov(fit_gs))) ## and its standard error
      ## Error model
      fit_gs = glm(formula = binX ~ Z,
                   family = binomial,
                   data = dat,
                   subset = binXstar == 1)
      sett_res[r, c("eta0_gs", "eta1_gs")] = coefficients(fit_gs) ## estimated log prevalence ratio
      sett_res[r, c("se_eta0_gs", "se_eta1_gs")] = sqrt(diag(vcov(fit_gs))) ## and its standard error

      # Fit the gold standard model
      fit_n = glm(formula = Cases ~ binXstar + Z,
                  family = poisson,
                  offset = log(P),
                  data = dat)
      sett_res[r, c("beta0_n", "beta1_n", "beta2_n")] = coefficients(fit_n) ## estimated log odds ratio
      sett_res[r, c("se_beta0_n", "se_beta1_n", "se_beta2_n")] = sqrt(diag(vcov(fit_n))) ## and its standard error

      # Query design 1: Select **random** subset of neighborhoods/rows for map-based measures
      query_rows = sample(x = 1:N,
                          size = ceiling(q * N),
                          replace = FALSE)
      
      # Make X NA/missing for rows not in selected subset (query_rows)
      dat[, "binX_srs"] = dat[, "binX"] ## initialize as binX
      dat[!(dat$id %in% query_rows), "binX_srs"] = NA ## make missing
      
      # Add message if X = 1 for all X* = 1 in queried subset
      queried_ppv = sum(dat$binX_srs == 1 & dat$binXstar == 1, na.rm = TRUE) /
        sum(dat$binXstar == 1 & !is.na(dat$binX_srs), na.rm = TRUE)
      sett_res[r, "queried_ppv"] = queried_ppv
      sett_res[r, "sim_msg"] = ifelse(test = queried_ppv == 1,
                                      yes = "100% PPV (Queried)",
                                      no = "")
      
      # Fit the complete case models
      ## Analysis model
      fit_cc = glm(formula = Cases ~ binX_srs + Z,
                   family = poisson,
                   offset = log(P),
                   data = dat)
      sett_res[r, c("beta0_cc", "beta1_cc", "beta2_cc")] = coefficients(fit_cc) ## estimated log prevalence ratio
      sett_res[r, c("se_beta0_cc", "se_beta1_cc", "se_beta2_cc")] = sqrt(diag(vcov(fit_cc))) ## and its standard error
      ## Error model
      # fit_cc = glm(formula = binX_srs ~ Z,
      #              family = binomial,
      #              data = dat,
      #              subset = binXstar == 1)
      # sett_res[r, c("eta0_cc", "eta1_cc")] = coefficients(fit_cc) ## estimated odds ratio
      # sett_res[r, c("se_eta0_cc", "se_eta1_cc")] = sqrt(diag(vcov(fit_cc))) ## and its standard error
      
      # Fit the MLE for both models at once
      fit_mle = suppressMessages(
        mlePossum(analysis_formula = Cases ~ binX_srs + Z + offset(log(P)),
                  error_formula = binX_srs ~ binXstar + Z,
                  data = dat,
                  noSE = FALSE)
      )
      sett_res[r, c("beta0_mle", "beta1_mle", "beta2_mle")] = fit_mle$coefficients$coeff ## estimated log prevalence ratio
      sett_res[r, c("se_beta0_mle", "se_beta1_mle", "se_beta2_mle")] = fit_mle$coefficients$se ## and its standard error
      sett_res[r, c("eta0_mle", "eta1_mle")] = fit_mle$misclass_coefficients$coeff ## estimated log odds ratio
      sett_res[r, c("se_eta0_mle", "se_eta1_mle")] = fit_mle$misclass_coefficients$se ## and its standard error
      sett_res[r, "mle_msg"] = fit_mle$converged_msg
      
      # Query design 1: Select **stratified random** subset of neighborhoods/rows for map-based measures
      ## Take as many with X* = 1 as possible
      query_rows = sample(x = which(dat$binXstar == 1), 
                          size = min(length(which(dat$binXstar == 1)), 
                                     ceiling(q * N)),
                          replace = FALSE)
      ## And if there are fewer than n, take the rest from X* = 0 
      ## Just to make it "fair" in each simulation (to have the same n)
      if (length(query_rows) < ceiling(q * N)) {
        query_rows = c(query_rows, 
                       sample(x = which(dat$binXstar == 0), 
                              size = (ceiling(q * N) - length(which(dat$binXstar == 1))),
                              replace = FALSE))
      }
      
      # Make X NA/missing for rows not in selected subset (query_rows)
      dat[, "binX_strat"] = dat[, "binX"] ## initialize as binX
      dat[!(dat$id %in% query_rows), "binX_strat"] = NA ## make missing
      
      # Add message if X = 1 for all X* = 1 in queried subset
      queried_ppv = sum(dat$binX_strat == 1 & dat$binXstar == 1, na.rm = TRUE) /
        sum(dat$binXstar == 1 & !is.na(dat$binX_strat), na.rm = TRUE)
      sett_res[r, "queried_strat_ppv"] = queried_ppv
      sett_res[r, "sim_strat_msg"] = ifelse(test = queried_ppv == 1,
                                            yes = "100% PPV (Queried)",
                                            no = "")
      
      if (queried_ppv < 1) {
        # Fit the MLE for both models at once
        fit_mle = suppressMessages(
          mlePossum(analysis_formula = Cases ~ binX_strat + Z + offset(log(P)),
                    error_formula = binX_strat ~ binXstar + Z,
                    data = dat,
                    noSE = FALSE)
        )
        sett_res[r, c("beta0_mle_strat", "beta1_mle_strat", "beta2_mle_strat")] = fit_mle$coefficients$coeff ## estimated log prevalence ratio
        sett_res[r, c("se_beta0_mle_strat", "se_beta1_mle_strat", "se_beta2_mle_strat")] = fit_mle$coefficients$se ## and its standard error
        sett_res[r, c("eta0_mle_strat", "eta1_mle_strat")] = fit_mle$misclass_coefficients$coeff ## estimated log odds ratio
        sett_res[r, c("se_eta0_mle_strat", "se_eta1_mle_strat")] = fit_mle$misclass_coefficients$se ## and its standard error
        sett_res[r, "mle_strat_msg"] = fit_mle$converged_msg  
      }

      # Save results
      write.csv(x = sett_res, 
                file = paste0("varyQ_N", N, "_q", 100 * q, "_seed", sim_seed, ".csv"),
                row.names = F)
    }
  }
  toc() ## End runtime for sims with current sample size N
}

# Timing from tictoc:
# Sims with N = 390: 518.39 sec elapsed
# Sims with N = 2200: 1872.032 sec elapsed

# Create combined results file
# Read simulation data in from GitHub repo
all_files = paste0(list.files())
files <- all_files[grepl("varyQ", all_files)]
res = do.call(what = rbind,
              args = lapply(X = files,
                            FUN = read.csv)
)

write.csv(x = res,
          file = paste0("one_sided_vary_Q.csv"),
          row.names = F)


