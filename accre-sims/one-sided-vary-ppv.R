
# Simulation Setting: Varied Error (One-Sided)

# Load Libraries
library(possum)#, lib.loc = "~/R/rlib") ## for MLE
library(tictoc) ## to calculate runtime

# Set working directory
local_string <- "/Users/ashleymullan/Documents/Food-Access/"
repo_string <- "food_access_misclassification/accre-sims/check-me/"

setwd(paste0(local_string, repo_string))

# Random seed to be used for each simulation setting
sim_seed <- 1031

# Number of replicates per simulation setting
num_reps <- 1000

# Set parameters that won't be varied in the loop
## These values will be set as the defaults in the sim_data() function for convenience
f_beta0 <- -2.2778541  ## outcome model intercept (leads to ~ 11% prevalence)
f_beta1 <- 0.1764544 ## log prevalence ratio for X on Y (given Z)
f_beta2 <- 0.1413503 ## log prevalence ratio for Z on Y (given X)
f_eta1 <- 0.3897812 ## log odds ratio for Z on X (among X* = 1)
f_q <- 0.1 ## proportion of neighborhoods queried

# ---------------------------------------------------------------------------------
# Function to simulate data (arguments defined as follows) ------------------------
## N = number of neighborhoods (sample size)
## beta0 = model intercept
## beta1 = log prevalence ratio of access (given land area) on the outcome
## beta2 = log prevalence ratio of land area (given access) on the outcome
## eta1 = log odds ratio of land area on access
# ---------------------------------------------------------------------------------
sim_data <- function(N, beta0 = f_beta0, beta1 = f_beta1, beta2 = f_beta2,
                    eta0, eta1 = f_eta1) {
  ## Simulate land area Z
  Z <- rgamma(n = N,
             shape = 0.6,
             scale = 0.2)

  ## Simulate straight-line access X*|Z
  binXstar <- rbinom(n = N,
                    size = 1,
                    prob = 1 / (1 + exp(- (1 - Z))))

  ## Simulate map-based access X|X*,Z
  binX <- binXstar
  binX[binXstar == 1] <- rbinom(n = sum(binXstar == 1),
                               size = 1,
                               prob = 1 / (1 + exp(- (eta0 + eta1 * Z[binXstar == 1]))))

  ## Simulate population O (offset)
  P <- rpois(n = N,
            lambda = 4165)

  ## Simulate cases of health outcome Y|X,Z,O
  lambda <- exp(beta0 + beta1 * binX + beta2 * Z)
  Cases <- rpois(n = N,
                lambda = P * lambda)

  ## Create dataset
  dat <- data.frame(id = 1:N, binX, binXstar, P, Cases, Z)

  ## Return dataset
  return(dat)
}

# Set parameters to vary in the loop
ppv_list <- seq(0.5, 0.9, 0.1) ## list of candidate positive predictive values
N_list <- c(390,2200) ## small (Piedmont Triad) and large (NC) full sample sizes

# Loop over all settings
for (N in N_list) { ## iterate through sample sizes
    for (ppv in ppv_list){ ## iterate through sample sizes

      tic(paste("Sims with N =", N,
                "and ppv =", ppv)) ## start counting runtime
      set.seed(sim_seed) ## set random seed (reproducibility rocks!)
      eta0 <- log(ppv/(1 - ppv)) ## compute intercept of error model

      results <- data.frame( ## initialize results container
        sim = paste0("seed",sim_seed, "-run",1:num_reps), ## run label
        N, beta0 = f_beta0, beta1 = f_beta1, beta2 = f_beta2, q = f_q, eta1 = f_eta1, ## fixed simulation setting
        eta0 = eta0, ppv = ppv, queried_ppv = NA, avg_prev = NA, fpr = NA, ## simulation setting
        beta0_gs = NA, se_beta0_gs = NA, beta1_gs = NA, se_beta1_gs = NA, beta2_gs = NA, se_beta2_gs = NA, ## gold standard analysis (outcome)
        beta0_n = NA, se_beta0_n = NA, beta1_n = NA, se_beta1_n = NA, beta2_n = NA, se_beta2_n = NA, ## naive analysis (outcome)
        beta0_cc = NA, se_beta0_cc = NA, beta1_cc = NA, se_beta1_cc = NA, beta2_cc = NA, se_beta2_cc = NA, ## complete case analysis (outcome)
        beta0_mle = NA, se_beta0_mle = NA, beta1_mle = NA, se_beta2_mle = NA, beta2_mle = NA, se_beta1_mle = NA, ## MLE analysis (outcome)
        eta0_mle = NA, se_eta0_mle = NA, eta1_mle = NA, se_eta1_mle = NA, mle_msg = "", ## MLE analysis (misclassification)
        sim_msg = "" ## convergence message
      )

      # Loop over replicates
      for (r in 1:num_reps) {

        # Generate data
        dat <- sim_data(N  = N, ## sample size
                     eta0 = eta0) ## intercept in misclassification mechanism

        # Save average neighborhood prevalence
        results$avg_prev[r] <- mean(dat$Cases / dat$P)

        # Save false positive rate
        results$fpr[r] <- sum(dat$binX == 0 & dat$binXstar == 1) /
           (sum(dat$binX == 0 & dat$binXstar == 1) +
           sum(dat$binX == 0 & dat$binXstar == 0))

        # Fit the gold standard models
        ## Analysis model
        fit_gs <- glm(formula = Cases ~ binX + Z,
                      family = poisson,
                      offset = log(P),
                      data = dat)
        results[r, c("beta0_gs", "beta1_gs", "beta2_gs")] = coefficients(fit_gs) ## estimated log prevalence ratio
        results[r, c("se_beta0_gs", "se_beta1_gs", "se_beta2_gs")] = sqrt(diag(vcov(fit_gs))) ## and its standard error

        # Fit the naive model
        fit_n = glm(formula = Cases ~ binXstar + Z,
                    family = poisson,
                    offset = log(P),
                    data = dat)
        results[r, c("beta0_n", "beta1_n", "beta2_n")] = coefficients(fit_n) ## estimated log odds ratio
        results[r, c("se_beta0_n", "se_beta1_n", "se_beta2_n")] = sqrt(diag(vcov(fit_n))) ## and its standard error

        # Reduce data to only rows with X* = 1
        dat_reduced <- dat[dat$binXstar == 1,]

        # Prepare for query
        query_rows = sample(x = dat_reduced$id, ## can only sample from IDs with X* = 1
                            size = ceiling(f_q * N),
                            replace = FALSE)

        # Make X NA/missing for rows not in selected subset (query_rows)
        dat[!(dat$id %in% query_rows), "binX"] = NA ## make missing

        # Add message if X = 1 for all X* = 1 in queried subset
        queried_ppv = sum(dat$binX == 1 & dat$binXstar == 1, na.rm = TRUE) /
        sum(dat$binXstar == 1 & !is.na(dat$binX), na.rm = TRUE)
        results[r, "queried_ppv"] = queried_ppv
        results[r, "sim_msg"] = ifelse(test = queried_ppv == 1,
                                        yes = "100% PPV (Queried)",
                                        no = "")

        # Fit the complete case models
        ## Analysis model
        fit_cc = glm(formula = Cases ~ binX + Z,
                     family = poisson,
                     offset = log(P),
                     data = dat) #glm auto-fits the complete case
        results[r, c("beta0_cc", "beta1_cc", "beta2_cc")] = coefficients(fit_cc) ## estimated log prevalence ratio
        results[r, c("se_beta0_cc", "se_beta1_cc", "se_beta2_cc")] = sqrt(diag(vcov(fit_cc))) ## and its standard error


        # Fit the MLE (temp if/else to catch queried ppv = 1 issue)
        if(queried_ppv != 1) {
          fit_mle <- suppressMessages(
            mlePossum(analysis_formula = Cases ~ binX + Z + offset(log(P)),
                      error_formula = binX ~ binXstar + Z,
                      beta_init = "Complete-data",
                      eta_init = "Complete-data",
                      data = dat,
                      noSE = FALSE,
                      alternative_SE = TRUE)
            )
          results[r, c("beta0_mle", "beta1_mle", "beta2_mle")] = fit_mle$coefficients[,"Estimate"] ## estimated log prevalence ratio
          results[r, c("se_beta0_mle", "se_beta1_mle", "se_beta2_mle")] = fit_mle$coefficients[,"Std. Error"] ## and its standard error
          results[r, c("eta0_mle", "eta1_mle")] = fit_mle$misclass_coefficients[,"Estimate"] ## estimated log odds ratio
          results[r, c("se_eta0_mle", "se_eta1_mle")] = fit_mle$misclass_coefficients[,"Std. Error"] ## and its standard error
          results[r, "mle_msg"] = fit_mle$converged_msg
        }
        else{
          results[r, c("beta0_mle", "beta1_mle", "beta2_mle")] = results[r, c("beta0_n", "beta1_n", "beta2_n")] ## estimated log prevalence ratio
          results[r, c("se_beta0_mle", "se_beta1_mle", "se_beta2_mle")] = results[r, c("se_beta0_n", "se_beta1_n", "se_beta2_n")] ## and its standard error
          results[r, "mle_msg"] = "naive used as replacement, edge case"
        }



      } ## end loop over the reps
        # Save results from all reps
        write.csv(x = results,
                  file = paste0("varyPPV_N", N, "_ppv", 100 * ppv, "_seed", sim_seed, ".csv"),
                  row.names = F)
        toc() ## End runtime for sims with current N and ppv
    } ## end loop over ppv
} ## end loop over N

# Timing from tictoc:
# Sims with N = 390: 712.096 sec elapsed
# Sims with N = 2200: 2504.638 sec elapsed

# Create combined results file
# Read simulation data in from GitHub repo
all_files = paste0(list.files())
files = all_files[grepl("varyPPV_", all_files)]
res = do.call(what = rbind,
              args = lapply(X = files,
                            FUN = read.csv)
)

write.csv(x = res,
          file = paste0("one-sided-vary-ppv.csv"),
          row.names = F)

