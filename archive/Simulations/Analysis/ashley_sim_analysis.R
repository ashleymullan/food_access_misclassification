#Load libraries
library(tidyverse)
library(kableExtra)

#Set working directory

##replace local_string with where you keep your github repo directory
local_string <- "/Users/ashleymullan/Documents/Food-Access/"
repo_string <- "food_access_misclassification/Simulations/Ashley_Sims"
setwd(paste0(local_string, repo_string))


#Create list of sim results
base <- "one_sided_vary_"
extension <- ".csv"
num_variations <- 4
varied <- c("ppv", "prevalence", "Q", "PR")

files <- rep(NA, times = num_variations)
for(i in 1:num_variations){
  files[i] <- paste0(base, varied[i], extension)
}

results_list <- list(ppv = read.csv(files[1]),
     prevalence = read.csv(files[2]),
     q = read.csv(files[3]),
     pr = read.csv(files[4]))

#compute coverage proportion
cp = function(est, se, truth) {
  mean((est - 1.96 * se) <= truth & truth <= (est + 1.96 * se))
}

#function to create summary table

## args:
### results_df is the data frame storing the results from the sim settings of interest
### varied_var is the name of the variable you're varying

beta1_summary <- function(results_df, varied_var){
  results_df |>
    select(-c(mle_msg, sim_msg)) |>
    drop_na() |>
    group_by(N, {{varied_var}}) |>
    summarize(bias_gs = mean((beta1_gs - beta1) / beta1), ese_gs = sd(beta1_gs),
              bias_n = mean((beta1_n - beta1) / beta1), ese_n = sd(beta1_n),
              bias_cc = mean((beta1_cc - beta1) / beta1), ese_cc = sd(beta1_cc),
              bias_mle = mean((beta1_mle - beta1) / beta1), ese_mle = sd(beta1_mle),
              ase_mle = mean(se_beta1_mle),
              cp_mle = cp(est = beta1_mle, se = se_beta1_mle, truth = beta1)
    ) |>
    dplyr::mutate(re_cc = (ese_gs ^ 2) / (ese_cc ^ 2),
                  re_mle = (ese_gs ^ 2) / (ese_mle ^ 2))
}

#test for varied ppv

results_list$ppv |> nrow() #10K replicates
results_list$ppv |>
  select(-c(mle_msg, sim_msg)) |>
  drop_na() |>
  nrow() #9990 reps (10 dropped)

results_list$ppv |>
  beta1_summary(ppv) |>
  round(digits = 4) |> View()

#test for varied q

results_list$q |> nrow() #8K replicates
results_list$q |>
  select(-c(mle_msg, sim_msg)) |>
  drop_na() |>
  nrow() #7999 reps (1 dropped)

results_list$q |>
  beta1_summary(q) |>
  round(digits = 4) |> View()

#test for varied prevalence

results_list$prevalence |> nrow() #8K replicates
results_list$prevalence |>
  select(-c(mle_msg, sim_msg)) |> #get rid of columns that we don't care about NA
  drop_na() |> nrow() #7997 reps (3 dropped)

results_list$prevalence |>
  beta1_summary(beta0) |>
  round(digits = 4) |> View()

#test for varied PR

results_list$pr |> nrow() #10K replicates
results_list$pr |>
  select(-c(mle_msg, sim_msg)) |>
  drop_na() |>
  nrow() #7995 reps (5 dropped)

results_list$pr |>
  beta1_summary(beta1) |>
  round(digits = 4) |> View()

##ASHLEY DEAL WITH LATER

#function to format numbers for LaTeX
format_num = function(num, digits = 3) {
  paste0("$", format(round(num, 3), nsmall = digits), "$")
}

## Format for LaTeX
full_result_summary = full_result_summary |>
  mutate_at(.vars = 2:18, .funs = format_num, digits = 3) |>
  mutate_at(.vars = 2:4, .funs = format_num, digits = 2)

#change col names
colnames(full_result_summary) = c("$\\pmb{N}$", "$\\pmb{Q}$",
                                  "FPR", "TPR",
                                  "\\pmb{\\beta_0}", "\\pmb{\\beta_1}",
                                  rep(c("Bias", "ESE"), times = 4),
                                  "ASE","CP", "RE", "RE")
#GGNNCCMM MMCM

full_result_summary |>
  kable(format = "latex",
        booktabs = TRUE,
        escape = FALSE,
        align = "cccccccccccccccccc") |>
  kable_styling() |>
  add_header_above(header = c(" " = 6, "Gold Standard" = 2,
                              "Naive" = 2, "Complete Case" = 2,
                              "MLE" = 4, "Complete Case" = 1, "MLE" = 1),
                   bold = TRUE) |>
  row_spec(row = 0, bold = TRUE)
## And a \multicolumn used to separate the three parameters





