#Load libraries
library(tidyverse)
library(kableExtra)

#Set working directory

##replace local_string with where you keep your github repo directory
local_string <- "/Users/ashleymullan/Documents/Food-Access/"
repo_string <- "food_access_misclassification/accre-sims/"
setwd(paste0(local_string, repo_string))


#Create list of sim results
base <- "one-sided-vary-"
extension <- ".csv"
num_variations <- 2 #change me as I run more sims
varied <- c("ppv", "q") #change me as I run more sims

files <- rep(NA, times = num_variations)
for(i in 1:num_variations){
  files[i] <- paste0(base, varied[i], extension)
}

results_list <- list(ppv = read.csv(files[1]),
                     q = read.csv(files[2]))

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
    dplyr::group_by(N, {{varied_var}}) |>
    dplyr::summarize(bias_gs = mean((beta1_gs - beta1) / beta1), ese_gs = sd(beta1_gs),
              bias_n = mean((beta1_n - beta1) / beta1), ese_n = sd(beta1_n),
              bias_cc = mean((beta1_cc - beta1) / beta1), ese_cc = sd(beta1_cc),
              bias_mle = mean((beta1_mle - beta1) / beta1), ese_mle = sd(beta1_mle),
              ase_mle = mean(se_beta1_mle),
              cp_mle = cp(est = beta1_mle, se = se_beta1_mle, truth = beta1)
    ) |>
    dplyr::mutate(re_cc = (ese_gs ^ 2) / (ese_cc ^ 2),
                  re_mle = (ese_gs ^ 2) / (ese_mle ^ 2))
}

#test dropped reps for varied ppv

results_list$ppv |> nrow() #10K replicates
results_list$ppv |>
  select(-c(mle_msg, sim_msg)) |>
  drop_na() |>
  nrow() #9976 reps (24 dropped)

#visual inspection: not needed
results_list$ppv |>
  beta1_summary(ppv) |>
  select(ese_mle, ase_mle, N, ppv) |>
  round(digits = 5) |> View()

#test dropped reps for varied q

results_list$q |> nrow() #8K replicates
results_list$q |>
  select(-c(mle_msg, sim_msg)) |>
  drop_na() |>
  nrow() #8000 reps (0 dropped)

results_list$q |>
  beta1_summary(q) |>
  round(digits = 5) |>
  select(ese_mle, ase_mle, N, q) |>
  View()


## save summary metrics
beta1_summary(results_list$q, q) |> head() |> View()


#function to format numbers for LaTeX
format_num = function(num, digits = 3) {
  paste0("$", format(round(num, 3), nsmall = digits), "$")
}

## Format for LaTeX
q_table <- results_list$q |>
  beta1_summary(q) |> #head() |> dim()
  mutate(across(c(3:13),
                ~ format_num(.x, digits = 3)))  |>
  mutate(across(c(2),
                ~ format_num(.x, digits = 2)))

ppv_table <- results_list$ppv |>
  beta1_summary(ppv) |>
  mutate(across(c(3:13),
                ~ format_num(.x, digits = 3)))  |>
  mutate(across(c(2),
                ~ format_num(.x, digits = 1)))


#change column names
colnames(q_table) = c("$\\pmb{N}$", "$\\pmb{Q}$", #GGNNCCMM MMCM
                                  rep(c("Bias", "ESE"), times = 4),
                                  "ASE","CP", "RE", "RE")
colnames(ppv_table) = c("$\\pmb{N}$", "$\\pmb{PPV}$", #GGNNCCMM MMCM
                      rep(c("Bias", "ESE"), times = 4),
                      "ASE","CP", "RE", "RE")

#move the relative efficiency of the complete case over
q_table <- q_table[,c(1:8,13,9:12,14)] #input to spit_tex_code
ppv_table <- ppv_table[,c(1:8,13,9:12,14)] #input to spit_tex_code

#spit out TeX code
spit_tex_code <- function(prepped_table){
  prepped_table |>
    kable(format = "latex",
          booktabs = TRUE,
          escape = FALSE,
          align = "ccrcrcrccrcccc") |>
    kable_styling() |>
    add_header_above(header = c(" " = 2,
                                "Gold Standard" = 2,
                                "Naive" = 2,
                                "Complete Case" = 3, #bias, ese, re
                                "MLE" = 5), #bias, ese, ase, re
                     bold = TRUE) |>
    row_spec(row = 0, bold = TRUE)
}

spit_tex_code(q_table)

spit_tex_code(ppv_table)