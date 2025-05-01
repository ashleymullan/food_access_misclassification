#Load libraries
library(tidyverse)
library(kableExtra)

#Set working directory

##replace local_string with where you keep your github repo directory locally
local_string <- "/Users/ashleymullan/Documents/Food-Access/"
repo_string <- "food_access_misclassification/Simulations/"
setwd(paste0(local_string, repo_string))


#Create list of sim results
base1 <- "one-sided-vary-"
base2 <- "two-sided-vary-"
extension <- ".csv"
num_variations1 <- 4 #one sided variations (change as needed)
num_variations2 <- 1 #two sided variations (change as needed)
total_variations <- num_variations1 + num_variations2
varied <- c("ppv", "q", "prev", "prevrat")

files <- rep(NA, times = total_variations) #set up the list of files to read

for(i in 1:num_variations1){ #build the file names for the one sided variations
  files[i] <- paste0(base1, varied[i], extension)
}
for(i in 1:num_variations2){ #build the file names for the two sided variations
  files[num_variations1 + i] <- paste0(base2, varied[i], extension)
}

#read all results files into one place
results_list <- vector(mode = "list", length = total_variations)
for(i in 1:total_variations) {
  results_list[[i]] <- read.csv(files[i])
  names(results_list)[i] <- files[i]}



#compute coverage proportion for Wald interval
cp <- function(est, se, truth) {
  mean((est - 1.96 * se) <= truth & truth <= (est + 1.96 * se))
}

#function to create summary table

## args:
### results_df is the data frame storing the results from the sim settings of interest
### varied_var is the name of the variable you're varying

beta1_summary <- function(results_df, varied_var){
  results_df |>
    dplyr::select(c(N, {{varied_var}}, beta1,
             beta1_gs, beta1_n, beta1_cc, beta1_mle,
             se_beta1_mle)) |>
    tidyr::drop_na() |>
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

#Function to check for dropped replications

## args:
### results_df is the data frame storing the results from the sim settings of interest
### varied_var is the name of the variable you're varying## args:

check_drops <- function(results_df, varied_var){

  full_nrow <- results_df |> nrow()
  check_nrow <- results_df |>
    dplyr::select(c(N, {{varied_var}}, beta1,
                    beta1_gs, beta1_n, beta1_cc, beta1_mle,
                    se_beta1_mle)) |>
    tidyr::drop_na() |> nrow()
  drops <- full_nrow - check_nrow
  perc <- round(drops / full_nrow * 100, 2)

  print(paste0("You originally had ", full_nrow, " replicates but dropped ",
               drops, ", which is ", perc, "%."))
}


#check the drops (not needed but helpful)
for(i in 1:total_variations){
  #grab variable being varied
  var <- substring(files[i], 16) |>
    stringi::stri_reverse() |>
    substring(5) |>
    stringi::stri_reverse()
  if(var == "prev") var <- "beta0"
  if(var == "prevrat") var <- "beta1"
  check_drops(results_list[[i]], var)
}

#function to format numbers for LaTeX
format_num = function(num, digits = 3) {
  paste0("$", format(round(num, 3), nsmall = digits), "$")
}

## Format for LaTeX
q_table <- results_list[[2]] |>
  beta1_summary(q) |>
  mutate(across(c(3:13),
                ~ format_num(.x, digits = 3)))  |>
  mutate(across(c(2),
                ~ format_num(.x, digits = 2)))

ppv_table <- results_list[[1]] |>
  beta1_summary(ppv) |>
  mutate(across(c(3:13),
                ~ format_num(.x, digits = 3)))  |>
  mutate(across(c(2),
                ~ format_num(.x, digits = 1)))

prev_table <- results_list[[3]] |>
  beta1_summary(beta0) |>
  mutate(across(c(3:13),
                ~ format_num(.x, digits = 3)))  |>
  mutate(across(c(2),
                ~ format_num(.x, digits = 1)))

prevrat_table <- results_list[[4]] |>
  beta1_summary(beta1) |>
  mutate(across(c(3:13),
                ~ format_num(.x, digits = 3)))  |>
  mutate(across(c(2),
                ~ format_num(.x, digits = 1)))

tppv_table <- results_list[[5]] |>
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
colnames(prev_table) = c("$\\pmb{N}$", "$\\pmb{\\beta_0}$", #GGNNCCMM MMCM
                        rep(c("Bias", "ESE"), times = 4),
                        "ASE","CP", "RE", "RE")
colnames(prevrat_table) = c("$\\pmb{N}$", "$\\pmb{\beta_1}$", #GGNNCCMM MMCM
                         rep(c("Bias", "ESE"), times = 4),
                         "ASE","CP", "RE", "RE")
colnames(tppv_table) = c("$\\pmb{N}$", "$\\pmb{PPV}$", #GGNNCCMM MMCM
                            rep(c("Bias", "ESE"), times = 4),
                            "ASE","CP", "RE", "RE")

#move the relative efficiency of the complete case over
q_table <- q_table[,c(1:8,13,9:12,14)] #input to spit_tex_code
ppv_table <- ppv_table[,c(1:8,13,9:12,14)] #input to spit_tex_code
prev_table <- prev_table[,c(1:8,13,9:12,14)] #input to spit_tex_code
prevrat_table <- prevrat_table[,c(1:8,13,9:12,14)] #input to spit_tex_code
tppv_table <- tppv_table[,c(1:8,13,9:12,14)] #input to spit_tex_code

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


#spit out LaTeX code
spit_tex_code(q_table)

spit_tex_code(ppv_table)

spit_tex_code(prev_table)

spit_tex_code(prevrat_table)

spit_tex_code(tppv_table)
