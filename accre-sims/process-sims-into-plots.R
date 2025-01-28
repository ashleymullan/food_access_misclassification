#Load libraries
library(tidyverse)
library(kableExtra)
library(latex2exp)
library(ggthemes)

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

#Pivot results for plotting
#write pivot function
pivot <- function(df){
  df_long <- df |>
    pivot_longer(cols = starts_with("beta1_"),
                 names_to = "method_type",
                 values_to = "betahat1")
  return(df_long)
}

#pivot the data
vppv_long <- pivot(results_list$ppv)
vq_long <- pivot(results_list$q)


#define base plot function
plot_base <- function(df_long, x_name, xlab, truth = FALSE){
  df <- data.frame(x = df_long |> pull(x_name) |> as.factor(),
                   y = df_long |> pull(betahat1),
                   n = df_long |> pull("N") |> as.character(),
                   method_type = case_when(
                     df_long$method_type == "beta1_gs" ~ "Gold Standard",
                     df_long$method_type == "beta1_cc" ~ "Complete Case",
                     df_long$method_type == "beta1_mle" ~ "MLE",
                     df_long$method_type == "beta1_n" ~ "Naive"))
  df <- df |>
    #filter(abs(y) < 3) |>
    mutate(method_type = factor(method_type,
                                levels = c("Naive", "Gold Standard",
                                           "MLE", "Complete Case")),
           n = factor(n, levels = c("390", "2200"),
                      labels = c("N = 390", "N = 2200")))
  p <- df |> ggplot(aes(x = x, y = y)) +
    geom_hline(yintercept = 0.1765, color = "magenta") +
    geom_boxplot(aes(fill = method_type)) +
    facet_wrap(vars(n)) +
    theme_minimal() +
    labs(x = xlab,
         y = TeX("$\\hat{\\beta_1}$"),
         title = "",
         fill = "Method") +
    scale_fill_colorblind() +
    theme(legend.position = "top")
  #add a reference line for all necessary cases (fixed prevalence ratio)
  #if(truth) {
  #  print("made it here")
  #  yint = df$beta1[1]
  #  p <- p + geom_hline(yintercept = as.numeric(yint),
  #                               color = "magenta",
  #                               linetype = "dashed",
  #                               linewidth = 5)}

  #little tweaks for specific axis labeling

  #if(x_name == "q") {
  #  df <- df |> mutate(q_labs = paste0(100 * q, "%"))
  #  p <- p + scale_x_discrete(labels = unique(df$q_labs))}
  #if(x_name == "b0") {p <- p + scale_x_discrete(labels = c('5%', '10%', '35%'))}
  #if(x_name == "b1") {p <- p + scale_x_discrete(labels = c('1.1', '1.17', '1.25'))}
  return(p)
}

#run the plots
vq_plot <- plot_base(vq_long, "q", "Query Percentage", truth = TRUE)
vppv_plot <- plot_base(vppv_long, "ppv", "Positive Predictive Value", FALSE)

#save the plots
ggsave(filename = "vq_plot.pdf",
       plot = vq_plot,
       width = 5,
       height = 3.5,
       units = "in")
ggsave(filename = "vpr_plot.pdf",
       plot = vpr_plot,
       width = 5,
       height = 3.5,
       units = "in")
ggsave(filename = "vb0_plot.pdf",
       plot = vb0_plot,
       width = 5,
       height = 3.5,
       units = "in")
ggsave(filename = "vb1_plot.pdf",
       plot = vb1_plot,
       width = 5,
       height = 3.5,
       units = "in")

