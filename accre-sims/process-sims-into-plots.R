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
plot_base <- function(df_long, x_name, xlab, truth = FALSE, legend_position = "NONE"){
  df <- data.frame(x = df_long |> pull(x_name) |> as.factor(),
                   y = df_long |> pull(betahat1) |> exp(),
                   n = df_long |> pull("N") |> as.character(),
                   method_type = case_when(
                     df_long$method_type == "beta1_gs" ~ "Gold Standard",
                     df_long$method_type == "beta1_cc" ~ "Complete Case",
                     df_long$method_type == "beta1_mle" ~ "MLE",
                     df_long$method_type == "beta1_n" ~ "Naive"))

  ## save true ratio if we have a fixed truth
  ### note: we picked max but it doesn't matter because they all match
  if(truth){true_hr <- df_long |> pull(beta1) |> max() |> exp() }

  df <- df |>
    mutate(method_type = factor(method_type,
                                levels = c("Naive", "Gold Standard",
                                           "MLE", "Complete Case")),
           n = factor(n, levels = c("390", "2200"),
                      labels = c("N = 390", "N = 2200")))

  ## start building the plot
  p <- df |> ggplot(aes(x = x, y = y))

  ## add reference line if we have the fixed truth
  if(truth) {
    p <- p + geom_hline(yintercept = true_hr,
                        color = "gray",
                        linetype = "dashed")
  }

  ## continue building plot
  p <- p +
    geom_boxplot(aes(fill = method_type)) +
    facet_wrap(vars(n)) +
    theme_bw() +
    labs(x = xlab,
         y = TeX("$exp(\\hat{\\beta_1}$)"),
         title = "",
         fill = "Method") +
    scale_fill_colorblind() +
    theme(legend.position = legend_position,
          strip.background = element_rect(fill = "#ffd3fa"),
          plot.margin = unit(c(5,5,5,5), "mm"))
  return(p)
}

#run the plots
vq_plot <- plot_base(vq_long, "q", "Query Percentage", truth = TRUE) +
  scale_x_discrete(labels = c("10%", "25%", "50%", "75%"))
vppv_plot <- plot_base(vppv_long |> ## filter to make graph readable, omits 28 points
                         filter(exp(betahat1) < 1.3) |> ## drops extreme outliers that make graph unreadable
                         filter(exp(betahat1) > 1.07), ## drops extreme outliers that make graph unreadable
                       "ppv", "Positive Predictive Value",
                       truth = TRUE) +
  scale_x_discrete(labels = c("50%", "60%","70%", "80%", "90%"))
vq_plot_no_legend <- plot_base(vq_long, "q", "Query Percentage", truth = TRUE) +
  scale_x_discrete(labels = c("10%", "25%", "50%", "75%"))
vppv_plot <- plot_base(vppv_long |> ## filter to make graph readable, omits 28 points
                         filter(exp(betahat1) < 1.3) |> ## drops extreme outliers that make graph unreadable
                         filter(exp(betahat1) > 1.07), ## drops extreme outliers that make graph unreadable
                       "ppv", "Positive Predictive Value",
                       truth = TRUE) +
  scale_x_discrete(labels = c("50%", "60%","70%", "80%", "90%"))

#throw naive out for a better visual
vq_plot_no_naive <- plot_base(vq_long |> filter(method_type != "beta1_n"),
                              "q", "Query Percentage", truth = TRUE) +
  scale_x_discrete(labels = c("10%", "25%", "50%", "75%")) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73"))
vppv_plot_no_naive <- plot_base(vppv_long |>
                                  filter(method_type != "beta1_n") |>
                                  filter(exp(betahat1) < 1.31) |>
                                  filter(exp(betahat1) > 1.08), ## takes out 28 extreme outliers for better visual
                              "ppv", "Positive Predictive Value", truth = TRUE) +
  scale_x_discrete(labels = c("50%", "60%", "70%", "80%", "90%")) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73"))

#save the plots

save_plot <- function(plot, path = paste0(getwd(), "/plots/")){
  fn <- paste0(deparse(substitute(plot)), ".pdf")
  ggsave(filename = fn,
         path = path,
         width = 5,
         height = 3.5,
         units = "in")
}

save_plot(vq_plot)
save_plot(vq_plot_no_naive)
save_plot(vq_plot_no_legend)
save_plot(vppv_plot)
save_plot(vppv_plot_no_naive)

