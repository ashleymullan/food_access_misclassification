# Load packages
library(ggplot2) ## for pretty plots
library(dplyr) ## for data wrangling

# Read simulation data in from GitHub repo
res = read.csv("https://raw.githubusercontent.com/ashleymullan/food_access_misclassification/main/Simulations/try_mlePossum2_with_Z_seed114.csv")

# Summarize convergence messages for MLE
res |> 
  group_by(N, ppv, mle_msg) |> 
  summarize(num = n()) |> 
  arrange(desc(mle_msg))

# Exclude few reps were probabilities numeric 0/1 occurred
res = res |> 
  filter(mle_msg != "Fitted probabilities numerically 0 or 1 at convergence")

# Create summary of results
res |> 
  mutate(eta0 = log(ppv / (1 - ppv))) |> 
  group_by(N, ppv) |> 
  summarize(avg_fpr = mean(fpr), 
            bias_beta0_gs = mean(beta0_gs - beta0), 
            ese_beta0_gs = sd(beta0_gs), 
            bias_beta1_gs = mean(beta1_gs - beta1), 
            ese_beta1_gs = sd(beta1_gs), 
            bias_beta2_gs = mean(beta2_gs - beta2), 
            ese_beta2_gs = sd(beta2_gs), 
            bias_eta0_gs = mean(eta0_gs - eta0), 
            bias_eta1_gs = mean(eta1_gs - eta1), 
            bias_beta0_n = mean(beta0_n - beta0), 
            ese_beta0_n = sd(beta0_n), 
            bias_beta1_n = mean(beta1_n - beta1), 
            ese_beta1_n = sd(beta1_n), 
            bias_beta2_n = mean(beta2_n - beta2), 
            ese_beta2_n = sd(beta2_n), 
            bias_beta0_cc = mean(beta0_cc - beta0), 
            ese_beta0_cc = sd(beta0_cc), 
            bias_beta1_cc = mean(beta1_cc - beta1), 
            ese_beta1_cc = sd(beta1_cc), 
            re_beta1_cc = var(beta1_gs) / var(beta1_cc),
            bias_beta2_cc = mean(beta2_cc - beta2), 
            ese_beta2_cc = sd(beta2_cc), 
            bias_eta0_cc = mean(eta0_cc - eta0), 
            bias_eta1_cc = mean(eta1_cc - eta1), 
            bias_beta0_mle = mean(beta0_mle - beta0), 
            ese_beta0_mle = sd(beta0_mle), 
            ase_beta0_mle = mean(se_beta0_mle, 
                                 na.rm = TRUE), 
            bias_beta1_mle = mean(beta1_mle - beta1), 
            ese_beta1_mle = sd(beta1_mle), 
            ase_beta1_mle = mean(se_beta1_mle, 
                                 na.rm = TRUE), 
            cp_beta1_mle = mean((beta1_mle - 1.96 * se_beta1_mle) <= beta1 & 
                                   beta1 <= (beta1_mle + 1.96 * se_beta1_mle), 
                                na.rm = TRUE), 
            re_beta1_mle = var(beta1_gs) / var(beta1_mle),
            bias_beta2_mle = mean(beta2_mle - beta2), 
            ese_beta2_mle = sd(beta2_mle, 
                               na.rm = TRUE), 
            ase_beta2_mle = mean(se_beta2_mle, 
                                 na.rm = TRUE), 
            bias_eta0_mle = mean(eta0_mle - eta0), 
            ese_eta0_mle = sd(eta0_mle), 
            ase_eta0_mle = mean(se_eta0_mle, 
                                na.rm = TRUE), 
            bias_eta1_mle = mean(eta1_mle - eta1), 
            ese_eta1_mle = sd(eta1_mle), 
            ase_eta1_mle = mean(se_eta1_mle, 
                                na.rm = TRUE)
            ) |> 
  write.csv("~/Documents/food_access_misclassification/Simulations/sims_try_mlePossum2_with_Z_summary.csv", 
            row.names = FALSE)

# Define color palette 
slide_colors = c("#F2958D", "#ADA264", "#EB5F3F", "#3A5724", "#FAB825", "#E67B33") 

# Create boxplot of estimated prevalence ratio
est_box = res |> 
  select(sim, N, beta0, beta1, beta2, ppv, 
                starts_with(c("beta0_", "beta1_", "beta2_"))) |> 
  tidyr::gather(key = "Coeff_Method", value = "Estimate", -c(1:6)) |> 
  mutate(
    Method = case_when(
      grepl(x = Coeff_Method, pattern ="gs") ~ "Gold\nStandard",
      grepl(x = Coeff_Method, pattern ="mle") ~ "MLE",
      grepl(x = Coeff_Method, pattern ="n") ~ "Naive",
      grepl(x = Coeff_Method, pattern ="cc") ~ "Complete\nCase"
    ), 
    Method = factor(x = Method, 
                    levels = c("Gold\nStandard", "MLE", "Naive", "Complete\nCase", "Truth")), 
    Truth = case_when(
      grepl(x = Coeff_Method, pattern ="beta0") ~ beta0,
      grepl(x = Coeff_Method, pattern ="beta1") ~ beta1,
      grepl(x = Coeff_Method, pattern ="beta2") ~ beta2
    ), 
    Coeff = case_when(
      grepl(x = Coeff_Method, pattern ="beta0") ~ "Intercept",
      grepl(x = Coeff_Method, pattern ="beta1") ~ "Log PR of X | Z",
      grepl(x = Coeff_Method, pattern ="beta2") ~ "Log PR of Z | X"
    ), 
    Coeff = factor(x = Coeff, 
                   levels = c("Intercept", "Log PR of X | Z", "Log PR of Z | X")),
    N = factor(x = N, 
               levels = c(390, 2200), 
               labels = c("N = 390", "N = 2200"))) |> 
  ggplot(aes(x = factor(ppv), y = exp(Estimate), fill = Method)) + 
  geom_hline(aes(yintercept = exp(Truth)), 
             linetype = 2, 
             color = "darkgrey") + 
  geom_boxplot() + 
  scale_fill_manual(values = slide_colors[c(-2)], 
                    name = "Method:") + 
  facet_wrap(N ~ Coeff, 
             scales = "free") + 
  theme_minimal(base_size = 20) + 
  xlab(label = latex2exp::TeX(input = "Positive Predictive Value of X*")) + 
  ylab(label = latex2exp::TeX(input = "Estimated Prevalence [Ratio]")) + 
  theme(legend.position = "right", 
        legend.box.margin = margin(t = 10, r = 5, b = -30, l = 5), 
        plot.margin = margin(t = 0, r = 10, b = 0, l = 10), 
        strip.background = element_rect(fill = slide_colors[2])) 

# Save it 
ggsave(filename = "~/Documents/food_access_misclassification/Figures/sim_box_pr_try_mlePossum2_with_Z.png",
       plot = est_box, 
       device = "png", 
       width = 7, 
       height = 5)
