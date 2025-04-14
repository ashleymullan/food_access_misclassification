# Load packages
library(ggplot2)

# Read simulation data in from GitHub repo
url_stem = "https://raw.githubusercontent.com/ashleymullan/food_access_misclassification/main/Simulations/try_mlePossum2/"
file_url = paste0(url_stem, 
                  "proximity_N", 
                  rep(c(390, 2200), each = 4), 
                  "_ppv", 
                  rep(c(50, 60, 80, 90), times = 2), 
                  "_seed11422.csv")
res = do.call(what = rbind, 
              args = lapply(X = file_url, 
                            FUN = read.csv)
)

# Create summary of results
res |> 
  dplyr::mutate(eta0 = log(ppv / (1 - ppv))) |> 
  dplyr::group_by(N, ppv) |> 
  dplyr::summarize(avg_fpr = mean(fpr), 
                   bias_beta0_gs = mean(beta0_gs - beta0), 
                   ese_beta0_gs = sd(beta0_gs), 
                   bias_beta1_gs = mean(beta1_gs - beta1), 
                   ese_beta1_gs = sd(beta1_gs), 
                   bias_beta0_n = mean(beta0_n - beta0), 
                   ese_beta0_n = sd(beta0_n), 
                   bias_beta1_n = mean(beta1_n - beta1), 
                   ese_beta1_n = sd(beta1_n), 
                   bias_beta0_cc = mean(beta0_cc - beta0), 
                   ese_beta0_cc = sd(beta0_cc), 
                   bias_beta1_cc = mean(beta1_cc - beta1), 
                   ese_beta1_cc = sd(beta1_cc), 
                   re_beta1_cc = var(beta1_gs) / var(beta1_cc),
                   bias_beta0_mle = mean(beta0_mle - beta0), 
                   ese_beta0_mle = sd(beta0_mle), 
                   ase_beta0_mle = mean(se_beta0_mle), 
                   bias_beta1_mle = mean(beta1_mle - beta1), 
                   ese_beta1_mle = sd(beta1_mle), 
                   ase_beta1_mle = mean(se_beta1_mle), 
                   cp_beta1_mle = mean((beta1_mle - 1.96 * se_beta1_mle) <= beta1 & 
                                          beta1 <= (beta1_mle + 1.96 * se_beta1_mle)), 
                   re_beta1_mle = var(beta1_gs) / var(beta1_mle),
                   bias_eta0_mle = mean(eta0_mle - eta0), 
                   ese_eta0_mle = sd(eta0_mle), 
                   ase_eta0_mle = mean(se_eta0_mle), 
                   ) |> 
  write.csv("~/Documents/food_access_misclassification/Simulations/sims_try_mlePossum2_summary.csv", 
            row.names = FALSE)

# Define color palette 
slide_colors = c("#F2958D", "#ADA264", "#EB5F3F", "#3A5724", "#FAB825", "#E67B33") 

# Create boxplot of estimated prevalence ratio
est_box = res |> 
  dplyr::select(sim, N, beta1, ppv, dplyr::starts_with("beta1_")) |> 
  tidyr::gather(key = "Method", value = "est_beta1", -c(1:4)) |> 
  dplyr::mutate(Method = factor(x = Method, 
                                levels = c("beta1_gs", "beta1_mle", "beta1_cc", "beta1_n"), 
                                labels = c("Gold\nStandard", "MLE", "Complete\nCase", "Naive")),
                N = factor(x = N, 
                           levels = c(390, 2200), 
                           labels = c("N = 390", "N = 2200"))) |> 
  ggplot(aes(x = factor(ppv), y = exp(est_beta1), fill = Method)) + 
  geom_hline(aes(yintercept = exp(beta1)), 
             linetype = 2, 
             color = "darkgrey") + 
  geom_boxplot() + 
  scale_fill_manual(values = slide_colors[c(-2)], 
                    name = "Method:") + 
  facet_wrap(~N) + 
  theme_minimal(base_size = 20) + 
  xlab(label = latex2exp::TeX(input = "Positive Predictive Value of X*")) + 
  ylab(label = latex2exp::TeX(input = "Estimated Prevalence Ratio")) + 
  theme(legend.position = "right", 
        legend.box.margin = margin(t = 10, r = 5, b = -30, l = 5), 
        plot.margin = margin(t = 0, r = 10, b = 0, l = 10), 
        strip.background = element_rect(fill = slide_colors[2])) 

# Save it 
ggsave(filename = "~/Documents/food_access_misclassification/Figures/sim_box_pr_try_mlePossum2.png",
       plot = est_box, 
       device = "png", 
       width = 7, 
       height = 5)
