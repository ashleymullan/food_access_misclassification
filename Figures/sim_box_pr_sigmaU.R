# Read simulation data in from GitHub repo
url_stem = "https://raw.githubusercontent.com/ashleymullan/food_access_misclassification/main/Simulations/vary_sigmaU/"
file_url = paste0(url_stem, 
                  "proximity_N", 
                  rep(c(390, 2200), each = 5), 
                  "_sigmaU", 
                  rep(c(10, 20, 40, 80, 100), times = 2), 
                  "_seed11422.csv")
res = do.call(what = rbind, 
              args = lapply(X = file_url, 
                            FUN = read.csv)
)

# Create summary of results
res |> 
  dplyr::group_by(N, sigmaU) |> 
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
                   bias_beta0_smle = mean(beta0_smle - beta0), 
                   ese_beta0_smle = sd(beta0_smle), 
                   ase_beta0_smle = mean(se_beta0_smle), 
                   bias_beta1_smle = mean(beta1_smle - beta1), 
                   ese_beta1_smle = sd(beta1_smle), 
                   ase_beta1_smle = mean(se_beta1_smle), 
                   cp_beta1_smle = mean((beta1_smle - 1.96 * se_beta1_smle) <= beta1 & 
                                          beta1 <= (beta1_smle + 1.96 * se_beta1_smle)), 
                   re_beta1_smle = var(beta1_gs) / var(beta1_smle)
  ) |> 
  write.csv("~/Documents/food_access_misclassification/Simulations/sims_vary_sigmaU_binX_summary.csv", 
            row.names = FALSE)

# Define color palette 
# slide_colors = c("#C3CFFA", "#C9B0B0", "#DD6D53", "#60CAD6", "#F0D290", "#fac3cf") 
slide_colors = c("#F2958D", "#ADA264", "#EB5F3F", "#3A5724", "#FAB825", "#E67B33") 

# Create boxplot of estimated 
est_box = res |> 
  dplyr::select(-dplyr::starts_with("se_"), 
                -dplyr::contains("beta0"), 
                -ppv_smle) |> 
  tidyr::gather(key = "Method", value = "est_beta1", -c(1:8)) |> 
  dplyr::mutate(Method = factor(x = Method, 
                                levels = c("beta1_gs", "beta1_smle", "beta1_cc", "beta1_n"), 
                                labels = c("Gold Standard", "SMLE", "Complete Case", "Naive")),
                N = factor(x = N, 
                           levels = c(390, 2200), 
                           labels = c("N = 390", "N = 2200"))) |> 
  ggplot(aes(x = factor(sigmaU), y = exp(est_beta1), fill = Method)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = exp(beta1)), linetype = 2, color = slide_colors[5]) + 
  scale_fill_manual(values = slide_colors, 
                    name = "Method:") + 
  facet_wrap(~N) + 
  theme_minimal(base_size = 40) + 
  xlab(label = latex2exp::TeX(input = "Additive Error Variance ($\\sigma_U$)")) + 
  ylab(label = latex2exp::TeX(input = "Estimated Prevalence Ratio ($\\exp(\\hat{\\beta}_1)$)")) + 
  theme(legend.position = "top", 
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(t = 10, r = 5, b = -40, l = 5), 
        plot.margin = margin(t = 0, r = 10, b = 0, l = 10))
# Save it 
ggsave(filename = "~/Documents/food_access_misclassification/Figures/sim_box_pr_sigmaU.png",
       plot = est_box, 
       device = "png", 
       width = 8, 
       height = 5)
