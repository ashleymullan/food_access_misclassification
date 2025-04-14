# Load packages
library(ggplot2)

# Read simulation data in from GitHub repo
url_stem = "https://raw.githubusercontent.com/ashleymullan/food_access_misclassification/main/Simulations/vary_access/"
file_url = paste0(url_stem, 
                  "proximity_N", 
                  rep(c(390, 2200), each = 4), 
                  "_", 
                  rep(c(0.5, 1, 5, 10), times = 2), 
                  "mile_seed11422.csv")
res = do.call(what = rbind, 
              args = lapply(X = file_url, 
                            FUN = read.csv)
)

# Create summary of results
res |> 
  dplyr::group_by(N, access) |> 
  dplyr::summarize(avg_fpr = mean(fpr, na.rm = TRUE),
                   bias_beta0_gs = mean(beta0_gs - beta0, na.rm = TRUE),
                   ese_beta0_gs = sd(beta0_gs, na.rm = TRUE),
                   bias_beta1_gs = mean(beta1_gs - beta1, na.rm = TRUE),
                   ese_beta1_gs = sd(beta1_gs, na.rm = TRUE),
                   bias_beta0_n = mean(beta0_n - beta0, na.rm = TRUE),
                   ese_beta0_n = sd(beta0_n, na.rm = TRUE),
                   bias_beta1_n = mean(beta1_n - beta1, na.rm = TRUE),
                   ese_beta1_n = sd(beta1_n, na.rm = TRUE),
                   bias_beta0_cc = mean(beta0_cc - beta0, na.rm = TRUE),
                   ese_beta0_cc = sd(beta0_cc, na.rm = TRUE),
                   bias_beta1_cc = mean(beta1_cc - beta1, na.rm = TRUE),
                   ese_beta1_cc = sd(beta1_cc, na.rm = TRUE),
                   re_beta1_cc = var(beta1_gs, na.rm = TRUE) / var(beta1_cc, na.rm = TRUE),
                   bias_beta0_smle = mean(beta0_smle - beta0, na.rm = TRUE),
                   ese_beta0_smle = sd(beta0_smle, na.rm = TRUE),
                   ase_beta0_smle = mean(se_beta0_smle, na.rm = TRUE),
                   bias_beta1_smle = mean(beta1_smle - beta1, na.rm = TRUE),
                   ese_beta1_smle = sd(beta1_smle, na.rm = TRUE),
                   ase_beta1_smle = mean(se_beta1_smle, na.rm = TRUE),
                   cp_beta1_smle = mean((beta1_smle - 1.96 * se_beta1_smle) <= beta1 & 
                                          beta1 <= (beta1_smle + 1.96 * se_beta1_smle), na.rm = TRUE),
                   re_beta1_smle = var(beta1_gs, na.rm = TRUE) / var(beta1_smle, na.rm = TRUE)
                   ) |> 
  write.csv("~/Documents/food_access_misclassification/Simulations/sims_vary_access_binX_summary.csv", 
            row.names = FALSE)

# Define color palette 
slide_colors = c("#F2958D", "#ADA264", "#EB5F3F", "#3A5724", "#FAB825", "#E67B33") 

# Create boxplot of estimated prevalence ratio
est_box = res |> 
  dplyr::filter(N == 390) |> 
  dplyr::select(-dplyr::starts_with("se_"), 
                -dplyr::contains("beta0"), 
                -ppv_smle, 
                -msg) |> 
  tidyr::gather(key = "Method", value = "est_beta1", -c(1:9)) |> 
  dplyr::mutate(Method = factor(x = Method, 
                                levels = c("beta1_gs", "beta1_smle", "beta1_cc", "beta1_n"), 
                                labels = c("Gold\nStandard", "SMLE", "Complete\nCase", "Naive")),
                N = factor(x = N, 
                           levels = c(390, 2200), 
                           labels = c("N = 390", "N = 2200")), 
                access = factor(x = access, 
                                levels = c(0.5, 1, 5, 10), 
                                labels = c("Within 0.5 Miles", 
                                           "Within 1 Mile", 
                                           "Within 5 Miles",
                                           "Within 10 Miles"))) |> 
  ggplot(aes(x = access, y = exp(est_beta1), fill = Method)) + 
  geom_hline(aes(yintercept = exp(beta1)), 
             linetype = 2, 
             color = "darkgrey") + 
  geom_boxplot() + 
  scale_fill_manual(values = slide_colors[c(-2)], 
                    name = "Method:") + 
  facet_wrap(~N) + 
  theme_minimal(base_size = 20) + 
  xlab(label = latex2exp::TeX(input = "Threshold for Food Access")) + 
  ylab(label = latex2exp::TeX(input = "Estimated Prevalence Ratio")) + 
  theme(legend.position = "right", 
        legend.box.margin = margin(t = 10, r = 5, b = -30, l = 5), 
        plot.margin = margin(t = 0, r = 10, b = 0, l = 10), 
        strip.background = element_rect(fill = slide_colors[2])) + 
  guides(fill = guide_legend(nrow=2, byrow=TRUE)) 

# Save it 
ggsave(filename = "~/Documents/food_access_misclassification/Figures/sim_box_pr_access_N390.png",
       plot = est_box, 
       device = "png", 
       width = 10, 
       height = 5)

est_box = res |> 
  dplyr::filter(N == 2200) |> 
  dplyr::select(-dplyr::starts_with("se_"), 
                -dplyr::contains("beta0"), 
                -ppv_smle, 
                -msg) |> 
  tidyr::gather(key = "Method", value = "est_beta1", -c(1:9)) |> 
  dplyr::mutate(Method = factor(x = Method, 
                                levels = c("beta1_gs", "beta1_smle", "beta1_cc", "beta1_n"), 
                                labels = c("Gold\nStandard", "SMLE", "Complete\nCase", "Naive")),
                N = factor(x = N, 
                           levels = c(390, 2200), 
                           labels = c("N = 390", "N = 2200")), 
                access = factor(x = access, 
                                levels = c(0.5, 1, 5, 10), 
                                labels = c("Within 0.5 Miles", 
                                           "Within 1 Mile", 
                                           "Within 5 Miles",
                                           "Within 10 Miles"))) |> 
  ggplot(aes(x = access, y = exp(est_beta1), fill = Method)) + 
  geom_hline(aes(yintercept = exp(beta1)), 
             linetype = 2, 
             color = "darkgrey") + 
  geom_boxplot() + 
  scale_fill_manual(values = slide_colors[c(-2)], 
                    name = "Method:") + 
  facet_wrap(~N) + 
  theme_minimal(base_size = 20) + 
  xlab(label = latex2exp::TeX(input = "Threshold for Food Access")) + 
  ylab(label = latex2exp::TeX(input = "Estimated Prevalence Ratio")) + 
  theme(legend.position = "right", 
        legend.box.margin = margin(t = 10, r = 5, b = -30, l = 5), 
        plot.margin = margin(t = 0, r = 10, b = 0, l = 10), 
        strip.background = element_rect(fill = slide_colors[2])) + 
  guides(fill = guide_legend(nrow=2, byrow=TRUE)) 

# Save it 
ggsave(filename = "~/Documents/food_access_misclassification/Figures/sim_box_pr_access_N2200.png",
       plot = est_box, 
       device = "png", 
       width = 10, 
       height = 5)

# Plot FPR / PPV
res |> 
  dplyr::mutate(N = factor(x = N, 
                           levels = c(390, 2200), 
                           labels = c("N = 390", "N = 2200")), 
                access = factor(x = access, 
                                levels = c(0.5, 1, 5, 10), 
                                labels = c("Within 0.5 Miles", 
                                           "Within 1 Mile", 
                                           "Within 5 Miles",
                                           "Within 10 Miles"))) |> 
  dplyr::group_by(N, access) |> 
  dplyr::summarize(avg_fpr = mean(fpr, na.rm = TRUE), 
                   avg_ppv = mean(ppv_smle, na.rm = TRUE)) |> 
  ggplot(aes(x = access, linetype = N, group = N)) + 
  geom_line(aes(y = avg_fpr), color = slide_colors[1]) + 
  geom_line(aes(y = avg_ppv), color = slide_colors[2]) + 
  scale_fill_manual(values = slide_colors[c(-2)], 
                    name = "Method:") + 
  theme_minimal(base_size = 20) + 
  xlab(label = latex2exp::TeX(input = "Threshold for Food Access")) + 
  ylab(label = latex2exp::TeX(input = "Estimated Prevalence Ratio")) + 
  theme(legend.position = "right", 
        legend.box.margin = margin(t = 10, r = 5, b = -30, l = 5), 
        plot.margin = margin(t = 0, r = 10, b = 0, l = 10), 
        strip.background = element_rect(fill = slide_colors[2])) + 
  guides(fill = guide_legend(nrow=2, byrow=TRUE)) 
