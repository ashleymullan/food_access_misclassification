# Load packages
library(ggplot2)

# Read simulation data in from GitHub repo
url_stem = "https://raw.githubusercontent.com/ashleymullan/food_access_misclassification/main/Simulations/vary_design/"
file_url = paste0(url_stem, 
                  "proximity_N", 
                  rep(c(390, 2200), each = 4), 
                  "_q", 
                  rep(c(10, 25, 50, 75), times = 2), 
                  "_seed11422.csv")
res = do.call(what = rbind, 
              args = lapply(X = file_url, 
                            FUN = read.csv)
)

# Create summary of results
res |> 
  dplyr::group_by(N, q) |> 
  dplyr::summarize(avg_fpr = mean(fpr), 
                   bias_beta0_srs = mean(beta0_srs - beta0), 
                   ese_beta0_srs = sd(beta0_srs), 
                   bias_beta1_srs = mean(beta1_srs - beta1), 
                   ese_beta1_srs = sd(beta1_srs), 
                   ase_beta1_srs = mean(se_beta1_srs), 
                   cp_beta1_srs = mean((beta1_srs - 1.96 * se_beta1_srs) <= beta1 & 
                                          beta1 <= (beta1_srs + 1.96 * se_beta1_srs)), 
                   re_beta1_srs = var(beta1_srs) / var(beta1_srs),
                   bias_beta0_bcc = mean(beta0_bcc - beta0), 
                   ese_beta0_bcc = sd(beta0_bcc), 
                   bias_beta1_bcc = mean(beta1_bcc - beta1), 
                   ese_beta1_bcc = sd(beta1_bcc), 
                   ase_beta1_bcc = mean(se_beta1_bcc), 
                   cp_beta1_bcc = mean((beta1_bcc - 1.96 * se_beta1_bcc) <= beta1 & 
                                         beta1 <= (beta1_bcc + 1.96 * se_beta1_bcc)), 
                   re_beta1_bcc = var(beta1_srs) / var(beta1_bcc),
                   bias_beta0_cc = mean(beta0_cc - beta0), 
                   ese_beta0_cc = sd(beta0_cc), 
                   bias_beta1_cc = mean(beta1_cc - beta1), 
                   ese_beta1_cc = sd(beta1_cc), 
                   ase_beta1_cc = mean(se_beta1_cc), 
                   cp_beta1_cc = mean((beta1_cc - 1.96 * se_beta1_cc) <= beta1 & 
                                         beta1 <= (beta1_cc + 1.96 * se_beta1_cc)), 
                   re_beta1_cc = var(beta1_srs) / var(beta1_cc),
                   bias_beta0_resid = mean(beta0_resid - beta0), 
                   ese_beta0_resid = sd(beta0_resid), 
                   ase_beta0_resid = mean(se_beta0_resid), 
                   bias_beta1_resid = mean(beta1_resid - beta1), 
                   ese_beta1_resid = sd(beta1_resid), 
                   ase_beta1_resid = mean(se_beta1_resid), 
                   cp_beta1_resid = mean((beta1_resid - 1.96 * se_beta1_resid) <= beta1 & 
                                          beta1 <= (beta1_resid + 1.96 * se_beta1_resid)), 
                   re_beta1_resid = var(beta1_srs) / var(beta1_resid)
                   ) |> 
  write.csv("~/Documents/food_access_misclassification/Simulations/sims_vary_design_binX_summary.csv", 
            row.names = FALSE)

# Define color palette 
slide_colors = c("#F2958D", "#ADA264", "#EB5F3F", "#3A5724", "#FAB825", "#E67B33") 

# Create boxplot of estimated prevalence ratio
est_box = res |> 
  dplyr::filter(N == 390) |> 
  dplyr::select(-dplyr::starts_with("se_"), 
                -dplyr::contains("beta0"), 
                -ppv) |> 
  tidyr::gather(key = "Design", value = "est_beta1", -c(1:8)) |> 
  dplyr::mutate(Design = factor(x = Design, 
                                levels = c("beta1_srs", "beta1_cc", "beta1_bcc", "beta1_resid"), 
                                labels = c("SRS", "50/50", "Stratified 50/50", "Residual")),
                N = factor(x = N, 
                           levels = c(390, 2200), 
                           labels = c("N = 390", "N = 2200"))) |> 
  ggplot(aes(x = factor(q), y = exp(est_beta1), fill = Design)) + 
  geom_hline(aes(yintercept = exp(beta1)), 
             linetype = 2, 
             color = "darkgrey") + 
  geom_boxplot() + 
  scale_fill_manual(values = slide_colors[c(-2)], 
                    name = "Method:") + 
  facet_wrap(~N) + 
  theme_minimal(base_size = 20) + 
  xlab(label = latex2exp::TeX(input = "Proportion of Neighborhoods Queried ($q$)")) + 
  ylab(label = latex2exp::TeX(input = "Estimated Prevalence Ratio ($\\exp(\\hat{\\beta}_1)$)")) + 
  theme(legend.position = "right", 
        legend.box.margin = margin(t = 10, r = 5, b = -30, l = 5), 
        plot.margin = margin(t = 0, r = 10, b = 0, l = 10), 
        strip.background = element_rect(fill = slide_colors[2])) + 
  guides(fill = guide_legend(nrow=2, byrow=TRUE)) 

# Save it 
ggsave(filename = "~/Documents/food_access_misclassification/Figures/sim_box_pr_design_N390.png",
       plot = est_box, 
       device = "png", 
       width = 10, 
       height = 5)

est_box = res |> 
  dplyr::filter(N == 2200) |> 
  dplyr::select(-dplyr::starts_with("se_"), 
                -dplyr::contains("beta0"), 
                -ppv) |> 
  tidyr::gather(key = "Design", value = "est_beta1", -c(1:8)) |> 
  dplyr::mutate(Design = factor(x = Design, 
                                levels = c("beta1_srs", "beta1_cc", "beta1_bcc", "beta1_resid"), 
                                labels = c("SRS", "50/50", "Stratified 50/50", "Residual")),
                N = factor(x = N, 
                           levels = c(390, 2200), 
                           labels = c("N = 390", "N = 2200"))) |> 
  ggplot(aes(x = factor(q), y = exp(est_beta1), fill = Design)) + 
  geom_hline(aes(yintercept = exp(beta1)), 
             linetype = 2, 
             color = "darkgrey") + 
  geom_boxplot() + 
  scale_fill_manual(values = slide_colors[c(-2)], 
                    name = "Method:") + 
  facet_wrap(~N) + 
  theme_minimal(base_size = 20) + 
  xlab(label = latex2exp::TeX(input = "Proportion of Neighborhoods Queried ($q$)")) + 
  ylab(label = latex2exp::TeX(input = "Estimated Prevalence Ratio ($\\exp(\\hat{\\beta}_1)$)")) + 
  theme(legend.position = "right", 
        legend.box.margin = margin(t = 10, r = 5, b = -30, l = 5), 
        plot.margin = margin(t = 0, r = 10, b = 0, l = 10), 
        strip.background = element_rect(fill = slide_colors[2])) + 
  guides(fill = guide_legend(nrow=2, byrow=TRUE)) 

# Save it 
ggsave(filename = "~/Documents/food_access_misclassification/Figures/sim_box_pr_design_N2200.png",
       plot = est_box, 
       device = "png", 
       width = 10, 
       height = 5)
