# Read simulation data in from GitHub repo
url_stem = "https://raw.githubusercontent.com/ashleymullan/food_access_misclassification/main/Simulations/vary_sigmaU/"
file_url = paste0("proximity_N", 
                  rep(c(390, 2200), each = 5), 
                  "_sigmaU", 
                  repc(c(10, 20, 40, 80, 100), times = 2), 
                  "_seed11422.csv")
res = do.call(what = rbind, 
              args = lapply(X = file_url), 
                            FUN = read.csv))

# Define color palette 
slide_colors = c("#C3CFFA", "#C9B0B0", "#DD6D53", "#60CAD6", "#F0D290", "#fac3cf") 

# Create line plot of average FPR per setting 
plot_fpr = res |> 
  dplyr::mutate(N = factor(x = N, 
                           levels = c(390, 2200), 
                           labels = c("N = 390", "N = 2200"))) |> 
  dplyr::group_by(N, sigmaU) |> 
  dplyr::summarize(avg_fpr = mean(fpr)) |> 
  ggplot(aes(x = sigmaU, y = avg_fpr, group = N, color = factor(N))) + 
  geom_point() + 
  geom_line() + 
  scale_color_manual(values = slide_colors, 
                     name = "Sample Size:") + 
  theme_minimal(base_size = 40) + 
  xlab(label = latex2exp::TeX(input = "Additive Error Variance ($\\sigma_U$)")) + 
  ylab(label = latex2exp::TeX(input = "Average False Positive Rate in $X^*$")) + 
  theme(legend.position = "top", 
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(t = 10, r = 5, b = -20, l = 5), 
        plot.margin = margin(t = 0, r = 10, b = 0, l = 10))

## Save it 
ggsave(filename = "~/Documents/food_access_misclassification/Figures/sim_fpr_sigmaU.png",
       plot = plot_fpr, 
       device = "png", 
       width = 8, 
       height = 5)

# Create boxplot of estimated 
est_box = res |> 
  dplyr::select(-dplyr::starts_with("se_"), 
                -dplyr::contains("beta0")) |> 
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
ggsave(filename = "~/Documents/food_access_misclassification/Figures/sim_box_sigmaU.png",
       plot = est_box, 
       device = "png", 
       width = 8, 
       height = 5)