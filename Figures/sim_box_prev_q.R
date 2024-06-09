# Read simulation data in from GitHub repo
url_stem = "https://raw.githubusercontent.com/ashleymullan/food_access_misclassification/main/Simulations/vary_q/"
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

# Define color palette 
# slide_colors = c("#C3CFFA", "#C9B0B0", "#DD6D53", "#60CAD6", "#F0D290", "#fac3cf") 
slide_colors = c("#F2958D", "#ADA264", "#EB5F3F", "#3A5724", "#FAB825", "#E67B33") 

# Create boxplot of estimated 
est_box = res |> 
  dplyr::select(-dplyr::starts_with("se_"), 
                -dplyr::contains("beta1")) |> 
  tidyr::gather(key = "Method", value = "est_beta0", -c(1:8)) |> 
  dplyr::mutate(Method = factor(x = Method, 
                                levels = c("beta0_gs", "beta0_smle", "beta0_cc", "beta0_n"), 
                                labels = c("Gold Standard", "SMLE", "Complete Case", "Naive")),
                N = factor(x = N, 
                           levels = c(390, 2200), 
                           labels = c("N = 390", "N = 2200"))) |> 
  ggplot(aes(x = factor(q), y = exp(est_beta0), fill = Method)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = exp(beta0)), linetype = 2, color = slide_colors[5]) + 
  scale_fill_manual(values = slide_colors, 
                    name = "Method:") + 
  facet_wrap(~N) + 
  theme_minimal(base_size = 16) + 
  scale_y_continuous(labels = scales::percent) + 
  xlab(label = latex2exp::TeX(input = "Proportion of Neighborhoods Queried ($q$)")) + 
  ylab(label = latex2exp::TeX(input = "Estimated Baseline Prevalence ($\\exp(\\hat{\\beta}_0)$)")) + 
  theme(legend.position = "top", 
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(t = 10, r = 5, b = 10, l = -10), 
        plot.margin = margin(t = 0, r = 10, b = 0, l = -10))
# Save it 
ggsave(filename = "~/Documents/food_access_misclassification/Figures/sim_box_prev_q.png",
       plot = est_box, 
       device = "png", 
       width = 8, 
       height = 5, 
       units = "in")
