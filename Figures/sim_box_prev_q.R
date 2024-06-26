# Load packages
library(ggplot2)

# Read simulation data in from GitHub repo
#url_stem = "https://raw.githubusercontent.com/ashleymullan/food_access_misclassification/main/Simulations/vary_q/"
url_stem = "~/Documents/food_access_misclassification/Simulations/vary_q/"
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
slide_colors = c("#F2958D", "#ADA264", "#EB5F3F", "#3A5724", "#FAB825", "#E67B33") 

# Create boxplot of estimated prevalence ratio
est_box = res |> 
  dplyr::select(sim, N, beta0, q, dplyr::starts_with("beta0_")) |> 
  tidyr::gather(key = "Method", value = "est_beta0", -c(1:4)) |> 
  dplyr::mutate(Method = factor(x = Method, 
                                levels = c("beta0_gs", "beta0_smle", "beta0_cc", "beta0_n"), 
                                labels = c("Gold\nStandard", "SMLE", "Complete\nCase", "Naive")),
                N = factor(x = N, 
                           levels = c(390, 2200), 
                           labels = c("N = 390", "N = 2200"))) |> 
  ggplot(aes(x = factor(q), y = exp(est_beta0), fill = Method)) + 
  geom_hline(aes(yintercept = exp(beta0)), 
             linetype = 2, 
             color = "darkgrey") + 
  geom_boxplot() + 
  scale_fill_manual(values = slide_colors[c(-2)], 
                    name = "Method:", 
                    guide = "none") + 
  facet_wrap(~N) + 
  theme_minimal(base_size = 20) + 
  xlab(label = latex2exp::TeX(input = "Proportion of Neighborhoods Queried")) + 
  ylab(label = latex2exp::TeX(input = "Estimated Baseline Prevalence")) + 
  theme(legend.position = "right", 
        legend.box.margin = margin(t = 10, r = 5, b = -30, l = 5), 
        plot.margin = margin(t = 0, r = 10, b = 0, l = 10), 
        strip.background = element_rect(fill = slide_colors[2])) + 
  guides(fill = "none") #guide_legend(nrow=2, byrow=TRUE)) 

# Save it 
ggsave(filename = "~/Documents/food_access_misclassification/Figures/sim_box_prev_q.png",
       plot = est_box, 
       device = "png", 
       width = 7, 
       height = 5)

# Create boxplot of estimated prevalence ratio
est_box = res |> 
  dplyr::filter(N == 390) |> 
  dplyr::select(sim, N, beta0, q, dplyr::starts_with("beta0_")) |> 
  tidyr::gather(key = "Method", value = "est_beta0", -c(1:4)) |> 
  dplyr::mutate(Method = factor(x = Method, 
                                levels = c("beta0_gs", "beta0_smle", "beta0_cc", "beta0_n"), 
                                labels = c("Gold\nStandard", "SMLE", "Complete\nCase", "Naive")),
                N = factor(x = N, 
                           levels = c(390, 2200), 
                           labels = c("N = 390", "N = 2200"))) |> 
  ggplot(aes(x = factor(q), y = exp(est_beta0), fill = Method)) + 
  geom_hline(aes(yintercept = exp(beta0)), 
             linetype = 2, 
             color = "darkgrey") + 
  geom_boxplot() + 
  scale_fill_manual(values = slide_colors[c(-2)], 
                    name = "Method:", 
                    guide = "none") + 
  facet_wrap(~N) + 
  theme_minimal(base_size = 20) + 
  xlab(label = latex2exp::TeX(input = "Proportion of Neighborhoods Queried")) + 
  ylab(label = latex2exp::TeX(input = "Estimated Baseline Prevalence")) + 
  theme(legend.position = "right", 
        legend.box.margin = margin(t = 10, r = 5, b = -30, l = 5), 
        plot.margin = margin(t = 0, r = 10, b = 0, l = 10), 
        strip.background = element_rect(fill = slide_colors[2])) + 
  guides(fill = "none") #guide_legend(nrow=2, byrow=TRUE)) 

# Save it 
ggsave(filename = "~/Documents/food_access_misclassification/Figures/sim_box_prev_q_N390.png",
       plot = est_box, 
       device = "png", 
       width = 7, 
       height = 5)