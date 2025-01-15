
sims <- read.csv("one_sided_vary_ppv.csv")

library(latex2exp)





sims_long <- sims |>
  pivot_longer(cols = contains("beta1_"),
               names_to = "method_type",
               values_to = "betahat1")

simplot <- sims_long |>
  filter(method_type %in% c("beta1_cc", "beta1_gs", "beta1_mle",
                            "beta1_mle_strat", "beta1_n")) |> #cut to coefs
  filter(ppv %in% c(0.5, 0.9)) |> #extract two extremes
  mutate(Nstr = ifelse(N == 390, "N = 390", "N = 2200"),
         PPVstr = ifelse(ppv == 0.5, "PPV = 0.5", "PPV = 0.9"),
         grp = paste0(Nstr, ", ", PPVstr),
         method_type = case_when(
           method_type == "beta1_cc" ~ "Complete Case",
           method_type == "beta1_gs" ~ "Gold Standard",
           method_type == "beta1_mle" ~ "MLE",
           method_type == "beta1_mle_strat" ~ "MLE*",
           method_type == "beta1_n" ~ "Naive"
         )) |> #labels
  ggplot(aes(x = as.factor(method_type), y = exp(betahat1))) +
  geom_boxplot(aes(fill = method_type)) +
  theme_minimal() +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
  theme(plot.title = element_text(size = 9.5, hjust = 0.5),
        plot.subtitle = element_text(size = 9.5),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none",
        strip.background = element_rect(fill = "#ffd697"),
        strip.text = element_text(color = "black", face = "bold"),
        axis.text.x = element_text(size = 8)) +
  scale_fill_manual(values = c("#ed511f", "#ffb849", "#0a097d","#f4d1b1", "#8484be")) +
  facet_wrap(vars(grp)) +
  labs(y = TeX("$exp(\\hat{\\beta_1})$"),
       fill = "Method",
       x = "Method",
       title = "Method Comparison in Simulations")

ggsave("simplot.png", simplot, height = 4, width = 6, unit = "in", bg = "white")
