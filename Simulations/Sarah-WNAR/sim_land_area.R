#piedmont_data = read.csv("https://raw.githubusercontent.com/ashleymullan/food_access_misclassification/main/Analysis/piedmont_data.csv")
piedmont_data = read.csv("~/Documents/food_access_misclassification/Analysis/piedmont_data.csv")

## Unadjusted model (gold standard)
mod_unadj = glm(formula = Y_DIABETES ~ binX_full + offset(log(O_POP)), 
                family = quasipoisson,
                data = piedmont_data)
# mod_unadj$coefficients
# (Intercept)   binX_full 
# -2.2456410   0.1549768 

## Adjusted model (gold standard)
mod_adj = glm(formula = Y_DIABETES ~ binX_full + I(LandArea/100) + offset(log(O_POP)), 
              family = quasipoisson,
              data = piedmont_data)
# mod_adj$coefficients
# (Intercept)      binX_full I(LandArea/100)
# -2.27785405     0.17645438     0.01413503

## Misclassification mechanism
mod_mc = glm(formula = binX_full ~ I(LandArea/10), 
             family = binomial,
             data = piedmont_data, 
             subset = binXstar == 1)
# mod_mc$coefficients
# (Intercept) I(LandArea/10) 
# 0.3924274      0.3897812

## Fit Gamma distribution to land area
### With X* = 0 --> Gamma(1, 0.26)
EnvStats::egamma(piedmont_data$LandArea[piedmont_data$binXstar == 0] / 100, method = "mle", ci = FALSE, 
                 ci.type = "two-sided", ci.method = "normal.approx", 
                 normal.approx.transform = "kulkarni.powar", conf.level = 0.95)
### With X* = 1 --> Gamma(1, 0.06)
EnvStats::egamma(piedmont_data$LandArea[piedmont_data$binXstar == 1] / 100, method = "mle", ci = FALSE, 
                 ci.type = "two-sided", ci.method = "normal.approx", 
                 normal.approx.transform = "kulkarni.powar", conf.level = 0.95)

## Compare fitted Gamma distribution to empirical CDF
library(ggplot2)
piedmont_data |> 
  ggplot(aes(x = LandArea / 100, 
             col = factor(binXstar), 
             group = factor(binXstar))) + 
  stat_ecdf() + 
  scale_color_manual(values = c("violetred1", 
                                "slateblue2")) + 
  stat_function(fun = pgamma, 
                args = c(shape = 1, scale = 0.26), 
                linetype = 2, 
                color = "violetred1")  + 
  stat_function(fun = pgamma, 
                args = c(shape = 1, scale = 0.06), 
                linetype = 2, 
                color = "slateblue2")