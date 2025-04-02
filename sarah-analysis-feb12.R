
library(dplyr)
library(possum)

#local_dir <- "/Users/ashleymullan/Documents/Food-Access/"
local_dir <- "~/Documents/"
data_path <- "food_access_misclassification/Analysis/"

data = read.csv("https://raw.githubusercontent.com/sarahlotspeich/food_access_imputation/refs/heads/main/piedmont-triad-data/analysis_data.csv") |> 
  dplyr::mutate(X1mile = as.numeric(X_full <= 1), 
                Xstar1mile = as.numeric(Xstar <= 1))
set.seed(212)
sample_ct = sample(x = data$GEOID[data$Xstar1mile == 1], size = 50)
data$X1mile_miss = data$X1mile
data$X1mile_miss[data$Xstar1mile == 1 & !(data$GEOID %in% sample_ct)] = NA

#gold standard
g1 <- summary(glm(formula = Y_DIABETES ~ X1mile * METRO + offset(log(O_POP)),
                 family = poisson(link = "log"),
                 data = data))$coefficients

#naive
n1 <- summary(glm(formula = Y_DIABETES ~ Xstar1mile * METRO + offset(log(O_POP)),
                  family = poisson(link = "log"),
                  data = data))$coefficients

#mle.
m1 <- mlePossum(analysis_formula = Y_DIABETES ~ X1mile_miss * METRO + offset(log(O_POP)),
                family = poisson,
                error_formula = X1mile_miss ~ Xstar1mile * METRO,
                data = data,
                beta_init = "Complete-data",
                eta_init = "Complete-data",
                noSE = FALSE,
                alternative_SE = FALSE,
                hN_scale = 0.5)

m1_nd <- mlePossum(analysis_formula = Y_DIABETES ~ X1mile_miss * METRO + offset(log(O_POP)),
                   family = poisson,
                   error_formula = X1mile_miss ~ Xstar1mile * METRO,
                   data = data,
                   beta_init = "Complete-data",
                   eta_init = "Complete-data",
                   noSE = FALSE,
                   alternative_SE = TRUE)

#cc
c1 <- summary(glm(formula = Y_DIABETES ~ X1mile_miss * METRO + offset(log(O_POP)),
              family = poisson(link = "log"),
              offset = logscaledpop,
              data = data_1))$coefficients

names(m1$coefficients) <- c("Estimate", "Std. Error")
names(m5$coefficients) <- c("Estimate", "Std. Error")

method <- c(rep("g", times = 3),
            rep("n", times = 3),
            rep("m", times = 3),
            rep("c", times = 3))

results1 <- rbind(g1[,1:2],n1[,1:2],m1$coefficients,c1[,1:2]) |> cbind(method)
results5 <- rbind(g5[,1:2],n5[,1:2],m5$coefficients,c5[,1:2]) |> cbind(method)
