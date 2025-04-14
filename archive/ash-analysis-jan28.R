
library(dplyr)
library(possum)

local_dir <- "/Users/ashleymullan/Documents/Food-Access/"
data_path <- "food_access_misclassification/Analysis/"

data_1 <- read.csv(paste0(local_dir,
                          data_path,
                          "piedmont_data.csv"))

data_5 <- read.csv(paste0(local_dir,
                          data_path,
                          "piedmont_data5.csv"))

data_1 <- data_1 |>
  mutate(logpop = log(POP),
         logLA = log(LandArea) - median(log(LandArea)),
         logscaledpop = log(POP / 1000)) #to change interpretation to rate per 10K people

data_5 <- data_5 |>
  mutate(logpop = log(POP),
         logLA = log(LandArea) - median(log(LandArea)),
         logscaledpop = log(POP / 1000)) #to change interpretation to rate per 10K people

#gold standard
g1 <- summary(glm(formula = DIABETES ~ binX_full + logLA,
                 family = poisson(link = "log"),
                 offset = logscaledpop,
                 data = data_1))$coefficients

g5 <- summary(glm(formula = DIABETES ~ binX_full + logLA,
                 family = poisson(link = "log"),
                 offset = logscaledpop,
                 data = data_5))$coefficients

#naive
n1 <- summary(glm(formula = DIABETES ~ binXstar + logLA,
                              family = poisson(link = "log"),
                              offset = logscaledpop,
                              data = data_1))$coefficients

n5 <- summary(glm(formula = DIABETES ~ binXstar + logLA,
                 family = poisson(link = "log"),
                 offset = logscaledpop,
                 data = data_5))$coefficients

#mle.
m1 <- mlePossum(analysis_formula = DIABETES ~ binX_partial + logLA + offset(logscaledpop),
          family = poisson,
          error_formula = binX_partial ~ binXstar + logLA,
          data = data_1,
          beta_init = "Complete-data",
          eta_init = "Complete-data",
          noSE = FALSE,
          alternative_SE = FALSE,
          hN_scale = 0.5)

m5 <- mlePossum(analysis_formula = DIABETES ~ binX_partial + logLA + offset(logscaledpop),
                family = poisson,
                error_formula = binX_partial ~ binXstar + logLA,
                data = data_1,
                beta_init = "Complete-data",
                eta_init = "Complete-data",
                noSE = FALSE,
                alternative_SE = FALSE,
                hN_scale = 0.5)

#cc
c1 <- summary(glm(formula = DIABETES ~ binX_partial + logLA,
            family = poisson(link = "log"),
            offset = logscaledpop,
            data = data_1))$coefficients

c5 <- summary(glm(formula = DIABETES ~ binX_partial + logLA,
                 family = poisson(link = "log"),
                 offset = logscaledpop,
                 data = data_5))$coefficients

names(m1$coefficients) <- c("Estimate", "Std. Error")
names(m5$coefficients) <- c("Estimate", "Std. Error")

method <- c(rep("g", times = 3),
            rep("n", times = 3),
            rep("m", times = 3),
            rep("c", times = 3))

results1 <- rbind(g1[,1:2],n1[,1:2],m1$coefficients,c1[,1:2]) |> cbind(method)
results5 <- rbind(g5[,1:2],n5[,1:2],m5$coefficients,c5[,1:2]) |> cbind(method)
