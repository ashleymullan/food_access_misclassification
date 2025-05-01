
# load libraries
library(dplyr) ## for data wrangling
library(ggplot2) ## for graph
library(possum) ## for mle
library(ggimage) ## for emojis in graph
library(ggthemes) ## for colorblind palette
library(tidycensus) ## for shapefiles

local_dir <- "/Users/ashleymullan/Documents/Food-Access/"
plot_path <- "food_access_misclassification/Analysis/plots/"
emoji_path <- "food_access_misclassification/Analysis/"

# read in (and save only the necessary) data
data_file <- paste0("https://raw.githubusercontent.com/sarahlotspeich/food_access_imputation/",
                      "refs/heads/main/piedmont-triad-data/analysis_data.csv")
data <- read.csv(data_file) |>
  dplyr::select("GEOID", "CountyName", "Xstar", "X_full",
                "O_POP", "Y_DIABETES", "METRO")

# read in emoji icon paths for the plots
emoji_files <- c("city_image.png", "house_image.png")
emojis <- paste0(local_dir, emoji_path, emoji_files)
names(emojis) <- c("metro", "nonmetro")

#set up access indicators and recode covariates
data <- data |>
  ## is closest store within a half mile?
  mutate(Xh = ifelse(X_full <= 0.5, 1, 0), ## route based version
         Xstarh = ifelse(Xstar <= 0.5 ,1, 0)) |> ## Haversine version
  ## is closest store within a mile?
  mutate(X1 = ifelse(X_full <= 1, 1, 0), ## route based version
         Xstar1 = ifelse(Xstar <= 1, 1, 0)) |> ## Haversine version
  ## is closest store within 5 miles?
  mutate(X5 = ifelse(X_full <= 5, 1, 0), ## route based version
         Xstar5 = ifelse(Xstar <= 5, 1, 0)) |> ## Haversine version
  mutate(METRO = ifelse(METRO, 1, 0))

# set up missingness mechanism

## fix constants
set.seed(1031) ## for reproducibility
N <- nrow(data)
query_percent <- 0.20
n <- floor(query_percent * N)

## pick rows to query
### note: we have one sided misclassification, so we limit sampling to X* = 1
query_candsh <- which(data$Xstarh == 1)
query_rowsh <- sample(query_candsh, size = n, replace = FALSE)
queried_idsh <- data$GEOID[query_rowsh]
query_cands1 <- which(data$Xstar1 == 1)
query_rows1 <- sample(query_cands1, size = n, replace = FALSE)
queried_ids1 <- data$GEOID[query_rows1]
query_cands5 <- which(data$Xstar5 == 1)
query_rows5 <- sample(query_cands5, size = n, replace = FALSE)
queried_ids5 <- data$GEOID[query_rows5]

## remove nonqueried X at each radius
data <- data |>
  mutate(Xh_partial = ifelse(GEOID %in% queried_idsh, Xh, NA),
         X1_partial = ifelse(GEOID %in% queried_ids1, X1, NA),
         X5_partial = ifelse(GEOID %in% queried_ids5, X5, NA))

# fit models and save:
# coefficients, standard errors, p-values, and variance covariance matrices

## half mile radius
### gold standard
gh_full <- summary(glm(formula = Y_DIABETES ~ Xh * METRO + offset(log(O_POP)),
                       family = poisson(link = "log"),
                       data = data))
gh <- gh_full$coefficients
gh_vcov <- vcov(gh_full)

### naive
nh_full <- summary(glm(formula = Y_DIABETES ~ Xstarh * METRO + offset(log(O_POP)),
                       family = poisson(link = "log"),
                       data = data))
nh <- nh_full$coefficients
nh_vcov <- vcov(nh_full)

### mle
mh_full <- mlePossum(analysis_formula = Y_DIABETES ~ Xh_partial * METRO + offset(log(O_POP)),
                     family = poisson,
                     error_formula = Xh_partial ~ Xstarh * METRO,
                     data = data,
                     beta_init = "Complete-data",
                     eta_init = "Complete-data",
                     noSE = FALSE,
                     alternative_SE = FALSE,
                     hN_scale = 0.5)
mh <- mh_full$coefficients
mh_vcov <- mh_full$vcov

### complete case
ch_full <- summary(glm(formula = Y_DIABETES ~ Xh_partial * METRO + offset(log(O_POP)),
                       family = poisson(link = "log"),
                       data = data))
ch <- ch_full$coefficients
ch_vcov <- vcov(ch_full)

## one mile radius

### gold standard
g1_full <- summary(glm(formula = Y_DIABETES ~ X1 * METRO + offset(log(O_POP)),
                  family = poisson(link = "log"),
                  data = data))
g1 <- g1_full$coefficients
g1_vcov <- vcov(g1_full)

### naive
n1_full <- summary(glm(formula = Y_DIABETES ~ Xstar1 * METRO + offset(log(O_POP)),
                  family = poisson(link = "log"),
                  data = data))
n1 <- n1_full$coefficients
n1_vcov <- vcov(n1_full)

### mle
m1_full <- mlePossum(analysis_formula = Y_DIABETES ~ X1_partial * METRO + offset(log(O_POP)),
                family = poisson,
                error_formula = X1_partial ~ Xstar1 * METRO,
                data = data,
                beta_init = "Complete-data",
                eta_init = "Complete-data",
                noSE = FALSE,
                alternative_SE = FALSE,
                hN_scale = 0.5)
m1 <- m1_full$coefficients
m1_vcov <- m1_full$vcov

### complete case
c1_full <- summary(glm(formula = Y_DIABETES ~ X1_partial * METRO + offset(log(O_POP)),
                  family = poisson(link = "log"),
                  data = data))
c1 <- c1_full$coefficients
c1_vcov <- vcov(c1_full)

## five mile radius

### gold standard
g5_full <- summary(glm(formula = Y_DIABETES ~ X5 * METRO + offset(log(O_POP)),
                  family = poisson(link = "log"),
                  data = data))
g5 <- g5_full$coefficients
g5_vcov <- vcov(g5_full)

### naive
n5_full <- summary(glm(formula = Y_DIABETES ~ Xstar5 * METRO + offset(log(O_POP)),
                  family = poisson(link = "log"),
                  data = data))
n5 <- n5_full$coefficients
n5_vcov <- vcov(n5_full)

### mle
m5_full <- mlePossum(analysis_formula = Y_DIABETES ~ X5_partial * METRO + offset(log(O_POP)),
                family = poisson,
                error_formula = X5_partial ~ Xstar5 * METRO,
                data = data,
                beta_init = "Complete-data",
                eta_init = "Complete-data",
                noSE = FALSE,
                alternative_SE = FALSE,
                hN_scale = 0.1)
m5 <- m5_full$coefficients
m5_vcov <- m5_full$vcov

### complete case
c5_full <- summary(glm(formula = Y_DIABETES ~ X5_partial * METRO + offset(log(O_POP)),
                  family = poisson(link = "log"),
                  data = data))
c5 <- c5_full$coefficients
c5_vcov <- vcov(c5_full)

# set up data to plot

## one mile data

### save slopes in metro group
g1_mest <- g1["X1", "Estimate"] + g1["X1:METRO", "Estimate"]
n1_mest <- n1["Xstar1", "Estimate"] + n1["Xstar1:METRO", "Estimate"]
m1_mest <- m1["X1_partial", "Estimate"] + m1["X1_partial:METRO", "Estimate"]
c1_mest <- c1["X1_partial", "Estimate"] + c1["X1_partial:METRO", "Estimate"]

### save metro group standard error

#### gold standard
g1_vb1 <- g1_vcov["METRO", "METRO"]
g1_vb3 <- g1_vcov["X1:METRO", "X1:METRO"]
g1_cov <- g1_vcov["METRO", "X1:METRO"]
g1_mg_se <- sqrt(g1_vb1 + g1_vb3 + 2 * g1_cov)

#### naive
n1_vb1 <- n1_vcov["METRO", "METRO"]
n1_vb3 <- n1_vcov["Xstar1:METRO", "Xstar1:METRO"]
n1_cov <- n1_vcov["METRO", "Xstar1:METRO"]
n1_mg_se <- sqrt(n1_vb1 + n1_vb3 + 2 * n1_cov)

#### mle
m1_vb1 <- m1_vcov["METRO", "METRO"]
m1_vb3 <- m1_vcov["X1_partial:METRO", "X1_partial:METRO"]
m1_cov <- m1_vcov["METRO", "X1_partial:METRO"]
m1_mg_se <- sqrt(m1_vb1 + m1_vb3 + 2 * m1_cov)

#### complete case
c1_vb1 <- c1_vcov["METRO", "METRO"]
c1_vb3 <- c1_vcov["X1_partial:METRO", "X1_partial:METRO"]
c1_cov <- c1_vcov["METRO", "X1_partial:METRO"]
c1_mg_se <- sqrt(c1_vb1 + c1_vb3 + 2 * c1_cov)

### append extra slope
g1 <- g1 |> rbind("mg" = c(g1_mest, g1_mg_se, NA, NA))
n1 <- n1 |> rbind("mg" = c(n1_mest, n1_mg_se, NA, NA))
m1 <- m1 |> rbind("mg" = c(m1_mest, m1_mg_se, NA, NA))
c1 <- c1 |> rbind("mg" = c(c1_mest, c1_mg_se, NA, NA))

one_mile <- rbind(g1, n1, m1, c1) |>
  mutate(model = c( ## add which model the estimate came from
    rep("g", times = nrow(g1)),
    rep("n", times = nrow(n1)),
    rep("m", times = nrow(m1)),
    rep("c", times = nrow(c1))))
one_mile <- one_mile |>
  mutate(coef = row.names(one_mile)) ## label each row with its parameter
one_mile <- one_mile |>
  ## recode model as factor
  mutate(model = factor(model,
                        levels = c("n", "g", "m", "c"),
                        labels = c("Naive", "Gold Standard",
                                   "MLE", "Complete Case"))) |>
  ## code indicators for interest in specific rows
  mutate(non_metro = coef %in% c("X1", "Xstar1", "X1_partial", "X1_partial1"),
         metro = coef %in% c("mg", "mg1", "mg2", "mg3")) |>
  mutate(img = ifelse(metro, emojis["metro"], emojis["nonmetro"])) |> ## grab images
  mutate(pr = exp(Estimate)) |> ## exponentiate for prevalence ratio (pr)
  mutate(ub = exp(Estimate + 1.96 * `Std. Error`), ## upper CI bound for pr
         lb = exp(Estimate - 1.96 * `Std. Error`)) ## lower CI bound for pr

## five mile data

### save slopes in metro group
g5_mest <- g5["X5", "Estimate"] + g5["X5:METRO", "Estimate"]
n5_mest <- n5["Xstar5", "Estimate"] + n5["Xstar5:METRO", "Estimate"]
m5_mest <- m5["X5_partial", "Estimate"] + m5["X5_partial:METRO", "Estimate"]
c5_mest <- c5["X5_partial", "Estimate"] + c5["X5_partial:METRO", "Estimate"]

### save metro group standard error

#### gold standard
g5_vb1 <- g5_vcov["METRO", "METRO"]
g5_vb3 <- g5_vcov["X5:METRO", "X5:METRO"]
g5_cov <- g5_vcov["METRO", "X5:METRO"]
g5_mg_se <- sqrt(g5_vb1 + g5_vb3 + 2 * g5_cov)

#### naive
n5_vb1 <- n5_vcov["METRO", "METRO"]
n5_vb3 <- n5_vcov["Xstar5:METRO", "Xstar5:METRO"]
n5_cov <- n5_vcov["METRO", "Xstar5:METRO"]
n5_mg_se <- sqrt(n5_vb1 + n5_vb3 + 2 * n5_cov)

#### mle
m5_vb1 <- m5_vcov["METRO", "METRO"]
m5_vb3 <- m5_vcov["X5_partial:METRO", "X5_partial:METRO"]
m5_cov <- m5_vcov["METRO", "X5_partial:METRO"]
m5_mg_se <- sqrt(m5_vb1 + m5_vb3 + 2 * m5_cov)

#### complete case
c5_vb1 <- c5_vcov["METRO", "METRO"]
c5_vb3 <- c5_vcov["X5_partial:METRO", "X5_partial:METRO"]
c5_cov <- c5_vcov["METRO", "X5_partial:METRO"]
c5_mg_se <- sqrt(c5_vb1 + c5_vb3 + 2 * c5_cov)

### append extra slope
g5 <- g5 |> rbind("mg" = c(g5_mest, g5_mg_se, NA, NA))
n5 <- n5 |> rbind("mg" = c(n5_mest, n5_mg_se, NA, NA))
m5 <- m5 |> rbind("mg" = c(m5_mest, m5_mg_se, NA, NA))
c5 <- c5 |> rbind("mg" = c(c5_mest, c5_mg_se, NA, NA))

five_mile <- rbind(g5, n5, m5, c5) |>
  mutate(model = c( ## add which model the estimate came from
    rep("g", times = nrow(g5)),
    rep("n", times = nrow(n5)),
    rep("m", times = nrow(m5)),
    rep("c", times = nrow(c5))))
five_mile <- five_mile |>
  mutate(coef = row.names(five_mile)) ## label each row with its parameter
five_mile <- five_mile |>
  ## recode model as factor
  mutate(model = factor(model,
                        levels = c("n", "g", "m", "c"),
                        labels = c("Naive", "Gold Standard",
                                   "MLE", "Complete Case"))) |>
  ## code indicators for interest in specific rows
  mutate(non_metro = coef %in% c("X5", "Xstar5", "X5_partial", "X5_partial1"),
         metro = coef %in% c("mg", "mg1", "mg2", "mg3")) |>
  mutate(img = ifelse(metro, emojis["metro"], emojis["nonmetro"])) |> ## grab images
  mutate(pr = exp(Estimate)) |> ## exponentiate for prevalence ratio (pr)
  mutate(ub = exp(Estimate + 1.96 * `Std. Error`), ## upper CI bound for pr
         lb = exp(Estimate - 1.96 * `Std. Error`)) ## lower CI bound for pr

## half mile data

### save slopes in metro group
gh_mest <- gh["Xh", "Estimate"] + gh["Xh:METRO", "Estimate"]
nh_mest <- nh["Xstarh", "Estimate"] + nh["Xstarh:METRO", "Estimate"]
mh_mest <- mh["Xh_partial", "Estimate"] + mh["Xh_partial:METRO", "Estimate"]
ch_mest <- ch["Xh_partial", "Estimate"] + ch["Xh_partial:METRO", "Estimate"]

### save metro group standard error

#### gold standard
gh_vb1 <- gh_vcov["METRO", "METRO"]
gh_vb3 <- gh_vcov["Xh:METRO", "Xh:METRO"]
gh_cov <- gh_vcov["METRO", "Xh:METRO"]
gh_mg_se <- sqrt(gh_vb1 + gh_vb3 + 2 * gh_cov)

#### naive
nh_vb1 <- nh_vcov["METRO", "METRO"]
nh_vb3 <- nh_vcov["Xstarh:METRO", "Xstarh:METRO"]
nh_cov <- nh_vcov["METRO", "Xstarh:METRO"]
nh_mg_se <- sqrt(nh_vb1 + nh_vb3 + 2 * nh_cov)

#### mle
mh_vb1 <- mh_vcov["METRO", "METRO"]
mh_vb3 <- mh_vcov["Xh_partial:METRO", "Xh_partial:METRO"]
mh_cov <- mh_vcov["METRO", "Xh_partial:METRO"]
mh_mg_se <- sqrt(mh_vb1 + mh_vb3 + 2 * mh_cov)

#### complete case
ch_vb1 <- ch_vcov["METRO", "METRO"]
ch_vb3 <- ch_vcov["Xh_partial:METRO", "Xh_partial:METRO"]
ch_cov <- ch_vcov["METRO", "Xh_partial:METRO"]
ch_mg_se <- sqrt(ch_vb1 + ch_vb3 + 2 * ch_cov)

### append extra slope
gh <- gh |> rbind("mg" = c(gh_mest, gh_mg_se, NA, NA))
nh <- nh |> rbind("mg" = c(nh_mest, nh_mg_se, NA, NA))
mh <- mh |> rbind("mg" = c(mh_mest, mh_mg_se, NA, NA))
ch <- ch |> rbind("mg" = c(ch_mest, ch_mg_se, NA, NA))

half_mile <- rbind(gh, nh, mh, ch) |>
  mutate(model = c( ## add which model the estimate came from
    rep("g", times = nrow(gh)),
    rep("n", times = nrow(nh)),
    rep("m", times = nrow(mh)),
    rep("c", times = nrow(ch))))
half_mile <- half_mile |>
  mutate(coef = row.names(half_mile)) ## label each row with its parameter
half_mile <- half_mile |>
  ## recode model as factor
  mutate(model = factor(model,
                        levels = c("n", "g", "m", "c"),
                        labels = c("Naive", "Gold Standard",
                                   "MLE", "Complete Case"))) |>
  ## code indicators for interest in specific rows
  mutate(non_metro = coef %in% c("Xh", "Xstarh", "Xh_partial", "Xh_partial1"),
         metro = coef %in% c("mg", "mg1", "mg2", "mg3")) |>
  mutate(img = ifelse(metro, emojis["metro"], emojis["nonmetro"])) |> ## grab images
  mutate(pr = exp(Estimate)) |> ## exponentiate for prevalence ratio (pr)
  mutate(ub = exp(Estimate + 1.96 * `Std. Error`), ## upper CI bound for pr
         lb = exp(Estimate - 1.96 * `Std. Error`)) ## lower CI bound for pr


# make (paper) plots!
facet_labels <- c(t = "Metro Tracts", f = "Non Metro Tracts")

one_mile_plot <- one_mile |>
  ## grab only the correct coefficients and set up to facet
  filter(metro | non_metro) |>
  ## build plot
  mutate(metro = ifelse(metro, "t", "f")) |>
  ggplot(aes(x = model, y = pr)) +
  geom_hline(yintercept = 1, color = 'gray',
             linewidth = 1, linetype = "dashed") +
  geom_point(color = "#56B4E9") +
  geom_errorbar(aes(ymax = ub, ymin = lb),
                width = 0.15,
                linewidth = 1,
                color = "#56B4E9") +
  facet_wrap(vars(metro), labeller = labeller(metro = facet_labels),
             axes = "all") +
  theme_classic() +
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1), "cm"),
        plot.title = element_text(hjust = 0.5),
        axis.line.y.left = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10),
                                    size = 8),
        axis.title.y = element_text(size = 8),
        axis.text = element_text(size = 8),
        strip.background = element_rect(fill = "gray")) +
  labs(title = "",
       x = "Analysis Method",
       y = "Prevalence Ratio")

half_mile_plot <- half_mile |>
  ## grab only the correct coefficients and set up to facet
  filter(metro | non_metro) |>
  ## build plot
  mutate(metro = ifelse(metro, "t", "f")) |>
  ggplot(aes(x = model, y = pr)) +
  geom_hline(yintercept = 1, color = 'gray',
             linewidth = 1, linetype = "dashed") +
  geom_point(color = "#56B4E9") +
  geom_errorbar(aes(ymax = ub, ymin = lb),
                width = 0.15,
                linewidth = 1,
                color = "#56B4E9") +
  facet_wrap(vars(metro), labeller = labeller(metro = facet_labels),
             axes = "all") +
  theme_classic() +
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1), "cm"),
        plot.title = element_text(hjust = 0.5),
        axis.line.y.left = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10),
                                    size = 8),
        axis.title.y = element_text(size = 8),
        axis.text = element_text(size = 8),
        strip.background = element_rect(fill = "gray")) +
  labs(title = "",
       x = "Analysis Method",
       y = "Prevalence Ratio")

# make (cutesy) plots!

one_mile_cutesy <- one_mile |>
  ggplot(aes(x = model, y = pr)) +
  geom_hline(yintercept = 1, color = "gray",
             linewidth = 1, linetype = "dashed") +
  ## display CIs for non-metro pr (on right)
  geom_errorbar(data = one_mile |> dplyr::filter(non_metro),
                aes(ymax = ub, ymin = lb),
                position = position_nudge(x = 0.15),
                width = 0.15,
                linewidth = 1,
                color = "#FFA105") +
  ## display point estimates for non-metro pr (on right)
  geom_image(data = one_mile |> dplyr::filter(non_metro),
             aes(image = img),
             position = position_nudge(x = 0.15),
             size = 0.04) +
  ## display CIs for metro pr (on left)
  geom_errorbar(data = one_mile |> dplyr::filter(metro),
                aes(ymax = ub, ymin = lb),
                position = position_nudge(x = -0.15),
                width = 0.15,
                linewidth = 1,
                color = "#FFA105") +
  ## display point estimates for metro pr (on left)
  geom_image(data = one_mile |> dplyr::filter(metro),
             aes(image = img),
             position = position_nudge(x = -0.15),
             size = 0.04) +
  theme_classic() +
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1), "cm"),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10),
                                    size = 20),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 20)) +
  labs(title = "",
       x = "Analysis Method",
       y = "Prevalence Ratio")


five_mile_cutesy <- five_mile |>
  ggplot(aes(x = model, y = pr)) +
  geom_hline(yintercept = 1, color = "gray",
             linewidth = 1, linetype = "dashed") +
  ## display CIs for non-metro pr (on right)
  geom_errorbar(data = five_mile |> dplyr::filter(non_metro),
                aes(ymax = ub, ymin = lb),
                position = position_nudge(x = 0.15),
                width = 0.15,
                linewidth = 1,
                color = "#FFA105") +
  ## display point estimates for non-metro pr (on right)
  geom_image(data = five_mile |> dplyr::filter(non_metro),
             aes(image = img),
             position = position_nudge(x = 0.15),
             size = 0.04) +
  ## display CIs for metro pr (on left)
  geom_errorbar(data = five_mile |> dplyr::filter(metro),
                aes(ymax = ub, ymin = lb),
                position = position_nudge(x = -0.15),
                width = 0.15,
                linewidth = 1,
                color = "#FFA105") +
  ## display point estimates for metro pr (on left)
  geom_image(data = five_mile |> dplyr::filter(metro),
             aes(image = img),
             position = position_nudge(x = -0.15),
             size = 0.04) +
  theme_classic() +
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1), "cm"),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10),
                                    size = 20),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 20)) +
  labs(title = "",
       x = "Analysis Method",
       y = "Prevalence Ratio")

# function to save plots
save_plot <- function(plot, path = paste0(local_dir, plot_path),
                      w = 10, h = 7){
  fn <- paste0(deparse(substitute(plot)), ".pdf")
  ggsave(filename = fn,
         plot = plot,
         path = path,
         width = w,
         height = h,
         units = "in")
}

#save plots
save_plot(one_mile_plot, w = 6, h = 4)
save_plot(five_mile_plot)
##############################################################################
###         MAPPING CODE                                                  ####
##############################################################################

## Define Piedmont Triad counties
piedmont_triad = c("SURRY", "STOKES", "ROCKINGHAM", "CASWELL",
                   "YADKIN", "FORSYTH", "GUILFORD", "ALAMANCE",
                   "DAVIE", "DAVIDSON", "RANDOLPH", "MONTGOMERY")
## load census tract shape files
tracts = get_acs(state = "NC",
                 geography = "tract",
                 county = piedmont_triad,
                 variables = "B19013_001",
                 geometry = TRUE,
                 year = 2015) |>
  mutate(GEOID = as.double(GEOID)) ### so we can join on the food access data
county = get_acs(state = "NC",
                 geography = "county",
                 county = piedmont_triad,
                 variables = "B19013_001",
                 geometry = TRUE,
                 year = 2015) |>
  mutate(GEOID = as.double(GEOID)) ### so we can join on the food access data

## add access coding
data <- data |>
  mutate(col_half = case_when(
    Xh == 0 & Xstarh == 1 ~ "Error-Prone Only",
    Xh == 1 & Xstarh == 0 ~ "True Only",
    Xh + Xstarh == 2 ~ "Both",
    Xh + Xstarh == 0 ~ "Neither"
  ),
  col_1 = case_when(
    X1 == 0 & Xstar1 == 1 ~ "Error-Prone Only",
    X1 == 1 & Xstar1 == 0 ~ "True Only",
    X1 + Xstar1 == 2 ~ "Both",
    X1 + Xstar1 == 0 ~ "Neither"
  ))


## create maps

#NOTE: I messed with this 3/20/25
half_mile_error_map <- data |> left_join(tracts) |>
  ggplot() +
  geom_sf(aes(geometry = geometry),
          color = "black", fill = NA, size = 25) +
  geom_sf(aes(fill = as.factor(col_half), geometry = geometry)) +
  coord_sf(default_crs = sf::st_crs(4326)) +
  ggspatial::annotation_scale(location = "br",
                              unit_category = "imperial",
                              pad_y = unit(2, "cm"),
                              pad_x = unit(6, "mm"),
                              bar_cols = c("slategray1", "slategray4")) +
  theme_void() +
  theme(legend.position = "none", #blue neither, gold wrong, black both
        plot.title = element_text(hjust = 0.5, size = 10)) +
  labs(fill = "", title = "") +
  theme(plot.margin = margin(2, 2, 2, 2, "mm")) +
  scale_fill_colorblind()



one_mile_error_map <- data |> left_join(tracts) |>
  ggplot() +
  geom_sf(aes(geometry = geometry),
          color = "black", fill = NA, size = 25) +
  geom_sf(aes(fill = as.factor(col_1), geometry = geometry)) +
  theme_void() +
  theme(legend.position = "none", #blue neither, gold wrong, black both
        plot.title = element_text(hjust = 0.5, size = 10)) +
  labs(fill = "", title = "") +
  theme(plot.margin = margin(2, 2, 2, 2, "mm")) +
  scale_fill_colorblind()


## save maps
save_plot(half_mile_error_map,
          path = paste0(local_dir,
                        "food_access_misclassification/Analysis/plots/Maps"))
save_plot(one_mile_error_map,
          path = paste0(local_dir,
                        "food_access_misclassification/Analysis/plots/Maps"))
