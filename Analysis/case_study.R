##############################################################################
###         LIBRARIES AND FILE PATHS                                      ####
##############################################################################

# load libraries
library(dplyr) ## for data wrangling
library(ggplot2) ## for graph
library(possum) ## for mle
library(ggthemes) ## for colorblind palette
library(tidycensus) ## for shapefiles
library(auditDesignR) ## for query study design

local_dir <- "/Users/ashleymullan/Documents/Food-Access/" ## change me as needed
plot_path <- "food_access_misclassification/Analysis/plots/"

##############################################################################
###         DATA SETUP CODE                                               ####
##############################################################################

# read in (and save only the necessary) data
data_file <- paste0("https://raw.githubusercontent.com/sarahlotspeich/food_access_imputation/",
                      "refs/heads/main/piedmont-triad-data/analysis_data.csv")
data <- read.csv(data_file) |>
  dplyr::select("GEOID", "CountyName", "Xstar", "X_full",
                "O_POP", "Y_DIABETES", "METRO")

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
data$query_indh <- 0 ### initialize all unqueried
data$query_indh[data$Xstarh == 1] <- sample_cc(dat = data[data$Xstarh == 1, ],
                                               phI = nrow(data[data$Xstarh == 1, ]),
                                               phII = n,
                                               sample_on = "METRO",
                                               prop_cases = 0.5)
data$query_ind1 <- 0 ### initialize all unqueried
data$query_ind1[data$Xstar1 == 1] <- sample_cc(dat = data[data$Xstar1 == 1, ],
                                               phI = nrow(data[data$Xstar1 == 1, ]),
                                               phII = n,
                                               sample_on = "METRO",
                                               prop_cases = 0.5)
data$query_ind5 <- 0 ### initialize all unqueried
data$query_ind5[data$Xstar5 == 1] <- sample_cc(dat = data[data$Xstar5 == 1, ],
                                               phI = nrow(data[data$Xstar5 == 1, ]),
                                               phII = n,
                                               sample_on = "METRO",
                                               prop_cases = 0.5)

data <- data |>
  mutate(
    Xh_partial = case_when(
      Xstarh == 1 & query_indh == 1 ~ Xh,
      Xstarh == 0 ~ 0,
      .default = NA),
    X1_partial = case_when(
      Xstar1 == 1 & query_ind1 == 1 ~ X1,
      Xstar1 == 0 ~ 0,
      .default = NA),
    X5_partial = case_when(
      Xstar5 == 1 & query_ind5 == 1 ~ X5,
      Xstar5 == 0 ~ 0,
      .default = NA)
    )

##############################################################################
###         MODELING CODE                                                 ####
##############################################################################

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

##############################################################################
###         PLOTTING CODE                                                 ####
##############################################################################

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
  mutate(pr = exp(Estimate)) |> ## exponentiate for prevalence ratio (pr)
  mutate(ub = exp(Estimate + 1.96 * `Std. Error`), ## upper CI bound for pr
         lb = exp(Estimate - 1.96 * `Std. Error`)) ## lower CI bound for pr


# make (paper) plots!
facet_labels <- c(t = "Metro Tracts", f = "Non Metro Tracts")

one_and_half_mile_plot <- half_mile |>
  ## grab only the correct coefficients and set up to facet
  filter(metro | non_metro) |>
  mutate(threshold = "Threshold: 1/2 Mile") |>
  bind_rows(
    one_mile |>
      ## grab only the correct coefficients and set up to facet
      filter(metro | non_metro) |>
      mutate(threshold = "Threshold: 1 Mile")
  ) |>
  ## build plot
  mutate(metro = ifelse(metro, "Metropolitan Tracts", "Non-Metropolitan Tracts"),
         threshold = factor(x = threshold,
                            levels = c("Threshold: 1/2 Mile", "Threshold: 1 Mile"))) |>
  ggplot(aes(x = model, y = pr, color = metro)) +
  geom_hline(yintercept = 1, color = 'gray',
             linewidth = 1, linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymax = ub, ymin = lb),
                position = position_dodge(width = 0.5),
                width = 0.15,
                linewidth = 1) +
  scale_color_manual(values = c("#56B4E9", "maroon"), name = "") +
  theme_classic() +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5),
        axis.line.y.left = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10),
                                    size = 8),
        axis.title.y = element_text(size = 8),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(color = "white"),
        legend.box.margin = unit(c(1,1,-0.5,1), "cm"),
        legend.background = element_blank()) +
  labs(title = "",
       x = "Analysis Method",
       y = "Prevalence Ratio (Access vs. No Access)") +
  facet_wrap(~threshold)

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

#save plot
save_plot(one_and_half_mile_plot, w = 6, h = 4)

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

half_mile_error_map <- data |> left_join(tracts) |>
  ggplot() +
  geom_sf(aes(geometry = geometry),
          color = "black", fill = NA, size = 25) +
  geom_sf(aes(fill = as.factor(col_half), geometry = geometry)) +
  coord_sf(default_crs = sf::st_crs(4326)) +
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
