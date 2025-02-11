################################################################################
# SETUP ########################################################################
################################################################################
## Load libraries
library(dplyr) ## to wrangle data

# Load data
## Dataset (1) Health outcomes from 2022 PLACES data
health = read.csv(file = "https://raw.githubusercontent.com/sarahlotspeich/food_access_imputation/main/piedmont-triad-data/disease_prevalences_2022.csv")

## Dataset (2) Proximity to health foods based on straight-line and map-based distances (census tracts)
food_access = read.csv(file = "https://raw.githubusercontent.com/sarahlotspeich/food_access_imputation/main/piedmont-triad-data/review_proximity_healthy_foods.csv") |>
  right_join(health,
             by = join_by(LocationID == TractFIPS))


#on Jan 29, use this to generate data (sent to other file)
#food_access <- read.csv(file = "https://raw.githubusercontent.com/sarahlotspeich/food_access_imputation/refs/heads/main/piedmont-triad-data/analysis_data.csv")

## Define access indicators (at 5 mile buffer)
buffer = 5
food_access = food_access |>
  mutate(binX_full = as.numeric(dist_closest_map <= buffer),
         binXstar = as.numeric(dist_closest_straight <= buffer))

## Dataset (3) Rural urban
ruca = read.csv("https://raw.githubusercontent.com/sarahlotspeich/food_access_imputation/main/piedmont-triad-data/ruca2010revised.csv")
ruca_labels = c("Metropolitan", "Metropolitan", "Metropolitan",
                "Micropolitan", "Micropolitan", "Micropolitan",
                "Small town", "Small town", "Small town", "Rural")
food_access = food_access |>
  dplyr::left_join(ruca, by = dplyr::join_by(LocationID == StateCountyTract)) |>
  dplyr::mutate(PrimaryRUCA = factor(PrimaryRUCA,
                                     labels = ruca_labels))

## Remove additional food access columns (not needed for primary analysis)
food_access = food_access |>
  dplyr::select(LocationID, CountyName, ### keep neighborhood identifiers,
                binX_full, binXstar, ### food access indicators,
                PopulationDensity, LandArea, ### additional covariates,
                POP, DIABETES) ### and outcome (with offset).

## Order by Location ID
food_access = food_access |>
  arrange(LocationID)

## Define query indicators (= 1 if map-based measures are available, = 0 otherwise)
set.seed(918) ### make the sampling reproducible
n = 180
queried_subset = food_access |>
  filter(binXstar == 1) |>
  pull(LocationID) |>
  sample(size = n, replace = FALSE)

### Create column for queried X from complete-case
food_access = food_access |>
  mutate(binX_partial = ifelse(test = LocationID %in% queried_subset, #$LocationID,
                               yes = binX_full,
                               no = NA))

# Save data

local_dir <- "/Users/ashleymullan/Documents/Food-Access/"#"local-dir-here/" #replace w/ local location of this repo
setwd(paste0(local_dir, "food_access_misclassification"))

food_access |>
  mutate(CountyName = factor(x = CountyName,
                             levels = c("Alamance", "Caswell", "Davidson", "Davie",
                                        "Forsyth", "Guilford", "Montgomery", "Randolph",
                                        "Rockingham", "Stokes", "Surry", "Yadkin" ),
                             labels = LETTERS[1:12])) |>
  write.csv("Analysis/piedmont_data5.csv",
            row.names = FALSE)

