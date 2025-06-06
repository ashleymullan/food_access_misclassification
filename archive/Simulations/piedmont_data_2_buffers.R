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
## Dataset (3) Rural urban
ruca = read.csv("https://raw.githubusercontent.com/sarahlotspeich/food_access_imputation/main/piedmont-triad-data/ruca2010revised.csv")
ruca_labels = c("Metropolitan", "Metropolitan", "Metropolitan",
                "Micropolitan", "Micropolitan", "Micropolitan",
                "Small town", "Small town", "Small town", "Rural")

for(buffer in c(1,5)){

## Define access indicators (at 1 or 5 mile buffer)
  food_access_bin = food_access |>
    mutate(binX_full = as.numeric(dist_closest_map <= buffer),
         binXstar = as.numeric(dist_closest_straight <= buffer))


  food_access_bin = food_access_bin |>
  dplyr::left_join(ruca, by = dplyr::join_by(LocationID == StateCountyTract)) |>
  dplyr::mutate(PrimaryRUCA = factor(PrimaryRUCA,
                                     labels = ruca_labels))

  ## Remove additional food access columns (not needed for primary analysis)
  food_access_bin = food_access_bin |>
  dplyr::select(LocationID, CountyName, ### keep neighborhood identifiers,
                binX_full, binXstar, ### food access indicators,
                LandArea, PopulationDensity, ### additional covariates,
                POP, DIABETES) ### and outcome (with offset).

  ## Order by Location ID
  food_access_bin = food_access_bin |>
    arrange(LocationID)

  ## Define query indicators (= 1 if map-based measures are available, = 0 otherwise)
  set.seed(918) ### make the sampling reproducible
  n = 4 ### set number of census tracts sampled from each county to be queried
  queried_subset = food_access_bin |>
    group_by(CountyName) |>
    sample_n(size = n, replace = FALSE)

  ### Create column for queried X from complete-case
  food_access_bin = food_access_bin |>
    mutate(binX_partial = ifelse(test = LocationID %in% queried_subset$LocationID,
                               yes = binX_full,
                               no = NA))

  # Save data
  food_access_bin |>
  mutate(CountyName = factor(x = CountyName,
                             levels = c("Alamance", "Caswell", "Davidson", "Davie",
                                        "Forsyth", "Guilford", "Montgomery", "Randolph",
                                        "Rockingham", "Stokes", "Surry", "Yadkin" ),
                             labels = LETTERS[1:12])) |>
  dplyr::rename(O_POP = POP,
                Y_DIABETES = DIABETES) |>
  write.csv(paste0("piedmont_data_", buffer, ".csv"),
            row.names = FALSE)
}
