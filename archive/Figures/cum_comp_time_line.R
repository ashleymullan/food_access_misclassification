############################################################################################
## SETUP ///////////////////////////////////////////////////////////////////////////////////
############################################################################################
## Load libraries
library(ggplot2) ## to make plots
library(dplyr) ## for data wrangling

############################################################################################
## LOAD FOOD ACCESS DATA FOR FORSYTH AND BORDERING COUNTIES' CENSUS TRACTS /////////////////
############################################################################################
## Proximity to health foods based on straight-line and map-based distances (census tracts)
food_access = read.csv("https://raw.githubusercontent.com/sarahlotspeich/food_access_imputation/main/piedmont-triad-data/analysis_data.csv")

# Define color palette 
slide_colors = c("#F2958D", "#ADA264", "#EB5F3F", "#3A5724", "#FAB825", "#E67B33") 

############################################################################################
## MAKE LINE GRAPH OF CUMULATIVE COMPUTING TIME ////////////////////////////////////////////
############################################################################################
food_access |> 
  select(LocationID, dplyr::ends_with("_time")) |> 
  mutate(row_id = 1:dplyr::n()) |> 
  tidyr::gather(key = "dist_calc", 
                value = "comp_time", 
                -c(1, 4)) |> 
  group_by(dist_calc) |> 
  mutate(cum_comp_time = cumsum(comp_time), 
         dist_calc = factor(x = dist_calc, 
                            levels = c("X_time", "Xstar_time"), 
                            labels = c("Map-Based", "Straight-Line"))) |> 
  ggplot(aes(x = row_id, 
             y = cum_comp_time, 
             color = dist_calc, 
             alpha = dist_calc)) + 
  geom_line(linewidth = 1) + 
  scale_color_manual(values = slide_colors[c(2:3)], 
                     name = "Distance Calculation:") + 
  scale_alpha_manual(values = c(0, 1), 
                     guide = "none") + 
  theme_minimal(base_size = 10) + 
  theme(plot.margin = margin(l=25, r=20, t=20, b=25), 
        legend.position = "top") + 
  labs(x = "Number of Neighborhoods' Distances Calculated",
       y = "Cumulative Computing Time (in Seconds)") 
ggsave(filename = "~/Documents/food_access_misclassification/Figures/cum_comp_time_straight_line.png", 
       device = "png", 
       width = 7, 
       height = 5, 
       units = "in")

food_access |> 
  select(LocationID, dplyr::ends_with("_time")) |> 
  mutate(row_id = 1:dplyr::n()) |> 
  tidyr::gather(key = "dist_calc", 
                value = "comp_time", 
                -c(1, 4)) |> 
  group_by(dist_calc) |> 
  mutate(cum_comp_time = cumsum(comp_time), 
         dist_calc = factor(x = dist_calc, 
                            levels = c("X_time", "Xstar_time"), 
                            labels = c("Map-Based", "Straight-Line"))) |> 
  ggplot(aes(x = row_id, 
             y = cum_comp_time, 
             color = dist_calc, 
             alpha = dist_calc)) + 
  geom_line(linewidth = 1) + 
  scale_color_manual(values = slide_colors[c(2:3)], 
                     name = "Distance Calculation:") + 
  scale_alpha_manual(values = c(1, 1), 
                     guide = "none") + 
  theme_minimal(base_size = 10) + 
  theme(plot.margin = margin(l=25, r=20, t=20, b=25), 
        legend.position = "top") + 
  labs(x = "Number of Neighborhoods' Distances Calculated",
       y = "Cumulative Computing Time (in Seconds)") 
ggsave(filename = "~/Documents/food_access_misclassification/Figures/cum_comp_time_both_line.png", 
       device = "png", 
       width = 7, 
       height = 5, 
       units = "in")
