
library(tidyverse)
library(MARSS)

# Read in and reformat data -------------------------------

avg_temps <- read.csv("Data/Merged_Average_Temperatures.csv")

# Convert to wide format time series
temps_wide <- avg_temps |>  
  mutate(Date = as.Date(Date, format = "%Y-%m-%d")) |>
  filter(Region != "All") |> 
  mutate(Region_Origin = paste(Region, Origin, sep = ":")) |> 
  dplyr::select(Date, Region_Origin, Mean_Temp) |> 
  pivot_wider(names_from = "Region_Origin", values_from = "Mean_Temp")

# Add in missing time steps
date_seq <- as.character(seq(min(temps_wide$Date), max(temps_wide$Date), by = "day"))

temps_wide <- temps_wide |> 
  mutate(Date = factor(as.character(Date), levels = date_seq)) |> 
  complete(Date) |>
  mutate(Date = as.Date(as.character(Date), format = "%Y-%m-%d")) |> 
  arrange(Date)

temps_wide <- ts(
  as.matrix(temps_wide[,-1]), 
  start = min(temps_wide$Date), end = max(temps_wide$Date), 
  frequency = 1
)

# Scale to unit variance to simplify model
temps_wide <- scale(temps_wide)

# Dynamic factor model ------------------------------------

## Model 1 ------------------------------------------------

# Here, latent factors are regions
# i.e., the model assumes each region has it's own trend 
# (with shared variability among them) and the in-situ and SSST
# data are measuring those trends imperfectly

# Time series and factor names
ts_names <- colnames(temps_wide)[grepl(":", colnames(temps_wide))]
fct_names <- c("F_Southern", "F_Western", "F_Northern")

# Loadings of regional factors on to time series
Z <- matrix(
  list(
    "z_si", 0, 0,
    0, "z_wi", 0,
    0, 0, "z_ni",
    "z_ss", 0, 0,
    0, "z_ws", 0,
    0, 0, "z_ns"
  ), 
  nrow = length(ts_names), 
  ncol = length(fct_names), 
  dimnames = list(ts_names, fct_names), 
  byrow = TRUE
)

# State equation - each latent factor is a random walk
B <- matrix(
  c(
    1, 0, 0,
    0, 1, 0, 
    0, 0, 1
  ), 
  nrow = length(fct_names), 
  ncol = length(fct_names), 
  dimnames = list(fct_names, fct_names), 
  byrow = TRUE
)

# Bind with other settings
marss_spec <- list(
  Z = Z, # Loadings of factors onto observations
  B = B, # Factor process model
  A = "zero", # offset / scaling 
  R = "diagonal and unequal", # Observation errors
  Q = "equalvarcov", # Factor variances
  U = "zero", # Drift in factor time series
  D = "zero", # Observation covariate effects
  d = "zero", # Ocbservation covariates
  C = "zero", # Factor covariate effects
  c = "zero", # Factor covariates
  x0 = "unequal" # Initial deviation from mean
)

# Run
marss_regional <- MARSS(t(temps_wide), model = marss_spec)

## Model 2 ------------------------------------------------

# Here, latent factors are data source
# i.e., the model assumes in-situ and SSST data have separate
# underlying trends (i.e., they are measuring different water bodies)
# and the observations at the region level are imperfect measurements
# of those trends

# Factor names
fct_names <- c("F_InSitu", "F_SSST")

# Loadings of dataset factors on to time series
Z <- matrix(
  list(
    "z_is", 0,
    "z_iw", 0,
    "z_in", 0,
    0, "z_ss",
    0, "z_sw",
    0, "z_sn"
  ), 
  nrow = length(ts_names), 
  ncol = length(fct_names), 
  dimnames = list(ts_names, fct_names), 
  byrow = TRUE
)

# State equation - each latent factor is a random walk
B <- matrix(
  c(
    1, 0,
    0, 1
  ), 
  nrow = length(fct_names), 
  ncol = length(fct_names), 
  dimnames = list(fct_names, fct_names), 
  byrow = TRUE
)

# Bind with other settings
marss_spec <- list(
  Z = Z, # Loadings of factors onto observations
  B = B, # Factor process model
  A = "zero", # offset / scaling 
  R = "diagonal and unequal", # Observation errors
  Q = "equalvarcov", # Factor variances
  U = "zero", # Drift in factor time series
  D = "zero", # Observation covariate effects
  d = "zero", # Ocbservation covariates
  C = "zero", # Factor covariate effects
  c = "zero", # Factor covariates
  x0 = "unequal" # Initial deviation from mean
)

# Run
marss_datasource <- MARSS(t(temps_wide), model = marss_spec, control = list(maxit = 3000))

# Variance decomposition ----------------------------------

# The model suggests in-situ data are more correlated between regions
# than they are to SSST data within regions - let's just confirm this
# with a quick correlation test

# Correlation between in-situ and SSST data within regions
cor_within <- cor(
  temps_wide[,c("Southern:In_Situ", "Western:In_Situ", "Northern:In_Situ")], 
  temps_wide[,c("Southern:SSST", "Western:SSST", "Northern:SSST")], 
  use = "pairwise.complete.obs"
)
cor_within <- mean(cor_within[lower.tri(cor_within)])

# Correlation between in-situ datasets across regions
cor_across <- cor(
  temps_wide[,c("Southern:In_Situ", "Western:In_Situ", "Northern:In_Situ")], 
  use = "pairwise.complete.obs"
)
cor_across <- mean(cor_across[lower.tri(cor_across)])

# Extract variance terms for data source model
V_global <- marss_datasource$par$Q["offdiag",] # Shared variance of in-situ and SSST trend across Palau
V_source <- marss_datasource$par$Q["diag",] - V_global
V_random <- marss_datasource$par$R[,1] # Observation error
V_total <- V_random + V_source + V_global # Total variance for each region

Vpart <- data.frame(
  Model = "data source (2 state)",
  Region = vapply(ts_names, \(x) strsplit(x, ":")[[1]][1], character(1)), 
  Source = vapply(ts_names, \(x) strsplit(x, ":")[[1]][2], character(1)), 
  Variance_Total = V_total,
  Prop_Shared = V_global / V_total,
  Prop_Source = V_source / V_total,
  Prop_Random = V_random / V_total
)

# Extract variance terms for regional model
V_global <- marss_regional$par$Q["offdiag",] # Shared variance of in-situ and SSST trend across Palau
V_region <- marss_regional$par$Q["diag",] - V_global
V_random <- marss_regional$par$R[,1] # Observation error
V_total <- V_random + V_region + V_global # Total variance for each region

Vpart <- bind_rows(Vpart, data.frame(
  Model = "region (3 state)",
  Region = vapply(ts_names, \(x) strsplit(x, ":")[[1]][1], character(1)), 
  Source = vapply(ts_names, \(x) strsplit(x, ":")[[1]][2], character(1)), 
  Variance_Total = V_total,
  Prop_Shared = V_global / V_total,
  Prop_Region = V_region / V_total,
  Prop_Random = V_random / V_total
))

write.csv(Vpart, "Tables/variance_partitioning.csv", row.names = FALSE)

# Plots ---------------------------------------------------

## Fits to data -------------------------------------------

plotdata <- autoplot(marss_datasource, "fitted.ytT")@data

plot_model_source <- plotdata |> 
  separate(.rownames, into = c("Region", "Source"), sep = ":") |> 
  mutate(Source = gsub("_", " ", Source)) |> 
  mutate(Date = as.Date(date_seq[t], format = "%Y-%m-%d")) |> 
  ggplot(aes(Date, y)) + 
  geom_point(aes(color = "observed"), na.rm = TRUE) + 
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.5) +
  geom_line(aes(y = .fitted, color = "fitted")) + 
  scale_color_manual(values = c(fitted = "black", observed = "dodgerblue3")) +
  facet_grid(Region ~ Source) + 
  ylab("Scaled Temperature") + 
  theme_bw() + 
  theme(strip.background = element_blank(), 
        strip.text = element_text(face = "bold"), 
        text = element_text(size = 14),
        legend.position = "bottom", 
        legend.title = element_blank())

ggsave("Figures/source_factor_model_fits.tif", plot_model_source, height = 5, width = 8, units = "in", dpi = 500)

plotdata <- autoplot(marss_regional, "fitted.ytT")@data

plot_model_region <- plotdata |> 
  separate(.rownames, into = c("Region", "Source"), sep = ":") |> 
  mutate(Source = gsub("_", " ", Source)) |> 
  mutate(Date = as.Date(date_seq[t], format = "%Y-%m-%d")) |> 
  ggplot(aes(Date, y)) + 
  geom_point(aes(color = "observed"), na.rm = TRUE) + 
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.5) +
  geom_line(aes(y = .fitted, color = "fitted")) + 
  scale_color_manual(values = c(fitted = "black", observed = "dodgerblue3")) +
  facet_grid(Region ~ Source) + 
  ylab("Scaled Temperature") + 
  theme_bw() + 
  theme(strip.background = element_blank(), 
        strip.text = element_text(face = "bold"), 
        text = element_text(size = 14),
        legend.position = "bottom", 
        legend.title = element_blank())

ggsave("Figures/region_factor_model_fits.tif", plot_model_region, height = 5, width = 8, units = "in", dpi = 500)

## Latent states ------------------------------------------

plotdata <- autoplot(marss_datasource, "xtT")@data

plot_states <- plotdata |> 
  mutate(
    State = ifelse(.rownames == "X1", "In Situ", "SSST"), 
    Date = as.Date(date_seq[t], format = "%Y-%m-%d")
  ) |> 
  ggplot(aes(Date, .estimate)) + 
  geom_ribbon(aes(ymin = .conf.low, ymax = .conf.up, fill = State), alpha = 0.5) +
  geom_line(aes(color = State)) + 
  scale_color_manual(values = c("In Situ" = "tomato4", "SSST" = "darkslategray4")) + 
  scale_fill_manual(values = c("In Situ" = "tomato4", "SSST" = "darkslategray4")) + 
  labs(y = "Scaled Temperature", color = "Latent State", fill = "Latent State") + 
  theme_bw() + 
  theme(strip.background = element_blank(), 
        strip.text = element_text(face = "bold"), 
        text = element_text(size = 14),
        legend.position = "bottom")

ggsave("Figures/source_factor_model_states.tif", plot_states, height = 4, width = 8, units = "in", dpi = 500)

## Variance partitioning by each dataset ------------------

plot_Vpart <- Vpart |> 
  filter(Model == "data source (2 state)") |> 
  dplyr::select(-Prop_Region) |> 
  pivot_longer(starts_with("Prop"), names_prefix = "Prop_", names_to = "type", values_to = "var") |>
  mutate(
    Source = gsub("_", " ", Source),
    type = factor(case_when(
      type == "Random" ~ "Obs. error", 
      type == "Source" ~ "Origin", 
      type == "Shared" ~ "Inter-origin"
    ), levels = c("Obs. error", "Inter-origin", "Origin"))
  ) |> 
  ggplot(aes(Region, var, fill = type)) + 
  geom_bar(stat = "identity") + 
  facet_wrap(~Source, ncol = 1) +
  coord_flip(expand = FALSE) + 
  labs(y = "Proportion Variance Explained", fill = "Source") + 
  scale_fill_manual(values = c("orangered3", "goldenrod3", "darkgreen")) + 
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_minimal() + 
  theme(strip.text = element_text(face = "bold", hjust = 0), 
        text = element_text(size = 14),
        legend.position = "bottom")

ggsave("Figures/variance_partitioning.tif", plot_Vpart, height = 4, width = 8, units = "in", dpi = 500)
