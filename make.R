
# Create figures and table folders
dir.create("Figures/Main_Figs", recursive = TRUE)
dir.create("Tables")

# Compare regions within in situ dataset
source("Time_Series_Functions.R")
source("Regional_Comparison_InSitu.R", echo = T)

# Download and compare regions within SSST dataset
source("SSST_Data_Download.R", echo = T)
source("Regional_Comparison_SSST.R", echo = T)

# Compare in situ and SSST thermal metrics
source("InSitu_vs_SSST.R", echo = T)
source("MMM_and_DHW_analyses.R", echo = T)
source("Variation_Analysis.R", echo = T)
source("Average_comparisons.R", echo = T)

# Create summary plot of metric comparisons
source("All_Correlations_Plot.R", echo = T)
source("Palau_Map.R", echo = T)

# Run state-space models to partition variance
source("State_Space_Models.R")
