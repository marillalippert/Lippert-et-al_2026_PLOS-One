# New analysis to see if any SSST metric is highly correlated with the daily (NOT nightly) variation the hobo's detected. 

library(tidyverse)
library(ggpubr)
library(GGally)

source("Table_Outputs_Function.R")
source("Time_Series_Functions.R")

full_insitu_data <- read.csv("Data/HOBO_Regional_Master.csv", stringsAsFactors = F)

full_insitu_data$Date <- as.Date(full_insitu_data$Date)
full_insitu_data$Region <- as.factor(full_insitu_data$Region)

diurnal_variation_data <- full_insitu_data |> group_by(Region, Date) |> 
  summarise(Day_Mean = mean(Temperature),
            Day_Min = min(Temperature),
            Day_Max = max(Temperature),
            Variation = max(Temperature) - min(Temperature),
            .groups = "drop")

diurnal_variation_data$date_match <- paste(diurnal_variation_data$Date, diurnal_variation_data$Region)
colnames(diurnal_variation_data) <- c("Region", "Date", "Mean_InSitu", "Min_InSitu", "Max_InSitu", "Day_Variation", "date_match")

# Need to re-match ssst dataset since we aren't working with night date anymore
ssst_data_full <- read.csv("Data/Average_SSST_Data.csv", stringsAsFactors = F)
ssst_data <- filter(ssst_data_full, Region != "All")
ssst_data$date <- as.Date(ssst_data$date)
ssst_data$date_match <- paste(ssst_data$date, ssst_data$Region)

ssst_data <- filter(ssst_data, date_match %in% diurnal_variation_data$date_match)
diurnal_variation_data <- filter(diurnal_variation_data, date_match %in% ssst_data$date_match)

ssst_region_avg <- ssst_data[,c(1:6)]
colnames(ssst_region_avg) <- c("Date", "Region", "Mean_SSST", "SD", "Min_SSST", "Max_SSST")

# Run Analysis
## Make long version of data with origin
ssst_region_avg$Origin <- rep("SSST", nrow(ssst_region_avg))
ssst_region_avg$Temp <- ssst_region_avg$Mean_SSST
diurnal_variation_data$Origin <- rep("In_Situ", nrow(diurnal_variation_data))
diurnal_variation_data$Temp <- diurnal_variation_data$Day_Variation

long_data <- rbind(ssst_region_avg[,c(7,2,1,8)], diurnal_variation_data[,c(8,1,2,9)])
long_data <- long_data |> mutate(time_step = 1 + as.integer(Date) - min(as.integer(Date)))
long_data$Region <- as.factor(long_data$Region)
var_mean_corr <- table_outputs(long_data, Temp)
write.csv(var_mean_corr, "Tables/Diurnal_Variation_vs_SSST_Mean_Outputs.csv", row.names = F)


#Combine datasets to plot
merged_data <- merge(diurnal_variation_data[,c(1:3,6)], ssst_region_avg[,c(1:3,5,6)], by = c("Region", "Date"))

#add coordinates to plot r values on plot into var_mean_corr dataframe
var_mean_corr$x <- 28
var_mean_corr$y <- 4.75

regions <- merged_data |> ggplot(aes(Mean_SSST, Day_Variation)) +
  geom_point() +
  stat_smooth( method = "lm", color = "magenta4") +  
  geom_text(data = filter(var_mean_corr, Region != "all"), aes(x = x, y = y, label = paste0("r = ", cor)), color = "magenta4") +
  stat_regline_equation(label.y = 4, color = "magenta4") +
  xlab("Mean SSST") + ylab("Diurnal Variation") + 
  facet_wrap(~Region, ncol = 1) +
  theme_bw() + theme(strip.background = element_blank())
#ggsave("Figures/Diurnal_Variation_vs_Mean_SSST.png", height = 10, width = 14, scale = 0.7)

# now create datasets and combine for all region to plot all panel
all_region_ssst <- ssst_data_full |>filter(Region == "All")
colnames(all_region_ssst) <- c("Date", "Region", "Mean_SSST", "SD", "Min_SSST", "Max_SSST")

all_region_insitu <- full_insitu_data |> group_by(Date) |> 
  summarise(Day_Mean = mean(Temperature),
            Day_Min = min(Temperature),
            Day_Max = max(Temperature),
            Variation = max(Temperature) - min(Temperature))
colnames(all_region_insitu) <- c("Date", "Mean_InSitu", "Min_InSitu", "Max_InSitu", "Day_Variation")

all_merged_data <- merge(all_region_insitu, all_region_ssst, by = "Date")

all_region <- all_merged_data |> ggplot(aes(Mean_SSST, Day_Variation)) +
  geom_point() +
  geom_abline(data = filter(var_mean_corr, Region == "all"), 
              aes(intercept = Estimate[1], slope = Estimate[2]), color = "magenta4", linewidth = 1.5) +
  geom_text(data = filter(var_mean_corr, Region == "all"), aes(x = x, y = y, label = paste0("r = ", cor)), color = "magenta4") +
  geom_text(data = filter(var_mean_corr, Region == "all"), aes(x = x + 0.2, y = y-0.25, label = paste0("y = ", Estimate[1], Estimate[2], "x")), color = "magenta4") +
  xlab("Mean SSST") + ylab("Diurnal Variation") + 
  labs(subtitle = "All") +
  theme_bw() + theme(strip.background = element_blank(), plot.subtitle = element_text(hjust = 0.5))

ggarrange(all_region, regions, nrow = 1)
ggsave("Figures/Main_Figs/Diurnal_Variation_Comp.tif", height = 6, width = 10, scale = 0.7, dpi = 300)

