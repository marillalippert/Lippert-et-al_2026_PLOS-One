# per request of Steve, we compared daily average in situ temperature (not nightly average)
# we will also create a supplemental figure showing how in situ daily average, in situ nightly average, and ssst nightly average compare

library(tidyverse)
library(ggpubr)
library(GGally)

source("Table_Outputs_Function.R")
source("Time_Series_Functions.R")

full_insitu_data <- read.csv("Data/HOBO_Regional_Master.csv", stringsAsFactors = F)

full_insitu_data$Date <- as.Date(full_insitu_data$Date)
full_insitu_data$Region <- as.factor(full_insitu_data$Region)

all <- full_insitu_data
all$Region <- "All"

full_insitu_data <- rbind(full_insitu_data, all)

diurnal_variation_data <- full_insitu_data |> group_by(Region, Date) |> 
  summarise(Day_Mean = mean(Temperature),
            Day_Min = min(Temperature),
            Day_Max = max(Temperature),
            Variation = max(Temperature) - min(Temperature),
            .groups = "drop")

diurnal_variation_data$date_match <- paste(diurnal_variation_data$Date, diurnal_variation_data$Region)
colnames(diurnal_variation_data) <- c("Region", "Date", "Mean_InSitu", "Min_InSitu", "Max_InSitu", "Day_Variation", "date_match")

# Need to re-match ssst dataset since we aren't working with night date anymore
ssst_data <- read.csv("Data/Average_SSST_Data.csv", stringsAsFactors = F)
#ssst_data <- filter(ssst_data, Region != "All")
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



colnames(diurnal_variation_data)[3] <- "Mean_Temp"
colnames(ssst_region_avg)[3] <- "Mean_Temp"
daily_means <- rbind(ssst_region_avg[,c(7,2,1,3)], diurnal_variation_data[c(8,1,2,3)])
daily_means <- daily_means |> filter(Region != "All")
daily_means <- daily_means |> mutate(time_step = 1 + as.integer(Date) - min(as.integer(Date)))
daily_mean_corr <- table_outputs(daily_means, Mean_Temp)
write.csv(daily_mean_corr, "Tables/Daily_Mean_Outputs.csv", row.names = F)
write.csv(daily_means, "Data/Daily_Means.csv", row.names = F)



## now to make the plot with all 3 comparisons in it
nightly <- read.csv("Data/Merged_average_Temperatures.csv", stringsAsFactors = F)

insitu_all_day_avg <- full_insitu_data |> group_by(Date) |> 
  summarise(Day_Mean = mean(Temperature),
            Day_Min = min(Temperature),
            Day_Max = max(Temperature),
            Variation = max(Temperature) - min(Temperature),
            .groups = "drop")


daily_means$Region <- as.factor(daily_means$Region)
daily_means$Date <- as.Date(daily_means$Date, format = "%Y-%m-%d")
nightly$Region <- as.factor(nightly$Region)
nightly$Date <- as.Date(nightly$Date, format = "%Y-%m-%d")

daily_insitu <- filter(diurnal_variation_data, Origin == "In_Situ" & Region == "All")
nightly_insitu <- filter(nightly, Origin == "In_Situ" & Region == "All")
nightly_ssst <- filter(nightly, Origin == "SSST" & Region == "All")

ggplot() +
  geom_point(data = daily_insitu, aes(x = Date, y = Mean_Temp), color = "orange2") +
  geom_point(data = nightly_insitu, aes(x = Date, y = Mean_Temp), color = "orange4") +
  geom_point(data = nightly_ssst, aes(x = Date, y = Mean_Temp), color = "turquoise4") +
  theme_bw()
  
dfs <- c(daily_insitu, nightly_insitu, nightly_ssst)
wide_data <- merge(daily_insitu, nightly_insitu, by = c("Date", "Region")) |> merge(nightly_ssst, by = c("Date", "Region"))

day_label <- paste0("r = ", round(daily_mean_corr[1,8], digits = 2))

nightly_mean <- read.csv("Tables/Mean_Temp_Outputs.csv", stringsAsFactors = F)
night_label <- paste0("r = ", round(nightly_mean[1,8], digits = 2))

wide_data  |> ggplot() +
  geom_point(aes(x = Mean_Temp, y = Mean_Temp.x, color = "Daily In Situ"), alpha = 0.6, size = 2) +
  geom_point(aes(x = Mean_Temp, y = Mean_Temp.y, color = "Nightly In Situ"), alpha = 0.6, size = 2) +
  scale_color_manual(values = c("orange2", "skyblue4")) +
  stat_smooth(method = "lm", aes(x = Mean_Temp, y = Mean_Temp.y), color = "skyblue4", linewidth = 1.5) +
  stat_smooth(method = "lm", aes(x = Mean_Temp, y = Mean_Temp.x), color = "orange2", linewidth = 1.5) +
  annotate("text", label = day_label, x = 30, y = 29, color = "orange2") +
  annotate("text", label = night_label, x = 28.5, y = 30.5, color = "skyblue4") +
  xlab("Nightly Mean SSST Temperature") + ylab("Mean In Situ Temperature") +
  theme_bw() + theme(legend.title = element_blank(), legend.position = "top", legend.text = element_text(size = 10))
ggsave("Figures/Daily_Nightly_Average_Comparison.tif", height = 8, width = 12, scale = 0.7, dpi = 300)
