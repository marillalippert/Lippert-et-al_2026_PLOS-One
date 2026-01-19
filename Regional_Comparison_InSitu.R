## Refine and select specific sets of HOBO data

library(tidyverse)
library(lubridate)
library(GGally)
library(ggpubr)
library(nlme)

source("Time_Series_Functions.R")

# Filter master HOBO  -------------
full_hobo_data <- read.csv("Data/HOBO_Master(Feb2023).csv", stringsAsFactors = F)
hobo_region <- filter(full_hobo_data, Patch_Reef %in% c(13:39))
hobo_region$Region <- ifelse(hobo_region$Patch_Reef %in% c(13:19), "Southern", 
                             ifelse(hobo_region$Patch_Reef %in% c(20:29), "Western", 
                                    ifelse(hobo_region$Patch_Reef %in% c(30:39), "Northern", NA)))
hobo_region$Region <- factor(hobo_region$Region, ordered = T, levels = c("Southern", "Western", "Northern"))
hobo_region <- filter(hobo_region, hobo_region$Patch_Reef != 21)
hobo_region$Date_Time_PWT <- as.POSIXct(hobo_region$Date_Time_PWT, format = "%Y-%m-%d %H:%M:%OS")

## To make sure we don't have duplicate reads for the same day/time/colony, group aqnd summarise
hobo_region <- hobo_region |> group_by(Date_Time_PWT, Colony, Region, Latitude, Longitude) |> 
  summarise(Temperature = mean(Temperature))
hobo_region <- as.data.frame(hobo_region)
## Separate out date/time variables
hobo_region$Date <- substr(hobo_region$Date_Time_PWT, 1, 10)
hobo_region$Hour <- substr(hobo_region$Date_Time_PWT, 12, 13)
hobo_region$Time <- substr(hobo_region$Date_Time_PWT, 12, 19)
## Colony 119 is troublesome so we filter it's dates manually
hobo_region <- hobo_region |> filter(!(Colony == 119 & Date_Time_PWT >= "2018-10-01 00:00:00 EST" & Date_Time_PWT < "2019-03-01 00:00:00 EST"))

## Make enat cutofff so at start of 2021
hobo_region <- filter(hobo_region, Date_Time_PWT < "2021-01-01 00:00:00 EST")
hobo_region$data_type <- rep("In-Situ", nrow(hobo_region))
hobo_region <- na.omit(hobo_region)
write.csv(hobo_region, "Data/HOBO_Regional_Master.csv", row.names = F)

# get average deployment time for hobo loggers
avg_deployment <- hobo_region |> group_by(Date, Region, Colony) |> 
  summarise(mean = mean(Temperature), stdev = sd(Temperature))

avg_deployment$date_number <- as.numeric(as.Date(avg_deployment$Date))
avg_deployment$date_diff <- rep(0, nrow(avg_deployment))

for(i in 2:nrow(avg_deployment)) {
    avg_deployment$date_diff[i-1] <- 
      avg_deployment$date_number[i-1] - avg_deployment$date_number[i]
}

avg_deployment$Colony <- as.factor(avg_deployment$Colony)
avg_deployment$Date <- as.Date(avg_deployment$Date, format = "%Y-%m-%d")
avg_deployment_ts <- avg_deployment |> 
  group_by(Colony) |> 
  mutate(ts_group = group_ts(Date))

all_avg_deployment <- avg_deployment_ts |> group_by(Colony, ts_group) |> 
  summarise(ts_length = length(ts_group))

average_deployment_time <- mean(all_avg_deployment$ts_length)

# 

hobo_region %>% group_by(Date, Region) %>% 
  summarise(mean = mean(Temperature), stdev = sd(Temperature)) %>% 
  ggplot(aes(Date, mean)) +
  geom_errorbar(aes(ymin = mean - stdev, ymax = mean + stdev)) +
  geom_point(aes(color = Region)) +
  theme_bw()

## Subset data to only encompass temperature reads between 22:00:00 and 02:00:00 
## to get average nighttime data for comparison
night <- c(22, 23, 00, 01, 02)
nighttime_insitu <- filter(hobo_region, hobo_region$Hour %in% night)
write.csv(nighttime_insitu, "Data/Nighttime_Full_InSitu.csv", row.names = F)

nighttime_insitu %>% group_by(Date, Region) %>%  
  summarise(mean = mean(Temperature), stdev = sd(Temperature)) %>% 
  ggplot(aes(Date, mean)) +
  geom_errorbar(aes(ymin = mean - stdev, ymax = mean + stdev)) +
  geom_point(aes(color = Region)) +
  theme_bw()

## get night averages for each region
regional_avg_night_insitu <- nighttime_insitu |> 
  group_by(Region, Colony) |>
  arrange(Date_Time_PWT, .by_group = TRUE) |> 
  mutate(time_block = cumsum(Time == "22:00:00")) |> 
  group_by(Region, Colony, time_block) |> 
  summarize(
    date = min(as.Date(Date_Time_PWT)),
    temp_mean = mean(Temperature),
    temp_min = min(Temperature),
    temp_max = max(Temperature),
    .groups = "drop"
  ) |> 
  group_by(Region, date) |> 
  summarize(mean_night_temp = mean(temp_mean),
            min_night_temp = min(temp_min),
            max_night_temp = max(temp_max),
            night_sd = sd(temp_mean),
            .groups = "drop")

## get night average for whole of Palau
total_avg_night_insitu <- nighttime_insitu |> 
  group_by(Region, Colony) |>
  arrange(Date_Time_PWT, .by_group = TRUE) |> 
  mutate(time_block = cumsum(Time == "22:00:00")) |> 
  group_by(Region, Colony, time_block) |> 
  summarize(
    date = min(as.Date(Date_Time_PWT)),
    temp_mean = mean(Temperature),
    temp_min = min(Temperature),
    temp_max = max(Temperature),
    .groups = "drop"
  ) |> 
  group_by(date) |> 
  summarize(mean_night_temp = mean(temp_mean),
            min_night_temp = min(temp_min),
            max_night_temp = max(temp_max),
            night_sd = sd(temp_mean),
            .groups = "drop")

total_avg_night_insitu$Region <- rep("All", nrow(total_avg_night_insitu))

## merge dataframes together so we have night averages of each region and total together
avg_night_temp_insitu <- rbind(total_avg_night_insitu, regional_avg_night_insitu)
write.csv(avg_night_temp_insitu, "Data/Average_Nighttime_Insitu_Data.csv", row.names = F)

# Regional Comparison ----------
## Average night time -----------
## Because our regions span different dates and have a different amount of days with data, 
## we need to filter the data so all regions cover the same dates
North_dates <- avg_night_temp_insitu$date[avg_night_temp_insitu$Region == "Northern"]
West_dates <- avg_night_temp_insitu$date[avg_night_temp_insitu$Region == "Western"]
South_dates <- avg_night_temp_insitu$date[avg_night_temp_insitu$Region == "Southern"]
regional_comparison_data <- avg_night_temp_insitu[avg_night_temp_insitu$date %in% North_dates & 
                                                    avg_night_temp_insitu$date %in% West_dates &
                                                    avg_night_temp_insitu$date %in% South_dates,]

regional_correlation_data <- regional_comparison_data[,c(1,2,6)] |> pivot_wider(names_from = Region, values_from = mean_night_temp)
insitu_regional_correlation <- cor(regional_correlation_data[,-1])
write.table(insitu_regional_correlation, "Tables/InSitu_Regional_Correlations_Outputs.csv", sep = ",", row.names = T, col.names = T)

ggpairs(regional_correlation_data[,-1], title = "InSitu Regional Correlations")
ggsave("Figures/InSitu_Regional_Correlations.tif", width = 10, height = 8, scale = 0.7, dpi = 300)

# Plot regional comparisons --------
## Comparing nighttime average temps - correlation plots

nw <- regional_correlation_data |> ggplot(aes(Northern, Western)) +
  geom_point() +
  geom_smooth(method = "lm", color = "palegreen4") +
  stat_cor(method = "pearson", aes(label = after_stat(r.label)), color = "palegreen4") +
  stat_regline_equation(label.y = 31, color = "palegreen4") +
  geom_abline(intercept = 0, slope = 1, color = "gray57") +
  labs(subtitle = "North-West Correlation") + 
  theme_bw()

ns <- regional_correlation_data |> ggplot(aes(Northern, Southern)) +
  geom_point() +
  geom_smooth(method = "lm", color = "palegreen4") +
  stat_cor(method = "pearson", aes(label = after_stat(r.label)), color = "palegreen4") +
  stat_regline_equation(label.y = 31, color = "palegreen4") +
  geom_abline(intercept = 0, slope = 1, color = "gray57") +
  labs(subtitle = "North-South Correlation") +
  theme_bw()

sw <- regional_correlation_data |> ggplot(aes(Southern, Western)) +
  geom_point() +
  geom_smooth(method = "lm", color = "palegreen4") +
  stat_cor(method = "pearson", aes(label = after_stat(r.label)), color = "palegreen4") +
  stat_regline_equation(label.y = 31, color = "palegreen4") +
  geom_abline(intercept = 0, slope = 1, color = "gray57") +
  labs(subtitle = "South-West Correlation") +
  theme_bw()

plotlist <- list(nw, ns, sw) 
ggarrange(plotlist = plotlist, common.legend = T, nrow = 1, legend = "bottom")
ggsave("Figures/Regional_Insitu_Correlation.tif", height = 6, width = 12, scale = 0.7, dpi = 300)

## Comparing offsets/differences between regional nighttime averages
regional_correlation_data$NW_diff <- regional_correlation_data$Northern - regional_correlation_data$Western
regional_correlation_data$NS_diff <- regional_correlation_data$Northern - regional_correlation_data$Southern
regional_correlation_data$WS_diff <- regional_correlation_data$Western - regional_correlation_data$Southern


## Get proportion of days where the difference is greater than zero (meaning region i is warmer than region j)
prop_W_warmer_than_N <- 100*round((nrow(filter(regional_correlation_data, NW_diff < 0)))/nrow(regional_correlation_data), 3)
prop_S_warmer_than_N <- 100*round((nrow(filter(regional_correlation_data, NS_diff < 0)))/nrow(regional_correlation_data), 3)
prop_S_warmer_than_W <- 100*round((nrow(filter(regional_correlation_data, WS_diff < 0)))/nrow(regional_correlation_data), 3)

prop_N_warmer_than_W <- 100*round((nrow(filter(regional_correlation_data, NW_diff > 0)))/nrow(regional_correlation_data), 3)
prop_N_warmer_than_S <- 100*round((nrow(filter(regional_correlation_data, NS_diff > 0)))/nrow(regional_correlation_data), 3)
prop_W_warmer_than_S <- 100*round((nrow(filter(regional_correlation_data, WS_diff > 0)))/nrow(regional_correlation_data), 3)


nw_diff <- regional_correlation_data |> ggplot(aes(NW_diff, Northern)) +
  geom_vline(xintercept = 0, color = "darkslategray4", linetype = "dashed", linewidth = 1.2) +
  geom_vline(aes(xintercept = mean(NW_diff)), color = "darkslategray4", linewidth = 1.2) +
  geom_point() +
  annotate("text", x = -1, y = 32.5, label = paste0("West Warmer (", prop_W_warmer_than_N, "% of points)")) +
  annotate("text", x = 0.6, y = 32.5, label = paste0("North Warmer (", prop_N_warmer_than_W, "% of points)")) +
  labs(title = "(a)", subtitle = "North - West Comparison") + 
  xlab (expression(North~"-"~West~Temperature~Difference~"("*degree*C~")")) +
  ylab(expression(Northern~Temperature~"("*degree*C~")")) +
  theme_bw() +
  scale_x_continuous(limits = c(-1.5, 1)) +
  theme(text = element_text(size = 10)) + 
  expand_limits(y = 32.75)

ns_diff <- regional_correlation_data |> ggplot(aes(NS_diff, Northern)) +
  geom_vline(xintercept = 0, color = "darkslategray4", linetype = "dashed", linewidth = 1.2) +
  geom_vline(aes(xintercept = mean(NS_diff)), color = "darkslategray4", linewidth = 1.2) +
  geom_point() +
  annotate("text", x = -1, y = 32.5, label = paste0("South Warmer (", prop_S_warmer_than_N, "% of points)")) +
  annotate("text", x = 0.6, y = 32.5, label = paste0("North Warmer (", prop_N_warmer_than_S, "% of points)")) +
  labs(subtitle = "North - South Comparison") +
  xlab (expression(North~"-"~South~Temperature~Difference~"("*degree*C~")")) +
  ylab(expression(Northern~Temperature~"("*degree*C~")")) +
  theme_bw() +
  scale_x_continuous(limits = c(-1.5, 1)) +
  theme(text = element_text(size = 10)) + 
  expand_limits(y = 32.75)

ws_diff <- regional_correlation_data |> ggplot(aes(WS_diff, Western)) +
  geom_vline(xintercept = 0, color = "darkslategray4", linetype = "dashed", linewidth = 1.2) +
  geom_vline(aes(xintercept = mean(WS_diff)), color = "darkslategray4", linewidth = 1.2) +
  geom_point() +
  annotate("text", x = -1, y = 32.5, label = paste0("South Warmer (", prop_S_warmer_than_W, "% of points)")) +
  annotate("text", x = 0.6, y = 32.5, label = paste0("West Warmer (", prop_W_warmer_than_S, "% of points)")) +
  labs(subtitle = "West - South Comparison") +
  xlab (expression(West~"-"~South~Temperature~Difference~"("*degree*C~")")) +
  ylab(expression(Western~Temperature~"("*degree*C~")")) +
  theme_bw() +
  scale_x_continuous(limits = c(-1.5, 1)) +
  theme(text = element_text(size = 10)) + 
  expand_limits(y = 32.75)

diff_plot <- ggarrange(plotlist = list(nw_diff, ns_diff, ws_diff), common.legend = T, nrow = 3)
ggsave("Figures/InSitu_Regional_Differences.tif", diff_plot, height = 10, width = 11, scale = 0.7, dpi = 300)

saveRDS(diff_plot, "Figures/InSitu_Regional_Differences.rds")
