library(tidyverse)
library(lubridate)
library(GGally)
library(ggpubr)
library(fishualize)
library(cowplot)
library(nlme)
library(ggtext)

source("Table_Outputs_Function.R")
source("Time_Series_Functions.R")

## Read in insitu data then ssst data and fix classes
insitu_data <- read.csv("Data/Average_Nighttime_InSitu_Data.csv", stringsAsFactors = F)
insitu_data$date <- as.Date(insitu_data$date)
insitu_data$Origin <- rep("In_Situ", nrow(insitu_data))
insitu_data$Variation <- insitu_data$max_night_temp - insitu_data$min_night_temp
colnames(insitu_data) <- c("Date", "Mean_Temp", "Min_Temp", "Max_Temp", "SD", "Region", "Origin", "Variation")

ssst_data <- read.csv("Data/Average_SSST_Data.csv", stringsAsFactors = F)
ssst_data$date <- as.Date(ssst_data$date)
ssst_data$Origin <- rep("SSST", nrow(ssst_data))
ssst_data$Variation <- ssst_data$max_temp - ssst_data$min_temp
colnames(ssst_data) <- c("Date", "Region", "Mean_Temp", "SD", "Min_Temp", "Max_Temp", "Origin", "Variation")

avg_temps <- rbind(insitu_data, ssst_data)
write.csv(avg_temps, "Data/Merged_Average_Temperatures.csv", row.names = F)

## Insitu vs ssst ----------

avg_temps$Region <- factor(avg_temps$Region, ordered = T, levels = c("Southern", "Western", "Northern", "All"))
avg_temps$Origin <- as.factor(avg_temps$Origin)

ts_gaps_temp_data <- avg_temps |> 
  group_by(Region) |> 
  mutate(ts_group = group_ts(Date))

#region_temp_data <- filter(ts_gaps_temp_data, Region != "All")
region_temp_data <- ts_gaps_temp_data

region_temp_data |> 
  ggplot(aes(Date, Mean_Temp)) + 
  geom_line(aes(group = ts_group, color = Origin), linewidth = 1.5, alpha = 0.5) + 
  scale_color_manual(values = c("darkslategray4", "tomato4")) +
  ylab(expression(Mean~Temperature~"("*degree*C~")")) +
  theme_bw() +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"), text = element_text(size = 10),
        legend.position = "bottom") #+
  #facet_wrap(~Region, ncol = 1, scale = "free_x)
ggsave("Figures/Main_Figs/InSitu_vs_SSST.tif", height = 6, width = 10, scale = 0.7, dpi = 300)


## Correlation tests, means, and ranges 
region_temps <- filter(avg_temps, Region != "All")
region_temps$Region <- as.character(region_temps$Region)
region_temps <- region_temps |> mutate(time_step = 1 + as.integer(Date) - min(as.integer(Date)))
region_temps <- region_temps[,c(1,6,7,2,3,4,5,8,9)]
mean_corr <- table_outputs(region_temps, Mean_Temp, compute_mae = TRUE)
write.csv(mean_corr, "Tables/Mean_Temp_Outputs.csv", row.names = FALSE)


# Differences in temps ----------
wide_ssst_insitu_data <- merge(insitu_data[,c(1:4,6,8)], ssst_data[,c(1:3,5,6,8)], by = c("Date", "Region"))
colnames(wide_ssst_insitu_data) <- c("Date", "Region", "InSitu_Mean_Temp", "InSitu_Min_Temp", "InSitu_Max_Temp", "InSitu_Variation", 
                                     "SSST_Mean_Temp", "SSST_Min_Temp", "SSST_Max_Temp", "SSST_Variation")
wide_ssst_insitu_data$Mean_Temp_Diff <- wide_ssst_insitu_data$InSitu_Mean_Temp - wide_ssst_insitu_data$SSST_Mean_Temp

wide_ssst_insitu_data |> filter(Region != "All") |> 
  ggplot(aes(InSitu_Mean_Temp, SSST_Mean_Temp)) +
  geom_point() +
  geom_smooth(method = "lm", color = "skyblue4", linewidth = 1.3) +
  facet_wrap(~Region) +
  stat_cor(method = "pearson", aes(label = after_stat(r.label)), color = "skyblue4") +
  stat_regline_equation(label.x = 29, label.y = 27, color = "skyblue4") +
  geom_abline(intercept = 0, slope = 1, color = "gray57") +
  theme_bw() +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 12), text = element_text(size = 12))
ggsave("Figures/InSitu_SST_Regional_Correlations.tif", height = 8, width = 12, scale = 0.7, dpi = 300)

## Create column for seasons and plot 
wide_ssst_insitu_data$Month <- lubridate::month(wide_ssst_insitu_data$Date)
wide_ssst_insitu_data$Season <- ifelse(wide_ssst_insitu_data$Month %in% c(01:03), "Winter", 
                        ifelse(wide_ssst_insitu_data$Month %in% c(04:06), "Spring", 
                               ifelse(wide_ssst_insitu_data$Month %in% c(07:09), "Summer", 
                                      ifelse(wide_ssst_insitu_data$Month %in% c(10:12), "Fall", NA))))

date_breaks <- seq(as.Date("2017-11-01"), as.Date("2020-12-01"), by = "2 months")

date_labs <- ifelse(
  month(date_breaks) == 1, 
  paste(as.character(month(date_breaks, label = TRUE)), paste0("**", year(date_breaks), "**"), sep = "<br>"), 
  as.character(month(date_breaks, label = TRUE))
)

wide_ssst_insitu_data |> 
  group_by(Region) |> 
  mutate(
    Diff_Smoothed = gaussian_smooth(Mean_Temp_Diff, Date, bw = 30)$x, 
    Block = group_ts(Date)
  ) |> 
  ungroup() |> 
  ggplot(aes(Date, Mean_Temp_Diff)) +
  geom_bar(stat = "identity", aes(fill = Season), width = 1) +
  scale_fill_manual(values = c("orangered3", "darkgreen", "goldenrod3", "dodgerblue3")) +
  ylab("Difference in Mean Temp (In Situ - SSST)") +
  scale_x_date(breaks = date_breaks, labels = date_labs) + 
  geom_hline(yintercept = 0, color = "black") +
  geom_line(aes(y = Diff_Smoothed, group = Block), lineend = "round", color = "white", alpha = 0.6, linewidth = 1.5) +
  geom_line(aes(y = Diff_Smoothed, group = Block), lineend = "round") +
  facet_wrap(~Region, ncol = 1) +
  theme_bw() +
  theme(
    text = element_text(size = 12), 
    axis.text.x = element_markdown(size = 10),
    strip.background = element_blank(), 
    strip.text = element_text(size = 12, hjust = 0)
  ) 

ggsave("Figures/Main_Figs/All_Regions_Diff_Over_Time.tif", height = 9, width = 12, scale = 0.7, dpi = 300)

## get percent of days when ssst reads higher than in situ data
all <- filter(wide_ssst_insitu_data, Region == "All")
percent_ssst_higher <- nrow(filter(all, Mean_Temp_Diff < 0))/nrow(all)

# Nightly Max temp comparison --------------

wide_ssst_insitu_data |> filter(Region != "All") |> 
  ggplot(aes(x = Date)) +
  geom_line(aes(y = SSST_Max_Temp, color = "SSST"), linewidth = 1.25) +
  geom_line(aes(y = InSitu_Max_Temp, color = "InSitu"), linewidth = 1.25) +
  scale_color_manual(values = c("tomato4", "darkslategray4")) +
  ylab(expression(Maximum~Temperature~"("*degree*C~")")) +
  facet_wrap(~Region, nrow = 1, scales = "free_x") +
  theme_bw() +
  theme(legend.title = element_blank(), strip.background = element_blank())
ggsave("Figures/Max_Temp_Comparison.tif", height = 6, width = 12, scale = 0.7)

## Max correlation
max_corr <- table_outputs(region_temps, Max_Temp, compute_mae = TRUE)
write.csv(max_corr, "Tables/Max_Temp_Outputs.csv", row.names = FALSE)

wide_ssst_insitu_data |> 
  ggplot(aes(InSitu_Max_Temp, SSST_Max_Temp)) +
  geom_point() +
  geom_smooth(method = "lm", color = "salmon4", linewidth = 1.3) +
  facet_wrap(~Region) +
  stat_cor(method = "pearson", aes(label = after_stat(r.label)), color = "salmon4") +
  stat_regline_equation(label.x = 29, label.y = 27, color = "salmon4") +
  geom_abline(intercept = 0, slope = 1, color = "gray57") +
  labs(main = "Maximum Temperatures") +
  theme_bw()
ggsave("Figures/Max_Temp_Correlation.tif", height = 8, width = 12, scale = 0.7, dpi = 300)


# Nightly Min temp comparison --------------

wide_ssst_insitu_data |> filter(Region != "All") |> 
  ggplot(aes(x = Date)) +
  geom_line(aes(y = SSST_Min_Temp, color = "SSST"), linewidth = 1.25) +
  geom_line(aes(y = InSitu_Min_Temp, color = "InSitu"), linewidth = 1.25) +
  scale_color_manual(values = c("tomato4", "darkslategray4")) +
  ylab(expression(Minimum~Temperature~"("*degree*C~")")) +
  facet_wrap(~Region, nrow = 1, scales = "free_x") +
  theme_bw() +
  theme(legend.title = element_blank(), strip.background = element_blank())
ggsave("Figures/Min_Temp_Comparison.tif", height = 6, width = 12, scale = 0.7, dpi = 300)

## Min correlation
min_corr <- table_outputs(region_temps, Min_Temp, compute_mae = TRUE)
write.csv(min_corr, "Tables/Min_Temp_Outputs.csv", row.names = FALSE)

wide_ssst_insitu_data |>  
  ggplot(aes(InSitu_Min_Temp, SSST_Min_Temp)) +
  geom_point() +
  geom_smooth(method = "lm", color = "steelblue4", linewidth = 1.3) +
  facet_wrap(~Region) +
  stat_cor(method = "pearson", aes(label = after_stat(r.label)), color = "steelblue4") +
  stat_regline_equation(label.x = 29, label.y = 27, color = "steelblue4") +
  geom_abline(intercept = 0, slope = 1, color = "gray57") +
  labs(main = "Minimum Temperatures") +
  theme_bw()
ggsave("Figures/Min_Temp_Correlation.tif", height = 8, width = 12, scale = 0.7, dpi = 300)
