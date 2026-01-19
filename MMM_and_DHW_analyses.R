# MMM for Palau (29.2309) pulled from NOAA's CRW site (https://coralreefwatch.noaa.gov/product/vs/data/palau.txt)
# MMM for in-situ data 
# Bleaching threshold for corals is MMM + 1 degree C

## start by getting month means and seeing average difference
## adjust NOAA MMM by average offset to get insitu MMM to calculate DHW with

library(tidyverse)
library(ggpubr)
library(GGally)
library(nlme)
library(cowplot)

source("Table_Outputs_Function.R")
source("Time_Series_Functions.R")

temp_data <- read.csv("Data/Merged_Average_Temperatures.csv", stringsAsFactors = F)

temp_data$Month <- lubridate::month(temp_data$Date)
temp_data$Year <- lubridate::year(temp_data$Date)
temp_data$Year_Month <- paste0(temp_data$Year, "-", temp_data$Month)

month_means <- temp_data |> group_by(Origin, Region, Year, Month) |> 
  summarise(Month_Mean = mean(Mean_Temp))

month_means$Year_Month <- paste0(month_means$Year, "-", month_means$Month)
month_means <- as.data.frame(month_means)
month_means$Date <- as.Date(paste0(month_means$Year_Month, "-01"), format = "%Y-%m-%d")

month_means_comp <- month_means |> 
  filter(Region != "All") |> 
  group_by(Region) |> 
  mutate(time_step = interval(min(Date), Date) %/% months(1)) |> 
  ungroup()

mm_comp <- table_outputs(month_means_comp, Month_Mean, compute_mae = TRUE, mae_interval = 1)
write.csv(mm_comp, "Tables/Month_Means_Output.csv", row.names = F)



 # https://coralreefwatch.noaa.gov/product/5km/methodology.php#dhw
## Website above says "long-term mean SST of the climatologically hottest month of year (often referred to as the Maximum Monthly Mean (MMM) SST climatology"
MMM_by_year <- month_means |> group_by(Origin, Year) |> summarise(Max_Month_Mean = max(Month_Mean))

insitu_mm <- filter(month_means, Origin == "In_Situ")
ssst_mm <- filter(month_means, Origin == "SSST")

mm_wide <- merge(insitu_mm[,c(2,6,5)], ssst_mm[,c(2,6,5)], by = c("Region", "Year_Month"))
colnames(mm_wide) <- c("Region", "Year_Month", "InSitu_Month_Mean", "SSST_Month_Mean")
mm_wide$Month_Mean_Diff <- mm_wide$InSitu_Month_Mean - mm_wide$SSST_Month_Mean

ssst_mmm <- 29.2309
## Since we have gaps in our insitu data, instead of calculating our own MMM we adjust the NOAA MMM by adding the average month mean difference (0.58)
insitu_mmm <- ssst_mmm + mean(mm_wide$Month_Mean_Diff)

ssst_bleaching_threshold <- ssst_mmm + 1
insitu_bleaching_threshold <- insitu_mmm + 1

#Calculate Bleaching stress -----------
temp_data$Over_Bleaching_Threshold <- ifelse(temp_data$Origin == "In_Situ" & temp_data$Mean_Temp >= insitu_bleaching_threshold, 
                                            temp_data$Mean_Temp - insitu_bleaching_threshold, 
                                            ifelse(temp_data$Origin == "SSST" & temp_data$Mean_Temp >= ssst_bleaching_threshold,
                                                   temp_data$Mean_Temp - ssst_bleaching_threshold, 0.00))
temp_data$Over_Bleaching_Threshold <- round(temp_data$Over_Bleaching_Threshold, 2)

# Heat accumulation 
temp_data$Is_Over <- ifelse(temp_data$Over_Bleaching_Threshold > 0, 1, 0)
temp_data$Date <- as.Date(temp_data$Date)
heat_accumulation <- filter(temp_data, Region != "All")
heat_accumulation <- heat_accumulation |> mutate(time_step = 1 + as.integer(Date) - min(as.integer(Date)))
ts_day_heat_acc <- heat_accumulation |> 
  group_by(Region) |> 
  mutate(ts_group = group_ts(Date))


heat_acc <- table_outputs(ts_day_heat_acc, Over_Bleaching_Threshold)
write.csv(heat_acc, "Tables/Over_Bleaching_Threshold_Outputs.csv", row.names = F)

# 2020 Heatwave ------------------
## Just take the "All" region because that's the only way to get the continuous time series in 2020
hw_data <- filter(temp_data, Year == "2020" & Region == "All")

over_bleaching_threshold_plot <- hw_data |> 
  ggplot(aes(Date, Over_Bleaching_Threshold)) + 
  geom_line(aes(color = Origin), linewidth = 1, lineend = "round") +
  scale_color_manual(values = c("tomato4", "darkslategray4")) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y") +
  ylab(expression(Over~Bleaching~Threshold~"("*degree*C~")")) + xlab("") + labs(subtitle = "(a)") +
  theme_bw() + 
  theme(
    legend.position = "none", 
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "lines")
  )

## sliding window of 12 weeks to calculate Degree heating weeks (DHW)
hw_data$Day <- as.numeric(hw_data$Date)
Window <- 1:length(unique(hw_data$Date))
Day <- unique(hw_data$Day)

DHW_12w <- rep(0, length(Window))

## Create a data frame for each origin to bind together later
ssst_DHW_data <- as.data.frame(cbind(Window, Day, DHW_12w))
ssst_DHW_data$Origin <- rep("SSST", nrow(ssst_DHW_data))
insitu_DHW_data <- as.data.frame(cbind(Window, Day, DHW_12w))
insitu_DHW_data$Origin <- rep("In_Situ", nrow(insitu_DHW_data))

## split into ssst and insitu datasets to calcuate DHW with shifting window
ssst_hw_data <- filter(hw_data, Origin == "SSST")
insitu_hw_data <- filter(hw_data, Origin == "In_Situ")

generate_DHW <- function(x, y) {    ### x is the hw dataset by origin, y is the empty DHW dataset by origin to fill
  
  start_date <- min(x$Day)
  end_date <- max(x$Day)
  
  for(i in (start_date + 84):end_date) {
    data <- filter(x, x$Day <= i & x$Day > (i - 84))
    data$TwelveWeek_Accumulation <- rep(NA, nrow(data))
    data$TwelveWeek_Accumulation[1] <- data$Over_Bleaching_Threshold[1]
    for(j in 2:nrow(data)) {
      data$TwelveWeek_Accumulation[j] <- data$TwelveWeek_Accumulation[j-1] + data$Over_Bleaching_Threshold[j]
    }
    
    data$DHW <- data$TwelveWeek_Accumulation/7
    
    last_row <- nrow(data)
    y[y$Day == i, "DHW_12w"] <- data$DHW[last_row]
  }
  return(y)
}


insitu_DHW_data <- generate_DHW(insitu_hw_data, insitu_DHW_data)
ssst_DHW_data <- generate_DHW(ssst_hw_data, ssst_DHW_data)

DHW_data <- rbind(insitu_DHW_data, ssst_DHW_data)
DHW_data$Day <- as.Date(DHW_data$Day)

dhw_plot <- DHW_data |> 
  mutate(Origin = gsub("_", " ", Origin)) |> 
  ggplot(aes(Day, DHW_12w)) +
  geom_line(aes(color = Origin),linewidth = 1, lineend = "round") +
  scale_color_manual(values = c("tomato4", "darkslategray4")) +
  ylab("Degree Heating Weeks (DHW)") + labs(subtitle = "(b)") + xlab("Date") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y") +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "lines")
  )

ggarrange(over_bleaching_threshold_plot, dhw_plot, nrow = 2, common.legend = T, legend = "bottom")
ggsave("Figures/Main_Figs/DHW_and_OverBT_Comparison_2020.tif", height = 12, width = 12, scale = 0.7, dpi = 300)


dhw_wide <- merge(insitu_DHW_data, ssst_DHW_data, by = c("Day", "Window"))
dhw_wide <- dhw_wide[,c(1,2,3,5)]
colnames(dhw_wide) <- c("Day", "Window", "DHW_12w_InSitu", "DHW_12w_SSST")
dhw_lm <- lm(DHW_12w_InSitu ~ DHW_12w_SSST, dhw_wide)
# time step = window?
# need to account for autocorrelation with Newey-West
# also need to adjust the null to a slope of 1

# Because this is a moving window, we need an ARMA(1,1) (autoregressive moving-average order 1)
# estimator for the covariance matrix instead of the AR(1) used by the Newey-West
# covariance matrix estimator for the other analyses
lms_autocor <- as.data.frame(unclass(
  coeftest(dhw_lm, NeweyWest(dhw_lm, order.by = dhw_wide$Window))
))

# Compute difference between new null of intercept = 0, slope = 1
lms_autocor_adj <- lms_autocor |> mutate(df = dhw_lm$df.residual, .after = "t value")
lms_autocor_adj$`t value` <- (lms_autocor_adj$Estimate - c(0, 1)) / lms_autocor_adj$`Std. Error`
lms_autocor_adj$`Pr(>|t|)` <- 2 * pt(abs(lms_autocor_adj$`t value`), lms_autocor_adj$df, lower.tail = FALSE)
dhw_cor <- cor(dhw_wide$DHW_12w_SSST, dhw_wide$DHW_12w_InSitu, method = "pearson")

dhw_comp <- dhw_wide |> 
  ggplot(aes(DHW_12w_SSST, DHW_12w_InSitu)) +
  geom_abline(slope = 1, intercept = 0, color = "grey60") +
  geom_point(size = 1.5) +
  geom_abline(slope = lms_autocor_adj[2,1], intercept = lms_autocor_adj[1,1], 
              color = "dodgerblue1", linewidth = 1) +
  coord_fixed() +
  labs(x = "SSST DHW", y = "In Situ DHW", subtitle = "(c)") + scale_x_continuous(limits = c(0.0, 1.5)) +
  theme_bw() + 
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "lines"))

ggsave("Figures/Main_Figs/DHW_insitu_vs_ssst_C.tif", scale = 0.7, height = 8, width = 8, dpi = 300)

legend <- cowplot::get_legend(dhw_plot)

cowplot::plot_grid(
  over_bleaching_threshold_plot, dhw_comp, dhw_plot + theme(legend.position = "none"), legend, 
  nrow = 2, ncol = 2, align = "hv", axis = "tblr", rel_widths = c(0.7, 0.3)
) + theme(plot.background = element_rect(fill = "white", color = NA))

ggsave("Figures/Main_Figs/Heat_Accumulation_Comp.tif", height = 7, width = 12, scale = 0.7, dpi = 300)


# create stats output table in format to match table outputs 
dhw_stats_outputs <- rownames_to_column(lms_autocor_adj)
dhw_stats_outputs$Region <- "all"
dhw_stats_outputs$cor <- dhw_cor
dhw_stats_outputs$diff <- round(mean(dhw_wide$DHW_12w_SSST - dhw_wide$DHW_12w_InSitu), 3)
dhw_stats_outputs <- dhw_stats_outputs[,c(7,1,2,3,4,5,6,8,9)]
colnames(dhw_stats_outputs) <- c("Region", "parameter", "Estimate", "Std..Error", "t.value", "df", "Pr...t..", 
                                 "cor", "diff")
write.csv(dhw_stats_outputs, "Tables/DHW_Outputs.csv", row.names = F)
