library(tidyverse)
library(ggpubr)
library(GGally)

region_sst <- read.csv("Data/Regional_CRW_SST_Full.csv", header = T, stringsAsFactors = F)
region_sst$layer <- gsub("X", "", region_sst$layer)
region_sst$layer <- as.POSIXct(region_sst$layer, format = "%Y.%m.%d.%H.%M.%OS")
colnames(region_sst) <- c("long", "lat", "date_time", "sst", "Region")
region_sst$date <- as.Date(format(region_sst$date_time, "%Y-%m-%d"))
region_sst$Region <- factor(region_sst$Region, ordered = T, levels = c("Southern", "Western", "Northern"))
region_sst$data_type <- rep("Satellite", nrow(region_sst))
write.csv(region_sst, "Data/Regional_CRW_SST_Full_Clean.csv", row.names = F)

## Read in insitu data so we can filter the sst data to match date-wise
insitu_data <- read.csv("Data/Average_Nighttime_InSitu_Data.csv", stringsAsFactors = F)
insitu_data$date_match <- paste(insitu_data$date, insitu_data$Region)

region_sst$date_match <- paste(region_sst$date, region_sst$Region)

region_sst <- filter(region_sst, date_match %in% insitu_data$date_match)

# Plot mean and st dev by region
avg_region_sst <- region_sst |> group_by(date, Region) |>  
  summarise(mean = mean(sst), 
            stdev = sd(sst), 
            min_temp = min(sst),
            max_temp = max(sst),
            .groups = "drop") 

avg_total_sst <- region_sst |> group_by(date) |>  
  summarise(mean = mean(sst), 
            stdev = sd(sst), 
            min_temp = min(sst),
            max_temp = max(sst),
            .groups = "drop") 
avg_total_sst$Region <- rep("All", nrow(avg_total_sst))

full_avg_region_sst <- rbind(avg_region_sst, avg_total_sst)
write.csv(full_avg_region_sst, "Data/Average_SSST_Data.csv", row.names = F)

avg_region_sst |> 
    ggplot(aes(date, mean)) +
      geom_ribbon(aes(ymin = mean - stdev, ymax = mean + stdev)) +  
      geom_line(aes(color = Region), linewidth = 1.25, alpha = 0.9) +
      scale_colour_manual(values = c("darkslategray2", "darkslategray4", "darkslategrey")) +
      theme(panel.background= element_rect(fill = "gray95"), panel.grid = element_blank()) +
      labs(title = "Mean SST") + theme_bw()


# Regional Comparison ----------

## Because our regions span different dates and have a different amount of days with data, 
## we need to filter the data so all regions cover the same dates
North_dates <- avg_region_sst$date[avg_region_sst$Region == "Northern"]
West_dates <- avg_region_sst$date[avg_region_sst$Region == "Western"]
South_dates <- avg_region_sst$date[avg_region_sst$Region == "Southern"]
regional_comparison_data <- avg_region_sst[avg_region_sst$date %in% North_dates & 
                                                    avg_region_sst$date %in% West_dates &
                                                    avg_region_sst$date %in% South_dates,]

regional_correlation_data <- regional_comparison_data[,c(1:3)] |> pivot_wider(names_from = Region, values_from = mean)
sst_regional_correlation <- cor(regional_correlation_data[,-1])
write.table(sst_regional_correlation, "Tables/SSST_Regional_Correlations_Outputs.csv", sep = ",", row.names = T, col.names = T)

ggpairs(regional_correlation_data[,-1], title = "SSST Regional Correlations")
ggsave("Figures/SSST_Regional_Correlations.tif", width = 10, height = 8, scale = 0.7, dpi = 300)


# Plot regional comparisons --------
  ## Comparing nighttime average temps - correlation plots
  
  nw <- regional_correlation_data |> ggplot(aes(Northern, Western)) +
  geom_point() +
  geom_smooth(method = "lm", color = "mediumpurple") +
  stat_cor(method = "pearson", aes(label = after_stat(r.label)), color = "mediumpurple") +
  stat_regline_equation(label.y = 31, color = "mediumpurple") +
  geom_abline(intercept = 0, slope = 1, color = "gray57") +
  labs(subtitle = "North-West Correlation") +
  theme_bw()

ns <- regional_correlation_data |> ggplot(aes(Northern, Southern)) +
  geom_point() +
  geom_smooth(method = "lm", color = "mediumpurple") +
  stat_cor(method = "pearson", aes(label = after_stat(r.label)), color = "mediumpurple") +
  stat_regline_equation(label.y = 31, color = "mediumpurple") +
  geom_abline(intercept = 0, slope = 1, color = "gray57") +
  labs(subtitle = "North-South Correlation") +
  theme_bw()

sw <- regional_correlation_data |> ggplot(aes(Southern, Western)) +
  geom_point() +
  geom_smooth(method = "lm", color = "mediumpurple") +
  stat_cor(method = "pearson", aes(label = after_stat(r.label)), color = "mediumpurple") +
  stat_regline_equation(label.y = 31, color = "mediumpurple") +
  geom_abline(intercept = 0, slope = 1, color = "gray57") +
  labs(subtitle = "South-West Correlation") +
  theme_bw()

plotlist <- list(nw, ns, sw) 
ggarrange(plotlist = plotlist, common.legend = T, nrow = 1, legend = "bottom")
ggsave("Figures/Regional_SSST_Correlation.tif", height = 6, width = 12, scale = 0.7, dpi = 300)

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
  geom_vline(xintercept = 0, color = "tomato4", linetype = "dashed", linewidth = 1.2) +
  geom_vline(aes(xintercept = mean(NW_diff)), color = "tomato4", linewidth = 1.2) +
  geom_point() +
  annotate("text", x = -1, y = 32.5, label = paste0("West Warmer (", prop_W_warmer_than_N, "% of points)")) +
  annotate("text", x = 0.6, y = 32.5, label = paste0("North Warmer (", prop_N_warmer_than_W, "% of points)")) +
  labs(title = "(b)", subtitle = "North - West Comparison") + 
  xlab (expression(North~"-"~West~Temperature~Difference~"("*degree*C~")")) +
  ylab(expression(Northern~Temperature~"("*degree*C~")")) +
  theme_bw() +
  scale_x_continuous(limits = c(-1.5, 1)) +
  theme(text = element_text(size = 10)) + 
  expand_limits(y = 32.75)

ns_diff <- regional_correlation_data |> ggplot(aes(NS_diff, Northern)) +
  geom_vline(xintercept = 0, color = "tomato4", linetype = "dashed", linewidth = 1.2) +
  geom_vline(aes(xintercept = mean(NS_diff)), color = "tomato4", linewidth = 1.2) +
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
  geom_vline(xintercept = 0, color = "tomato4", linetype = "dashed", linewidth = 1.2) +
  geom_vline(aes(xintercept = mean(WS_diff)), color = "tomato4", linewidth = 1.2) +
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

ssst_diff_plot <- ggarrange(plotlist = list(nw_diff, ns_diff, ws_diff), common.legend = T, nrow = 3)
ggsave("Figures/SSST_Regional_Differences.tif", height = 10, width = 11, scale = 0.7, dpi = 300)

insitu_diff_plot <- readRDS("Figures/InSitu_Regional_Differences.rds")

diff_plot_joined <- ggarrange(insitu_diff_plot, ssst_diff_plot, nrow = 1)
ggsave("Figures/Main_Figs/Regional_Differences.tif", diff_plot_joined, 
       height = 11, width = 18, scale = 0.7, units = "in", dpi = 300)
