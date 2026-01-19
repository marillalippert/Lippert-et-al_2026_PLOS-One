## Make a plot of correlations for all metrics to summarize.

library(tidyverse)
library(ggpubr)
library(cowplot)

# Read in and combine outputs from mean, min, max, and month mean comparisons
mean_cor_data <- read.csv("Tables/Mean_Temp_Outputs.csv", stringsAsFactors = F)
mean_cor_data$Metric <- rep("Mean_Temp", nrow(mean_cor_data))

max_cor_data <- read.csv("Tables/Max_Temp_Outputs.csv", stringsAsFactors = F)
max_cor_data$Metric <- rep("Max_Temp", nrow(max_cor_data))

min_cor_data <- read.csv("Tables/Min_Temp_Outputs.csv", stringsAsFactors = F)
min_cor_data$Metric <- rep("Min_Temp", nrow(min_cor_data))

month_mean_cor_data <- read.csv("Tables/Month_Means_Output.csv", stringsAsFactors = F)
month_mean_cor_data$Metric <- rep("Month_Mean", nrow(month_mean_cor_data))

over_bt_cor_data <- read.csv("Tables/Over_Bleaching_Threshold_Outputs.csv", stringsAsFactors = F)
over_bt_cor_data$Metric <- rep("Over_BT", nrow(over_bt_cor_data))

var_cor_data <- read.csv("Tables/Diurnal_Variation_vs_SSST_Mean_Outputs.csv", stringsAsFactors = F)
var_cor_data$Metric <- rep("Diurnal_Var", nrow(over_bt_cor_data))

dhw_cor_data <- read.csv("Tables/DHW_Outputs.csv", stringsAsFactors = F)
dhw_cor_data$Metric <- rep("DHW", nrow(dhw_cor_data))
dhw_cor_data$rmse <- dhw_cor_data$rsq <- dhw_cor_data$mae <- NA

daily_temp_data <- read.csv("Tables/Daily_Mean_Outputs.csv", stringsAsFactors = F)
daily_temp_data$Metric <- rep("Daily_Avg", nrow(daily_temp_data))

corr_data <- rbind(mean_cor_data, max_cor_data, min_cor_data, month_mean_cor_data, over_bt_cor_data, var_cor_data, dhw_cor_data, daily_temp_data)
colnames(corr_data) <- c("Region","Parameter", "Estimate", "Standard_Error", "T-value", "DoF", "P-value", 
                         "Correlation", "Mean_Difference(SSST-InSitu)", "RMSE", "R_sq", "MAE", "Metric")
write.csv(corr_data, "Tables/Compiled_Correlations_Table.csv", row.names = F)


corr_data$Metric <- as.factor(corr_data$Metric)
corr_data$Region <- as.factor(corr_data$Region)
#corr_data$Significance <- ifelse(corr_data$`P-value` < 0.05, "Significant (<0.05)", "Not Significant")
plot_corr_data <- pivot_wider(corr_data, names_from = "Parameter", 
                 values_from = c("Estimate", "Standard_Error", "P-value", "T-value"))

plot_corr_data$Intercept_Significance <- ifelse(plot_corr_data$`P-value_(Intercept)` < 0.05, 
                                                "Significant (<0.05)", "Not Significant")
plot_corr_data$Slope_Significance <- ifelse(plot_corr_data$`P-value_SSST` < 0.05, 
                                                "Significant (<0.05)", "Not Significant")



# Upper left quadrant = In-situ more variable, SST biased low
# Upper right = In situ more variable, SST biased high
# Lower left = SST more variable, SST biased low
# Lower right = SST more variable, SST biased high
plot_corr_data <- filter(plot_corr_data, Metric != "Diurnal_Var" & Metric != "Over_BT" & Metric != "DHW" & Metric != "Daily_Avg")

text_df <- data.frame(
  x = c(-3, -3, 2.5, 2.5, -3, -3, 2.5, 2.5), 
  y = c(1.25, 1.22, 1.25, 1.22, 0.92, 0.89, 0.92, 0.89),
  label = c("In Situ more variable", "SSST biased high", "In Situ more variable", "SSST biased low", 
            "SSST more variable", "SSST biased high", "SSST more variable", "SSST biased low"),
  Metric = "Mean_Temp"
)

# create a separate dataframe based on significant estimates (slope or intercept) with coords for where the * marker should go
significance <- data.frame(
  x = c(-0.488, -0.516),
  y = c(1.053, 1.075),
  label = "*",
  Metric = c("Mean_Temp", "Max_Temp")
)

facet_labels <- c('Mean_Temp' = "Nightly Average Temperature", 
                  'Month_Mean' = "Monthly Mean",
                  'Max_Temp' = "Nightly Maximum Temperature", 
                  'Min_Temp' = "Nightly Minimum Temperature")

coef_plot <- plot_corr_data |> 
  ggplot(aes(x = `Estimate_(Intercept)`, y = Estimate_SSST)) +
  geom_hline(yintercept = 1, linetype = "dashed") + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_linerange(aes(
    xmin = qnorm(0.025, mean = `Estimate_(Intercept)`, sd = `Standard_Error_(Intercept)`), 
    xmax = qnorm(0.975, mean = `Estimate_(Intercept)`, sd = `Standard_Error_(Intercept)`),
    color = Region
  ), linewidth = 0.7, lineend = "round") +
  geom_linerange(aes(
    ymin = qnorm(0.025, mean = Estimate_SSST, sd = Standard_Error_SSST), 
    ymax = qnorm(0.975, mean = Estimate_SSST, sd = Standard_Error_SSST), 
    color = Region
  ), linewidth = 0.7, lineend = "round") +
  geom_point(aes(color = Region), size = 3) +
  scale_linetype_discrete(guide = "none") + 
  geom_text(aes(x = x, y = y, label = label), data = text_df, size = 2.5, fontface = "bold") + 
  geom_text(aes(x = x, y = y * 1.025, label = label), data = significance, size = 8, colour = "mediumpurple4") + 
  facet_wrap(~factor(Metric, c("Mean_Temp", "Month_Mean", "Max_Temp", "Min_Temp")), 
             labeller = as_labeller(facet_labels)) +
  scale_color_manual(values = c("mediumpurple4", "goldenrod3", "darkgreen", "magenta4")) +
  coord_cartesian(clip = "off") + 
  theme_bw() + 
  theme(strip.background = element_blank(), legend.position = "bottom") + 
  scale_y_continuous(limits = c(0.89, 1.4)) +
  labs(x = "Intercept", y = "Slope")

ggsave("Figures/Main_Figs/All_Metrics_Coefs.tif", coef_plot, height = 8, width = 12, dpi = 300, scale = 0.7)

stat_plot <- plot_corr_data |> 
  mutate(diff = - `Mean_Difference(SSST-InSitu)`) |> 
  dplyr::select(Region, Metric, RMSE, R_sq, MAE, diff) |>
  pivot_longer(c(RMSE, R_sq, MAE, diff), names_to = "Statistic") |> 
  mutate(
    Statistic = factor(case_when(
      Statistic == "MAE" ~ "One-Month~Mean~Absolute~Error~(degree*C)", 
      Statistic == "RMSE" ~ "Model~Root~Mean~Square~Error~(degree*C)", 
      Statistic == "R_sq" ~ "Model~R^2", 
      Statistic == "diff" ~ "Mean~In~Situ~-~SSST~difference~(degree*C)"
    ), 
    levels = c(
      "Mean~In~Situ~-~SSST~difference~(degree*C)", 
      "Model~R^2", 
      "Model~Root~Mean~Square~Error~(degree*C)", 
      "One-Month~Mean~Absolute~Error~(degree*C)"
    ))
  ) |> 
  ggplot(aes(Metric, value, color = Region)) +
  geom_segment(aes(yend = 0), position = position_dodge(width = 0.6), linewidth = 0.8, alpha = 0.5) +
  geom_point(aes(group = Region), position = position_dodge(width = 0.6), size = 2.8, color = "white") + 
  geom_point(position = position_dodge(width = 0.6), size = 1.8) + 
  scale_x_discrete(labels = gsub(" Temperature", "", facet_labels)) +
  scale_y_continuous(expand = c(0, 0), position = "right", labels = function(x) ifelse(x == 0, "0", x)) + 
  geom_blank(aes(y = value * 1.1)) +
  facet_wrap(~Statistic, nrow = 2, scales = "free_x", labeller = "label_parsed") + 
  scale_color_manual(values = c("mediumpurple4", "goldenrod3", "darkgreen", "magenta4")) +
  coord_flip(clip = "off") +
  theme_bw() + 
  theme(
    strip.background = element_blank(), 
    legend.position = "bottom", 
    axis.title = element_blank(),
    strip.placement = "outside", 
    panel.spacing.y = unit(1, "lines"), 
    panel.spacing.x = unit(0.5, "lines"), 
    panel.grid.major.y = element_blank()
  ) 

ggsave("Figures/Main_Figs/All_Metrics_Performance.tif", stat_plot, height = 8, width = 12, dpi = 300, scale = 0.7)

