##6/12/2023
# Table output function

#' Create table of correlation and regression values by region
#' @param x data frame with time series of thermal metric in long format. Requires "Origin" column
#' with "SSST" and "In_Situ" values, Region with "Northern", "Western", and "Southern" values (NOT "All"),
#' and integer time step column with number of days since first observation.  
#' @param var  name of the variable to buil comparisons for (without quotes)
#' @param print if false (default), return dataframe of values. Otherwise, print summary to console.
#' @param compute_mae whether to compute mean absolute error on a rolling window
#' @param mae_interval number of time steps to use for MAE rolling window (e.g., 30 for daily data, 1 for monthly data)
#'
#' @return Console output or data frame
table_outputs <- function(x, var, print = FALSE, compute_mae = FALSE, mae_interval = 30) {
  
  require("tidyverse")
  require("lmtest")
  require("sandwich")
  
  var <- enquo(var)
  
  # Function to print correlation, mean difference, and range of difference
  print_summary <- function(y1, y2) {
    cat(paste("correlation:", round(cor(y1, y2), 3)), "\n")
    mean_difference <- y2 - y1
    cat(paste("mean:", round(mean(mean_difference), 3)), "\n")
    cat(paste("range:", paste(round(range(mean_difference), 3), collapse = ", ")), "\n")
  }
  
  # Wide data frame for all regions
  data <- x |>
    dplyr::select(Origin, Region, time_step, !!var) |> 
   # mutate(Origin = gsub("-", "_", Origin)) |> 
    pivot_wider(names_from = "Origin", values_from = quo_name(var))
  
  # Data frames for each region
  data <- append(list(all = data), split(data, data$Region))
  
  # Fit linear models
  models <- lapply(data, function(x) lm(In_Situ ~ SSST, data = x))
  
  # Apply autocorrelation corrections and extract coefficient tables
  lm_tabs <- append(
    list(all = as.data.frame(coeftest(models$all, vcovPC(models$all, cluster = data$all$Region, order.by = data$all$time_step, pairwise = TRUE))[,])),
    purrr::map2(models[2:4], data[2:4], function(x, y) {
      as.data.frame(coeftest(x, NeweyWest(x, order.by = y$time_step))[,])
    }))
  
  # Change slope null to 1, adjust for multiple comparisons
  lm_tabs <- purrr::map2(lm_tabs, models, function(x, y) {
    
    x <- x |> mutate(df = y$df.residual, .after = "t value")
    x$`t value` <- (x$Estimate - c(0, 1)) / x$`Std. Error`
    x$`Pr(>|t|)` <- 2 * pt(abs(x$`t value`), x$df, lower.tail = FALSE)
    x$`Pr(>|t|)` <- p.adjust(x$`Pr(>|t|)`, method = "bonferroni", n = length(lm_tabs))
    x <- round(x, 3)
    x$`Pr(>|t|)` <- as.character(ifelse(x$`Pr(>|t|)` < 0.001, "< 0.001", x$`Pr(>|t|)`))
    x
    
  })
  
  # Compute mean absolute error (MAE)
  # Measured here by making 1-month-ahead predictions starting after the second time block
  if (isTRUE(compute_mae)) {
    
    data <- lapply(data, \(x) mutate(x, block = group_ts(time_step)))
    mae <- setNames(vector("list", length(data)), names(data))
    
    for (i in seq_along(data)) {
      
      for (b in 3:max(data[[i]]$block)) {
        
        data_fit <- data[[i]][data[[i]]$block < b,]
        data_test <- data[[i]][data[[i]]$block == b,]
        
        months_ib <- as.integer(cut(
          data_test$time_step, 
          breaks = seq(min(data_test$time_step), max(data_test$time_step) + mae_interval, by = mae_interval), 
          include.lowest = TRUE
        ))
        
        data_fit <- rbind(data_fit, data_test[months_ib == 1,])
        
        for (m in 2:max(months_ib)) {
          
          fit <- lm(In_Situ ~ SSST, data = data_fit)
          
          obs <- data_test$In_Situ[months_ib == m]
          fcst <- predict(fit, newdata = data_test[months_ib == m,])
          
          mae[[i]] <- c(mae[[i]], mean(abs(obs - fcst)))
          
          data_fit <- rbind(data_fit, data_test[months_ib == m,])
          
        }
        
      } 
      
    }
    
    mae <- round(vapply(mae, mean, numeric(1)), 3)
    
  }
  
  
  if (print) {
    
    for (i in 1:length(data)) {
      cat(paste0("\n", str_pad(paste0(names(data)[i], " "), 50, "right", pad = "-"), "\n\n"))
      with(data[[i]], print_summary(In_Situ,SSST))
      cat(paste0("\nLinear model:\n\n"))
      print(lm_tabs[[names(data)[i]]])
    }
    
    cat("\nTests are against null of intercept = 0, slope = 1")
    cat(paste("\nBonferroni correction applied for a family of 4 tests", length(lm_tabs), "tests"))
    cat(paste("\nNewey-West / Beck & Katz correction applied to regional / full SEs"))
    
  } else {
    
    lm_tabs <- lm_tabs |> lapply(rownames_to_column, "parameter") |> bind_rows(.id = "Region")
    
    cors <- vapply(data, \(x) round(cor(x$SSST, x$In_Situ), 3), numeric(1))
    lm_tabs$cor <- cors[lm_tabs$Region]
    
    diffs <- vapply(data, \(x) round(mean(x$SSST - x$In_Situ), 3), numeric(1))
    lm_tabs$diff <- diffs[lm_tabs$Region]
    
    rmse <- vapply(models, \(x) round(sqrt(mean(resid(x)^2)), 3), numeric(1))
    lm_tabs$rmse <- rmse[lm_tabs$Region]
    
    # Redundant with lm_tabs$cor^2 but storing anyway to avoid rounding error
    rsq <- vapply(models, \(x) round(summary(x)$r.squared, 3), numeric(1))
    lm_tabs$rsq <- rsq[lm_tabs$Region]
    
    if (isTRUE(compute_mae)) {
      lm_tabs$mae <- mae[lm_tabs$Region]
    } else {
      lm_tabs$mae <- NA
    }
    
    return(lm_tabs)
    
  }
  
}
