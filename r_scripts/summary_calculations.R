summary_calculations <- function(column) {
  column <- na.omit(column)  # remove NA values
  n <- length(column)
  mean <- mean(column)
  sd <- sd(column)
  se <- sd / sqrt(n)
  median <- median(column)
  iqr <- IQR(column)
  conf_low <- mean - qt(0.975, df=n-1)*se
  conf_high <- mean + qt(0.975, df=n-1)*se
  list(mean = mean, sd = sd, se = se, median = median, iqr = iqr, conf_low = conf_low, conf_high = conf_high)
}

