library(ggplot2)

# simple scatter plot with regression line
simpleScatter <- function(x, y, title = "", pointSize = 0.5){
  df <- cbind.data.frame(x = x, y = y)
  lm <- lm(y ~ x, df)
  rsq <- round(summary(lm)$r.squared, 4)
  m <- round(lm$coefficients[2], 3)
  b <- round(lm$coefficients[1], 3)
  ggplot(df, aes(x = x, y = y)) +
    geom_point(size = pointSize) +
    theme_bw() +
    geom_abline(slope = m, intercept = b) +
    labs(title = title,
         subtitle = paste0("Rsq = ", rsq, "\n",
                           "y = ", m, "x + ", b))
}