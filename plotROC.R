library(ggplot2)
library(pROC)

plotROC <- function(roc){
  df <- data.frame(Sensitivity = roc$sensitivities,
                   Specificity = roc$specificities)
  df <- df[order(df$Sensitivity),]
  ggplot(df, aes(x = Specificity, y = Sensitivity)) +
    geom_line() +
    xlim(1,0) +
    geom_abline(intercept = 1, slope = 1, linetype="dotted") +
    labs(subtitle = paste0("AUC = ", round(roc$auc, 3))) +
    theme_bw()
}