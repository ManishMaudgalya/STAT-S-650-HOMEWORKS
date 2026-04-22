#Load Data
data <- read.csv("/Users/manishsm/Desktop/Spring 2026/Time Series Analysis/Assignments/quantmod_dataset.csv")
data$Date <- as.Date(data$Date)


#Log Returns
# log return: r_t = log(P_t / P_{t-1})
data$ret_GSPC <- c(NA, diff(log(data$GSPC.Close)))
data$ret_IXIC <- c(NA, diff(log(data$IXIC.Close)))

# Also compute log returns for IXIC OHLC (used in Model 3)
data$ret_IXIC_Open <- c(NA, diff(log(data$IXIC.Open)))
data$ret_IXIC_High <- c(NA, diff(log(data$IXIC.High)))
data$ret_IXIC_Low  <- c(NA, diff(log(data$IXIC.Low)))

#2. Add Lagged Output Variables (ARX feedback terms)
data$ret_GSPC_lag1 <- c(NA, head(data$ret_GSPC, -1))
data$ret_GSPC_lag2 <- c(NA, NA, head(data$ret_GSPC, -2))
data$ret_GSPC_lag3 <- c(NA, NA, NA, head(data$ret_GSPC, -3))
data_clean <- na.omit(data)
cat("Rows after computing log returns and lags:", nrow(data_clean), "\n")
cat("Date range:", as.character(min(data_clean$Date)),
    "to", as.character(max(data_clean$Date)), "\n\n")

# Quick summary of log returns
cat("Log Return Summary\n")
cat(sprintf("GSPC log return: Mean=%.5f, SD=%.5f\n",
            mean(data_clean$ret_GSPC), sd(data_clean$ret_GSPC)))
cat(sprintf("IXIC log return: Mean=%.5f, SD=%.5f\n\n",
            mean(data_clean$ret_IXIC), sd(data_clean$ret_IXIC)))


#  Fit Models white noise / iid error assumed via OLS
model1 <- lm(ret_GSPC ~ ret_IXIC, data = data_clean)
model2 <- lm(ret_GSPC ~ ret_IXIC + ret_GSPC_lag1 + ret_GSPC_lag2,
             data = data_clean)
model3 <- lm(ret_GSPC ~ ret_IXIC_Open + ret_IXIC_High + ret_IXIC_Low,
             data = data_clean)


cat("MODEL 1 SUMMARY (Simple - nested)\n")
cat("  ret_GSPC ~ ret_IXIC\n")
print(summary(model1))


cat("MODEL 2 SUMMARY (Complex - contains Model 1)\n")
cat("  ret_GSPC ~ ret_IXIC + ret_GSPC_lag1 + ret_GSPC_lag2\n")
print(summary(model2))


cat("MODEL 3 SUMMARY (Non-nested)\n")
cat("  ret_GSPC ~ ret_IXIC_Open + ret_IXIC_High + ret_IXIC_Low\n")
print(summary(model3))


#  4: F-Test — Simple (Model 1) vs Complex (Model 2)
# H0: lag1 and lag2 coefficients are jointly zero (Model 1 is sufficient)
# H1: at least one lag term adds explanatory power


cat("F-TEST — Model 1 vs Model 2\n")
ftest <- anova(model1, model2)
print(ftest)
cat("If p < 0.05 -> reject H0 -> Model 2 significantly better than Model 1\n\n")


#  5: AIC and BIC — Compare All Three Models
cat(" 5: AIC vs BIC \n")

aic_vals <- AIC(model1, model2, model3)
bic_vals <- BIC(model1, model2, model3)

comparison_table <- data.frame(
  Model   = c("Model1 (Simple)", "Model2 (Complex)", "Model3 (Non-nested)"),
  Params  = aic_vals$df,
  AIC     = round(aic_vals$AIC, 2),
  BIC     = round(bic_vals$BIC, 2),
  R2      = round(c(summary(model1)$r.squared,
                    summary(model2)$r.squared,
                    summary(model3)$r.squared), 4)
)
print(comparison_table)
cat("\nLower AIC/BIC = better fit penalized for complexity.\n\n")


#  6: Choose Best Model
cat(" 6: BEST MODEL SELECTION\n")
cat("Best by AIC:", comparison_table$Model[which.min(comparison_table$AIC)], "\n")
cat("Best by BIC:", comparison_table$Model[which.min(comparison_table$BIC)], "\n\n")
cat("Justification: The best model balances goodness-of-fit (R2) and\n")
cat("parsimony (AIC/BIC penalty). If Model 2 wins on AIC/BIC AND the\n")
cat("F-test is significant, it is the clear best choice.\n\n")

# Use Model 2 as best model going forward
best_model <- model2


#  7: Residual White Noise Diagnostics

residuals_best <- residuals(best_model)


cat(" 7: WHITE NOISE CHECK ON RESIDUALS\n")
# Ljung-Box test
lb_test <- Box.test(residuals_best, lag = 20, type = "Ljung-Box")
cat("Ljung-Box Test (lag = 20):\n")
print(lb_test)
cat("H0: residuals are white noise\n")
cat("If p < 0.05 -> residuals are NOT white noise -> go to  8\n\n")

# Plots
par(mfrow = c(2, 2))

acf(residuals_best, lag.max = 40,
    main = "ACF of Residuals (Model 2, log returns)",
    xlab = "Lag", ylab = "Autocorrelation")

pacf(residuals_best, lag.max = 40,
     main = "PACF of Residuals (Model 2, log returns)",
     xlab = "Lag", ylab = "Partial Autocorrelation")

spectrum(residuals_best, method = "pgram", taper = 0.1,
         main = "Spectral Density of Residuals",
         xlab = "Frequency", ylab = "Spectral Density", log = "no")

plot(data_clean$Date, residuals_best, type = "l", col = "blue",
     main = "Residuals Over Time",
     xlab = "Date", ylab = "Residual")
abline(h = 0, col = "red", lty = 2)

par(mfrow = c(1, 1))


#  8: Extend to lag p=3 if residuals not white noise


cat(" 8: EXTENDED MODEL WITH lag p = 3\n")


model2_ext <- lm(ret_GSPC ~ ret_IXIC + ret_GSPC_lag1 +
                   ret_GSPC_lag2 + ret_GSPC_lag3,
                 data = data_clean)
print(summary(model2_ext))

residuals_ext <- residuals(model2_ext)
lb_ext <- Box.test(residuals_ext, lag = 20, type = "Ljung-Box")
cat("Ljung-Box Test on Extended Model (lag = 20):\n")
print(lb_ext)

cat("\nAIC/BIC: Model2 (p=2) vs Extended (p=3):\n")
ext_compare <- data.frame(
  Model = c("Model2 (lag p=2)", "Model2_ext (lag p=3)"),
  AIC   = round(c(AIC(best_model), AIC(model2_ext)), 2),
  BIC   = round(c(BIC(best_model), BIC(model2_ext)), 2)
)
print(ext_compare)
cat("\n")


#  9: Interpret Parameters + Stability Check


cat(" 9: PARAMETER INTERPRETATION\n")


coef_tbl <- summary(best_model)$coefficients
print(round(coef_tbl, 6))
cat("(Intercept)    : Average S&P 500 log return when all predictors = 0\n")
cat("ret_IXIC       : Beta — 1% NASDAQ return => this % change in S&P 500\n")
cat("ret_GSPC_lag1  : Momentum/feedback from S&P 500 return yesterday\n")
cat("ret_GSPC_lag2  : Feedback from S&P 500 return two days ago\n\n")

# Stability: for ARX model y_t = a1*y_{t-1} + a2*y_{t-2} + input terms
# Roots of characteristic polynomial z^2 - a1*z - a2 must have |root| > 1
a1 <- coef(best_model)["ret_GSPC_lag1"]
a2 <- coef(best_model)["ret_GSPC_lag2"]

char_roots <- polyroot(c(-a2, -a1, 1))
cat("Characteristic polynomial roots:\n")
print(char_roots)
cat("Moduli:\n")
print(round(Mod(char_roots), 4))
cat("Model is stable if all moduli > 1\n")
cat("Stable?", all(Mod(char_roots) > 1), "\n\n")


#  10: Plot Predictions vs Observed
cat(" 10: PREDICTIONS vs OBSERVED\n")
fitted_vals <- fitted(best_model)
obs_vals    <- data_clean$ret_GSPC

rmse <- sqrt(mean((obs_vals - fitted_vals)^2))
cat(sprintf("RMSE: %.6f\n", rmse))
cat(sprintf("R-squared: %.4f\n\n", summary(best_model)$r.squared))

par(mfrow = c(2, 1))

# Log return comparison
plot(data_clean$Date, obs_vals, type = "l", col = "blue", lwd = 0.8,
     main = "Model 2 — Observed vs Fitted Log Returns (S&P 500)",
     xlab = "Date", ylab = "Log Return")
lines(data_clean$Date, fitted_vals, col = "red", lwd = 1, lty = 2)
legend("topright", legend = c("Observed", "Fitted"),
       col = c("blue", "red"), lty = c(1, 2))

# Scatter: observed vs fitted
plot(fitted_vals, obs_vals, pch = 16, cex = 0.3,
     col = rgb(0.2, 0.4, 0.8, 0.4),
     main = "Observed vs Fitted (scatter)",
     xlab = "Fitted Log Return", ylab = "Observed Log Return")
abline(0, 1, col = "red", lwd = 2)
legend("topleft", legend = "45° line (perfect fit)", col = "red", lwd = 2)

par(mfrow = c(1, 1))