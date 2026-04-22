#3. Plot the time series two of your variables.
#Detrend the time series by fitting a polynomial of order 2 (or more if needed)
#and then plot the residuals.
library(quantmod)

# Two variables: S&P 500 (^GSPC) and NASDAQ Composite (^IXIC)
symbols <- c("^GSPC", "^IXIC")
getSymbols(symbols, src = "yahoo", from = "2010-01-01",
to = "2026-01-01", auto.assign = TRUE, warnings = FALSE)

sp500 <- Ad(GSPC)
nasdaq <- Ad(IXIC)
data <- na.omit(merge(sp500, nasdaq))
colnames(data) <- c("SP500", "NASDAQ")

#original time series
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
plot(index(data), data$SP500, type = "l",
main = "SP500 Adjusted Close", 
xlab = "Date", 
ylab = "Level")
plot(index(data), data$NASDAQ, type = "l",
main = "NASDAQ Adjusted Close",
xlab = "Date",
ylab = "Level")

# Detrend with 2nd degree polynomial
t <- seq_len(NROW(data))

fit_sp500 <- lm(as.numeric(data$SP500) ~ poly(t, 2, raw = TRUE))
fit_nasdaq <- lm(as.numeric(data$NASDAQ) ~ poly(t, 2, raw = TRUE))

resid_sp500 <- residuals(fit_sp500)
resid_nasdaq <- residuals(fit_nasdaq)

resid_xts <- xts(cbind(resid_sp500, resid_nasdaq), order.by = index(data))
colnames(resid_xts) <- c("SP500_resid", "NASDAQ_resid")

# Plot residuals (detrended series)
plot(index(resid_xts), resid_xts$SP500_resid, type = "l",
main = "SP500 Residuals 2nd deg poly detrend",
xlab = "Date",
ylab = "Residual")
abline(h = 0, lty = 2)

plot(index(resid_xts), resid_xts$NASDAQ_resid, type = "l",
main = "NASDAQ Residuals 2nd deg poly detrend",
xlab = "Date",
ylab = "Residual")
abline(h = 0, lty = 2)

#4. Compute a moving average estimate of the variance for each of the residuals of the two variables. Does the variance look stationary? If not, difference the series once and see if the difference is stationary.
window_size <- 126  # 6 months of trading days

# Rolling moving-average estimate of variance: E[e_t^2] over a moving window
rollvar_sp500 <- zoo::rollapply(resid_xts$SP500_resid^2,
width = window_size,FUN = mean,align = "right",fill = NA)

rollvar_nasdaq <- zoo::rollapply(resid_xts$NASDAQ_resid^2,
width = window_size,FUN = mean,align = "right",fill = NA)

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
plot(index(rollvar_sp500), rollvar_sp500, type = "l",
main = "Rolling Variance: SP500 Residuals",
xlab = "Date", ylab = "Variance (MA)")
abline(h = mean(na.omit(as.numeric(rollvar_sp500))), lty = 2)

plot(index(rollvar_nasdaq), rollvar_nasdaq, type = "l",
main = "Rolling Variance: NASDAQ Residuals",
xlab = "Date", ylab = "Variance (MA)")
abline(h = mean(na.omit(as.numeric(rollvar_nasdaq))), lty = 2)

# SP500 rolling variance stationarity check
rv_sp500 <- na.omit(as.numeric(rollvar_sp500))
fit_sp500_var <- lm(rv_sp500 ~ seq_along(rv_sp500))
p_sp500_var <- summary(fit_sp500_var)$coefficients[2, 4]
stationary_sp500 <- p_sp500_var >= 0.05
cat("SP500 rolling-variance trend test p-value =", signif(p_sp500_var, 4),
ifelse(stationary_sp500,"variance looks stationary.\n","variance looks non-stationary.\n"))

# NASDAQ rolling variance stationarity check
rv_nasdaq <- na.omit(as.numeric(rollvar_nasdaq))
fit_nasdaq_var <- lm(rv_nasdaq ~ seq_along(rv_nasdaq))
p_nasdaq_var <- summary(fit_nasdaq_var)$coefficients[2, 4]
stationary_nasdaq <- p_nasdaq_var >= 0.05
cat("NASDAQ rolling-variance trend test p-value =", signif(p_nasdaq_var, 4),
ifelse(stationary_nasdaq,"variance looks stationary.\n","variance looks non-stationary.\n"))

# Differentiation of the residuals
d_sp500 <- diff(resid_xts$SP500_resid)
d_nasdaq <- diff(resid_xts$NASDAQ_resid)

d_rollvar_sp500 <- zoo::rollapply(d_sp500^2,width = window_size,
FUN = mean,align = "right",fill = NA)

d_rollvar_nasdaq <- zoo::rollapply(d_nasdaq^2,width = window_size,
FUN = mean,align = "right",fill = NA)

plot(index(d_rollvar_sp500), d_rollvar_sp500, type = "l",
main = "Rolling Variance: Differenced SP500 Residuals",
xlab = "Date", ylab = "Variance (MA)")
abline(h = mean(na.omit(as.numeric(d_rollvar_sp500))), lty = 2)

plot(index(d_rollvar_nasdaq), d_rollvar_nasdaq, type = "l",
main = "Rolling Variance: Differenced NASDAQ Residuals",
xlab = "Date", ylab = "Variance (MA)")
abline(h = mean(na.omit(as.numeric(d_rollvar_nasdaq))), lty = 2)

# Differenced SP500 rolling variance stationarity check
rv_d_sp500 <- na.omit(as.numeric(d_rollvar_sp500))
fit_d_sp500_var <- lm(rv_d_sp500 ~ seq_along(rv_d_sp500))
p_d_sp500_var <- summary(fit_d_sp500_var)$coefficients[2, 4]
stationary_d_sp500 <- p_d_sp500_var >= 0.05
cat("Differenced SP500 rolling-variance trend test p-value =", signif(p_d_sp500_var, 4),
ifelse(stationary_d_sp500,"variance looks stationary.\n","variance looks non-stationary.\n"))

# Differenced NASDAQ rolling variance stationarity check
rv_d_nasdaq <- na.omit(as.numeric(d_rollvar_nasdaq))
fit_d_nasdaq_var <- lm(rv_d_nasdaq ~ seq_along(rv_d_nasdaq))
p_d_nasdaq_var <- summary(fit_d_nasdaq_var)$coefficients[2, 4]
stationary_d_nasdaq <- p_d_nasdaq_var >= 0.05
cat("Differenced NASDAQ rolling-variance trend test p-value =", signif(p_d_nasdaq_var, 4),
ifelse(stationary_d_nasdaq,"variance looks stationary.\n","variance looks non-stationary.\n"))

#5. Compute the autocorrelation functions for the residuals of two of your variables. Discuss the results. What type of noise process does the pattern suggest? AR type or MA type or a mixture? Does the autcorrelation have a cycle indicating a cyclical trend?
acf_lag_max <- 60

acf_sp500 <- acf(as.numeric(resid_xts$SP500_resid), lag.max = acf_lag_max, plot = FALSE)
acf_nasdaq <- acf(as.numeric(resid_xts$NASDAQ_resid), lag.max = acf_lag_max, plot = FALSE)

par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
plot(acf_sp500, main = "ACF: SP500 Detrended Residuals")
plot(acf_nasdaq, main = "ACF: NASDAQ Detrended Residuals")

# 95% significance threshold for ACF bars
n_sp500 <- length(na.omit(as.numeric(resid_xts$SP500_resid)))
n_nasdaq <- length(na.omit(as.numeric(resid_xts$NASDAQ_resid)))
bound_sp500 <- 1.96 / sqrt(n_sp500)
bound_nasdaq <- 1.96 / sqrt(n_nasdaq)

acf_vals_sp500 <- as.numeric(acf_sp500$acf)[-1]
acf_vals_nasdaq <- as.numeric(acf_nasdaq$acf)[-1]

sig_sp500 <- abs(acf_vals_sp500) > bound_sp500
sig_nasdaq <- abs(acf_vals_nasdaq) > bound_nasdaq

sig_count_sp500 <- sum(sig_sp500)
sig_count_nasdaq <- sum(sig_nasdaq)

lag1_sp500_sig <- sig_sp500[1]
lag1_nasdaq_sig <- sig_nasdaq[1]

tail_sig_sp500 <- sum(sig_sp500[2:10])
tail_sig_nasdaq <- sum(sig_nasdaq[2:10])

# SP500 interpretation
if (sig_count_sp500 <= 2) {
	sp500_noise_type <- "approximately white noise"
} else if (lag1_sp500_sig && tail_sig_sp500 == 0) {
	sp500_noise_type <- "MA-like (short memory, cutoff after lag 1)"
} else if (lag1_sp500_sig && tail_sig_sp500 >= 3) {
	sp500_noise_type <- "AR-like (persistent decay in ACF)"
} else {
	sp500_noise_type <- "mixed ARMA-type"
}

sp500_sign_change_count <- sum(diff(sign(acf_vals_sp500[1:20])) != 0, na.rm = TRUE)
sp500_cycle <- sp500_sign_change_count >= 6

cat("\nQ5 Interpretation SP500 Residuals:\n")
cat("Significant lags (|ACF| >", signif(bound_sp500, 3), "):", sig_count_sp500, "out of", acf_lag_max, "lags.\n")
cat("Suggested noise pattern:", sp500_noise_type, "\n")
cat("Cyclical indication:", ifelse(sp500_cycle,"possible oscillatory/cyclical component.","no strong cyclical pattern in ACF."), "\n")

# NASDAQ interpretation
if (sig_count_nasdaq <= 2) {
	nasdaq_noise_type <- "approximately white noise"
} else if (lag1_nasdaq_sig && tail_sig_nasdaq == 0) {
	nasdaq_noise_type <- "MA-like (short memory, cutoff after lag 1)"
} else if (lag1_nasdaq_sig && tail_sig_nasdaq >= 3) {
	nasdaq_noise_type <- "AR-like (persistent decay in ACF)"
} else {
	nasdaq_noise_type <- "mixed ARMA-type"
}

nasdaq_sign_change_count <- sum(diff(sign(acf_vals_nasdaq[1:20])) != 0, na.rm = TRUE)
nasdaq_cycle <- nasdaq_sign_change_count >= 6

cat("\nQ5 Interpretation NASDAQ Residuals:\n")
cat("Significant lags (|ACF| >", signif(bound_nasdaq, 3), "):", sig_count_nasdaq, "out of", acf_lag_max, "lags.\n")
cat("Suggested noise pattern:", nasdaq_noise_type, "\n")
cat("Cyclical indication:", ifelse(nasdaq_cycle,"possible oscillatory/cyclical component.","no strong cyclical pattern in ACF."), "\n")

#6. Plot the crosscorrelation function for the residuals for a pair of your variables. Interpret this result. Which variable seems to be anticpating the other?
ccf_lag_max <- 60

ccf_sp500_nasdaq <- ccf(
	as.numeric(resid_xts$SP500_resid),
	as.numeric(resid_xts$NASDAQ_resid),
	lag.max = ccf_lag_max,
	plot = FALSE
)

par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
plot(ccf_sp500_nasdaq, main = "CCF: SP500 vs NASDAQ Residuals")

ccf_vals <- as.numeric(ccf_sp500_nasdaq$acf)
ccf_lags <- as.numeric(ccf_sp500_nasdaq$lag)
n_ccf <- length(na.omit(as.numeric(resid_xts$SP500_resid)))
ccf_bound <- 1.96 / sqrt(n_ccf)

sig_ccf <- abs(ccf_vals) > ccf_bound
sig_count_ccf <- sum(sig_ccf)

idx_star <- which.max(abs(ccf_vals))
lag_star <- ccf_lags[idx_star]
ccf_star <- ccf_vals[idx_star]

direction_text <- c(
	"SP500 appears to anticipate NASDAQ",
	"strongest relationship is contemporaneous (no clear leader)",
	"NASDAQ appears to anticipate SP500"
)[sign(lag_star) + 2]

cat("\nQ6 Interpretation (SP500 vs NASDAQ residuals):\n")
cat("Significant CCF lags (|CCF| >", signif(ccf_bound, 3), "):", sig_count_ccf, "\n")
cat("Strongest cross-correlation at lag", lag_star,
		"with value", signif(ccf_star, 4), "\n")
cat("Manual rule: positive lag = NASDAQ leads; negative lag = SP500 leads; zero lag = contemporaneous.\n")
cat("Direction:", direction_text, "by", abs(lag_star), "period(s).\n")