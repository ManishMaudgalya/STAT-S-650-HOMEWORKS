# Spectral Analysis: S&P 500 (GSPC) vs NASDAQ (IXIC)
# Steps 3, 4, 5: Periodogram, Spectral Density, Coherence

library(ggplot2)
library(dplyr)
library(tidyr)
library(astsa)
df <- read.csv("/Users/manishsm/Desktop/Spring 2026/Time Series Analysis/Assignments/quantmod_dataset.csv")
df$Date <- as.Date(df$Date)
df <- df[order(df$Date), ]

# Use log returns to achieve stationarity
gspc <- diff(log(df$GSPC.Adjusted))
ixic <- diff(log(df$IXIC.Adjusted))

n <- length(gspc)
cat("Number of observations of log returns:", n, "\n")

# STEP 3: Periodogram
# Compute raw periodograms
pg_gspc <- spec.pgram(gspc, taper = 0, log = "no", plot = FALSE)
pg_ixic <- spec.pgram(ixic, taper = 0, log = "no", plot = FALSE)

# Convert frequency to cycles per year (252 trading days/year)
freq_year_gspc <- pg_gspc$freq * 252
freq_year_ixic <- pg_ixic$freq * 252

# Plot GSPC Periodogram
par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))

plot(freq_year_gspc, pg_gspc$spec,
     type = "l", col = "steelblue",
     xlab = "Frequency",
     ylab = "Spectral Density",
     main = "Raw Periodogram S&P500 Log Returns")
abline(v = c(1, 2, 4, 12, 52, 252), col = "red", lty = 2, lwd = 0.8)
legend("topright",
       legend = c("Periodogram", "Key frequencies"),
       col    = c("steelblue", "red"),
       lty    = c(1, 2), bty = "n")

#Plot IXIC Periodogram
plot(freq_year_ixic, pg_ixic$spec,
     type = "l", col = "darkorange",
     xlab = "Frequency",
     ylab = "Spectral Density",
     main = "Raw Periodogram NASDAQLog Returns")
abline(v = c(1, 2, 4, 12, 52, 252), col = "red", lty = 2, lwd = 0.8)
legend("topright",
       legend = c("Periodogram", "Key frequencies"),
       col    = c("darkorange", "red"),
       lty    = c(1, 2), bty = "n")

#Check for cycles
# Find top 5 peaks in each periodogram
top_gspc <- order(pg_gspc$spec, decreasing = TRUE)[1:5]
top_ixic <- order(pg_ixic$spec, decreasing = TRUE)[1:5]

cat("\nTop 5 dominant frequencies S&P 500in cycles/year:\n")
print(data.frame(
  freq_per_year = round(pg_gspc$freq[top_gspc] * 252, 3),
  period_days   = round(1 / pg_gspc$freq[top_gspc], 1),
  power         = round(pg_gspc$spec[top_gspc], 6)
))

cat("\nTop 5 dominant frequencies NASDAQ in cycles/year:\n")
print(data.frame(
  freq_per_year = round(pg_ixic$freq[top_ixic] * 252, 3),
  period_days   = round(1 / pg_ixic$freq[top_ixic], 1),
  power         = round(pg_ixic$spec[top_ixic], 6)
))


remove_cycle <- function(x, freq_dominant) {
  # freq_dominant: frequency in cycles per observation (not per year)
  t <- seq_along(x)
  cos_term <- cos(2 * pi * freq_dominant * t)
  sin_term <- sin(2 * pi * freq_dominant * t)
  fit <- lm(x ~ cos_term + sin_term)
  cat("Regression summary for cycle removal:\n")
  print(summary(fit))
  residuals(fit)
}

# gspc_detrended <- remove_cycle(gspc, freq_dominant = pg_gspc$freq[top_gspc[1]])
# ixic_detrended <- remove_cycle(ixic, freq_dominant = pg_ixic$freq[top_ixic[1]])

# 4 Smoothed Spectral Density (different window sizes)
par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))

spans_list <- list(
  "Narrow (m=3)"  = c(3),
  "Medium (m=9)"  = c(9),
  "Wide (m=21)"   = c(21)
)

for (nm in names(spans_list)) {
  sp <- spans_list[[nm]]

  sm_gspc <- spec.pgram(gspc, spans = sp, taper = 0.1, log = "no", plot = FALSE)
  sm_ixic <- spec.pgram(ixic, spans = sp, taper = 0.1, log = "no", plot = FALSE)

  plot(sm_gspc$freq * 252, sm_gspc$spec,
       type = "l", col = "steelblue",
       xlab = "Frequency (cycles/year)", ylab = "Spectral Density",
       main = paste("GSPC –", nm))
  abline(v = c(1, 5, 52), col = "red", lty = 2, lwd = 0.8)

  plot(sm_ixic$freq * 252, sm_ixic$spec,
       type = "l", col = "darkorange",
       xlab = "Frequency (cycles/year)", ylab = "Spectral Density",
       main = paste("IXIC –", nm))
  abline(v = c(1, 5, 52), col = "red", lty = 2, lwd = 0.8)
}

#CHOSEN window: Medium
cat("\nFinal chosen spectral estimates (spans = 9, taper = 10%)\n")
final_gspc <- spec.pgram(gspc, spans = 9, taper = 0.1, log = "no", plot = FALSE)
final_ixic <- spec.pgram(ixic, spans = 9, taper = 0.1, log = "no", plot = FALSE)

par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))

plot(final_gspc$freq * 252, final_gspc$spec,
     type = "l", col = "steelblue", lwd = 1.5,
     xlab = "Frequency (cycles/year)", ylab = "Spectral Density",
     main = "Smoothed Spectral Density – GSPC (spans=9, taper=10%)")
abline(v = c(1, 5, 12, 52), col = "red", lty = 2)
legend("topright",
       legend = c("Spectral density", "Key annual freqs"),
       col = c("steelblue", "red"), lty = c(1, 2), bty = "n")

plot(final_ixic$freq * 252, final_ixic$spec,
     type = "l", col = "darkorange", lwd = 1.5,
     xlab = "Frequency (cycles/year)", ylab = "Spectral Density",
     main = "Smoothed Spectral Density – IXIC (spans=9, taper=10%)")
abline(v = c(1, 5, 12, 52), col = "red", lty = 2)
legend("topright",
       legend = c("Spectral density", "Key annual freqs"),
       col = c("darkorange", "red"), lty = c(1, 2), bty = "n")

# 5 Coherence Between GSPC and IXIC

par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))

coh <- spec.pgram(cbind(gspc, ixic),
                  spans   = 9,
                  taper   = 0.1,
                  log     = "no",
                  plot    = FALSE)

# Squared coherence
sq_coh <- Mod(coh$coh)^2

plot(coh$freq * 252, sq_coh,
     type = "l", col = "purple", lwd = 1.5,
     ylim = c(0, 1),
     xlab = "Frequency (cycles/year)",
     ylab = expression(hat(C)^2 ~ "(Squared Coherence)"),
     main = "Squared Coherence: GSPC vs IXIC (spans=9, taper=10%)")


m_df <- length(spans_list[["Medium (m=9)"]])
# Degrees of freedom for spans=9: df ≈ 2*L where L=spans
L <- 9
threshold_95 <- 1 - (0.05)^(1 / (L - 1))
abline(h = threshold_95, col = "red", lty = 2, lwd = 1.5)
legend("topright",
       legend = c("Sq. Coherence", paste0("95% threshold ", round(threshold_95, 3), ")")),
       col = c("purple", "red"), lty = c(1, 2), bty = "n")

cat("\n95% coherence significance threshold:", round(threshold_95, 4), "\n")
cat("Mean squared coherence across frequencies:", round(mean(sq_coh), 4), "\n")
cat("Frequencies (cycles/year) where coherence > threshold:\n")
sig_freqs <- coh$freq[sq_coh > threshold_95] * 252
cat(round(range(sig_freqs), 2), "(range)\n")
cat("Proportion of frequencies with significant coherence:",
    round(mean(sq_coh > threshold_95), 3), "\n")