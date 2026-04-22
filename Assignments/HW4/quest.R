file_arg <- "--file="
script_args <- commandArgs(trailingOnly = FALSE)
script_match <- grep(file_arg, script_args)

if (length(script_match) == 0) {
	script_dir <- getwd()
} else {
	script_path <- sub(file_arg, "", script_args[script_match[1]])
	script_dir <- dirname(normalizePath(script_path))
}

trading_days <- 252
output_dir <- script_dir
data_path <- file.path(dirname(script_dir), "quantmod_dataset.csv")

if (!file.exists(data_path)) {
	stop("Could not find quantmod_dataset.csv next to the HW folders.")
}

df <- read.csv(data_path, stringsAsFactors = FALSE)
df$Date <- as.Date(df$Date)
df <- df[complete.cases(df[, c("GSPC.Adjusted", "IXIC.Adjusted")]), ]

time_index <- seq_len(nrow(df))
series_names <- c("GSPC", "IXIC")
series_values <- list(
	GSPC = df$GSPC.Adjusted,
	IXIC = df$IXIC.Adjusted
)

detrended <- list()
for (name in series_names) {
	fit <- lm(series_values[[name]] ~ poly(time_index, 2, raw = TRUE))
	detrended[[name]] <- residuals(fit)
}

cycle_threshold_ratio <- 15
periodogram_results <- list()

for (name in series_names) {
	x <- detrended[[name]]
	raw_spec_obj <- spec.pgram(
		ts(x, frequency = trading_days),
		taper = 0.1,
		fast = TRUE,
		demean = TRUE,
		detrend = FALSE,
		plot = FALSE,
		log = "no"
	)
	raw_positive <- raw_spec_obj$freq > 0
	raw_freq <- raw_spec_obj$freq[raw_positive]
	raw_spec <- raw_spec_obj$spec[raw_positive]
	raw_peak_idx <- which.max(raw_spec)
	raw_peak_freq <- raw_freq[raw_peak_idx]
	raw_peak_period_days <- trading_days / raw_peak_freq
	raw_peak_ratio <- raw_spec[raw_peak_idx] / median(raw_spec)

	cycle_removed <- raw_peak_ratio >= cycle_threshold_ratio
	adjusted_series_current <- x

	if (cycle_removed) {
		harmonic_data <- data.frame(
			x = x,
			sin_term = sin(2 * pi * raw_peak_freq * time_index / trading_days),
			cos_term = cos(2 * pi * raw_peak_freq * time_index / trading_days)
		)
		harmonic_fit <- lm(x ~ sin_term + cos_term, data = harmonic_data)
		adjusted_series_current <- residuals(harmonic_fit)
	}

	adjusted_spec_obj <- spec.pgram(
		ts(adjusted_series_current, frequency = trading_days),
		taper = 0.1,
		fast = TRUE,
		demean = TRUE,
		detrend = FALSE,
		plot = FALSE,
		log = "no"
	)
	adjusted_positive <- adjusted_spec_obj$freq > 0
	adjusted_freq <- adjusted_spec_obj$freq[adjusted_positive]
	adjusted_spec <- adjusted_spec_obj$spec[adjusted_positive]
	adjusted_peak_idx <- which.max(adjusted_spec)
	adjusted_peak_freq <- adjusted_freq[adjusted_peak_idx]
	adjusted_peak_period_days <- trading_days / adjusted_peak_freq

	periodogram_results[[name]] <- list(
		raw_freq = raw_freq,
		raw_spec = raw_spec,
		raw_peak_freq = raw_peak_freq,
		raw_peak_period_days = raw_peak_period_days,
		raw_peak_ratio = raw_peak_ratio,
		cycle_removed = cycle_removed,
		adjusted_series = adjusted_series_current,
		adjusted_freq = adjusted_freq,
		adjusted_spec = adjusted_spec,
		adjusted_peak_freq = adjusted_peak_freq,
		adjusted_peak_period_days = adjusted_peak_period_days
	)
}

png(file.path(output_dir, "periodograms.png"), width = 1600, height = 1000, res = 160)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

for (name in series_names) {
	result <- periodogram_results[[name]]

	plot(
		result$raw_freq,
		result$raw_spec,
		type = "l",
		col = "steelblue",
		lwd = 1.4,
		xlab = "Frequency (cycles per year)",
		ylab = "Periodogram",
		main = paste0(name, " raw periodogram")
	)
	abline(v = result$raw_peak_freq, col = "firebrick", lty = 2)
	mtext(
		paste0("Peak period ~= ", round(result$raw_peak_period_days, 1), " trading days"),
		side = 3,
		line = 0.1,
		cex = 0.8
	)

	if (result$cycle_removed) {
		adjusted_title <- paste0(name, " after harmonic regression")
	} else {
		adjusted_title <- paste0(name, " unchanged (no dominant cycle)")
	}

	plot(
		result$adjusted_freq,
		result$adjusted_spec,
		type = "l",
		col = "darkgreen",
		lwd = 1.4,
		xlab = "Frequency (cycles per year)",
		ylab = "Periodogram",
		main = adjusted_title
	)
	abline(v = result$adjusted_peak_freq, col = "firebrick", lty = 2)
}

dev.off()

cat("Question 3: Periodograms\n")
for (name in series_names) {
	result <- periodogram_results[[name]]
	cat(
		sprintf(
			"%s: strongest peak at %.2f cycles/year (period %.2f trading days), peak/median ratio %.2f.\n",
			name,
			result$raw_peak_freq,
			result$raw_peak_period_days,
			result$raw_peak_ratio
		)
	)

	if (result$cycle_removed) {
		cat(
			sprintf(
				"%s: removed that dominant harmonic with sin/cos regression; new strongest period is %.2f trading days.\n",
				name,
				result$adjusted_peak_period_days
			)
		)
	} else {
		cat(sprintf("%s: no single cycle stood out enough to justify regression removal.\n", name))
	}
}

adjusted_series <- list(
	GSPC = periodogram_results$GSPC$adjusted_series,
	IXIC = periodogram_results$IXIC$adjusted_series
)

spans_candidates <- list(c(5, 5), c(9, 9), c(13, 13))
candidate_labels <- c("(5,5)", "(9,9)", "(13,13)")
chosen_spans <- c(9, 9)
chosen_label <- "(9,9)"
colors <- c("gray55", "black", "darkorange3")

spectral_distribution_results <- list()
for (name in series_names) {
	spectral_distribution_results[[name]] <- list()

	for (idx in seq_along(spans_candidates)) {
		spans_value <- spans_candidates[[idx]]
		label <- candidate_labels[idx]
		smoothed_spec_obj <- spec.pgram(
			ts(adjusted_series[[name]], frequency = trading_days),
			spans = spans_value,
			taper = 0.1,
			fast = TRUE,
			demean = TRUE,
			detrend = FALSE,
			plot = FALSE,
			log = "no"
		)
		smoothed_positive <- smoothed_spec_obj$freq > 0
		smoothed_freq <- smoothed_spec_obj$freq[smoothed_positive]
		smoothed_spec <- smoothed_spec_obj$spec[smoothed_positive]
		cum_dist <- cumsum(smoothed_spec) / sum(smoothed_spec)
		peak_idx <- which.max(smoothed_spec)
		low_frequency_share <- sum(smoothed_spec[smoothed_freq <= 4]) / sum(smoothed_spec)

		spectral_distribution_results[[name]][[label]] <- list(
			freq = smoothed_freq,
			cum_dist = cum_dist,
			peak_period_days = trading_days / smoothed_freq[peak_idx],
			low_frequency_share = low_frequency_share
		)
	}
}

png(file.path(output_dir, "spectral_distributions.png"), width = 1400, height = 700, res = 160)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

for (name in series_names) {
	chosen_result <- spectral_distribution_results[[name]][[chosen_label]]

	plot(
		chosen_result$freq,
		chosen_result$cum_dist,
		type = "n",
		ylim = c(0, 1),
		xlab = "Frequency (cycles per year)",
		ylab = "Estimated spectral distribution",
		main = paste0(name, " cumulative spectrum")
	)

	for (idx in seq_along(candidate_labels)) {
		label <- candidate_labels[idx]
		line_width <- 1.2
		line_type <- 2
		if (label == chosen_label) {
			line_width <- 2.2
			line_type <- 1
		}

		lines(
			spectral_distribution_results[[name]][[label]]$freq,
			spectral_distribution_results[[name]][[label]]$cum_dist,
			col = colors[idx],
			lwd = line_width,
			lty = line_type
		)
	}

	legend(
		"bottomright",
		legend = paste("spans", candidate_labels),
		col = colors,
		lwd = c(1.2, 2.2, 1.2),
		lty = c(2, 1, 2),
		bty = "n",
		cex = 0.9
	)
	}

dev.off()

cat("\nQuestion 4: Estimated spectral distributions\n")
cat(
	sprintf(
		"Tried Daniell windows with spans %s, %s, and %s; used spans %s for interpretation because it smooths the jagged raw periodogram without flattening the dominant medium-run peak.\n",
		candidate_labels[1],
		candidate_labels[2],
		candidate_labels[3],
		chosen_label
	)
)

for (name in series_names) {
	chosen_result <- spectral_distribution_results[[name]][[chosen_label]]
	cat(
		sprintf(
			"%s: chosen smoothed peak at %.2f trading days; low-frequency share (<= 4 cycles/year) = %.3f. Most power remains in low to medium frequencies, so the dominant behavior is broad market swings with one recurring medium-run cycle rather than a sharp seasonal pattern.\n",
			name,
			chosen_result$peak_period_days,
			chosen_result$low_frequency_share
		)
	)
}

coherence_input <- ts(
	cbind(GSPC = adjusted_series$GSPC, IXIC = adjusted_series$IXIC),
	frequency = trading_days
)

coherence_spec <- spec.pgram(
	coherence_input,
	spans = chosen_spans,
	taper = 0.1,
	fast = TRUE,
	demean = TRUE,
	detrend = FALSE,
	plot = FALSE,
	log = "no"
)

positive <- coherence_spec$freq > 0
coh_freq <- coherence_spec$freq[positive]
coh_values <- coherence_spec$coh[positive, 1]
phase_values <- coherence_spec$phase[positive, 1]
top_idx <- order(coh_values, decreasing = TRUE)[1:5]
dominant_shared_freq <- mean(c(periodogram_results$GSPC$raw_peak_freq, periodogram_results$IXIC$raw_peak_freq))
dominant_coh_idx <- which.min(abs(coh_freq - dominant_shared_freq))
low_frequency_band <- coh_freq <= 4

png(file.path(output_dir, "coherence.png"), width = 1200, height = 700, res = 160)
par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
plot(
	coh_freq,
	coh_values,
	type = "l",
	col = "purple4",
	lwd = 1.6,
	ylim = c(0, 1),
	xlab = "Frequency (cycles per year)",
	ylab = "Coherence",
	main = paste0("Coherence: GSPC vs IXIC (spans ", chosen_label, ")")
)
abline(h = c(0.5, 0.8), col = c("gray70", "gray50"), lty = 3)

dev.off()

cat("\nQuestion 5: Coherence\n")
cat(sprintf("Mean coherence across frequencies: %.3f.\n", mean(coh_values)))
cat(sprintf("Mean coherence over low frequencies (<= 4 cycles/year): %.3f.\n", mean(coh_values[low_frequency_band])))
cat(
	sprintf(
		"Coherence at the shared dominant cycle (%.2f cycles/year, period %.2f trading days): %.3f.\n",
		coh_freq[dominant_coh_idx],
		trading_days / coh_freq[dominant_coh_idx],
		coh_values[dominant_coh_idx]
	)
)

for (idx in top_idx) {
	cat(
		sprintf(
			"High coherence at %.2f cycles/year (period %.2f trading days): coherence %.3f, phase %.3f radians.\n",
			coh_freq[idx],
			trading_days / coh_freq[idx],
			coh_values[idx],
			phase_values[idx]
		)
	)
}

cat("Interpretation: the two indices move together very strongly in the frequency domain. The low-frequency coherence is high, which means broad market swings are largely shared across the two indices, and the coherence at the common dominant cycle is also strong. The highest coherence spikes at very short periods show that some fast fluctuations are synchronized too, but those narrow peaks are less economically important than the broad low-frequency comovement.\n")

cat("\nSaved plots to:\n")
cat(file.path(output_dir, "periodograms.png"), "\n")
cat(file.path(output_dir, "spectral_distributions.png"), "\n")
cat(file.path(output_dir, "coherence.png"), "\n")
