set.seed(42)
T <- 200
#Generate x(t)
x <- ifelse(1:T <= 100, 1:T, 0)
#Generate y(t)
w <- rnorm(T, mean = 0, sd = 1)
 
y <- numeric(T)
y[1] <- 2
y[2] <- 4

for (t in 3:T) {
  y[t] <- 10 + 0.75 * y[t-1] - 0.65 * y[t-2] + 0.15 * x[t-1] + 5 * w[t]
}

# Fitting the model for t > 2:
#y(t) = b0 + b1*y(t-1) + b2*y(t-2) + b3*x(t-1) + error
# For values t = 3 to 200 
t_idx <- 3:T
df <- data.frame(y_t   = y[t_idx],y_lag1 = y[t_idx - 1],y_lag2 = y[t_idx - 2],
                 x_lag1 = x[t_idx - 1]
)

fit <- lm(y_t ~ y_lag1 + y_lag2 + x_lag1, data = df)
summary(fit)


# a) Comparing estimated vs true parameters

cat("Comparing Estimated and True Values \n")
true_params <- c(b0 = 10, b1 = 0.75, b2 = -0.65, b3 = 0.15, sqrt_MSE = 5)
est         <- coef(fit)
mse         <- mean(fit$residuals^2)
sqrt_mse    <- sqrt(mse)

cat("b0","\nTrue = ",true_params["b0"],"\nEstimated = ",est["(Intercept)"])
cat("b1","\nTrue = ",true_params["b1"],"\nEstimated = ",est["y_lag1"])
cat("b2","\nTrue = ",true_params["b2"],"\nEstimated = ",est["y_lag2"])
cat("b3","\nTrue = ",true_params["b3"],"\nEstimated = ",est["x_lag1"])
cat("sqrt(MSE)","\nTrue = ", true_params["sqrt_MSE"], "\nEstimated = ",sqrt_mse)


# R-squared
r2 <- summary(fit)$r.squared
cat("R-squared = \n", r2)


# Plot predicted vs observed

y_hat <- fitted(fit)
y_obs <- df$y_t

par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

# Plot 1: Time series overlay
plot(t_idx, y_obs, type = "l", col = "blue", lwd = 1.5,
     xlab = "t", ylab = "y(t)",
     main = "Observed vs Predicted over Time")
lines(t_idx, y_hat, col = "red", lwd = 1.5, lty = 2)
legend("topright", legend = c("Observed", "Predicted"),
       col = c("blue", "red"), lty = c(1, 2), lwd = 1.5, cex = 0.8)

# Plot 2: Scatter predicted vs observed
plot(y_obs, y_hat, pch = 16, cex = 0.5, col = "blue",
     xlab = "Observed y(t)", ylab = "Predicted y(t)",
     main = "Predicted vs Observed Scatter plot")
abline(0, 1, col = "red", lwd = 2)
legend("topleft", legend = sprintf("R² = %.4f", r2), bty = "n", cex = 0.9)
