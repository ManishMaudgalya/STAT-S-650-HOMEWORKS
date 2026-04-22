# ARMA in state-space form
# Model
# y_t = beta0 + beta1 * x_{t-1} + alpha1 * y_{t-1} + alpha2 * y_{t-2} + w_t + theta1 * w_{t-1} + theta2 * w_{t-2}

rm(list = ls())

#1. Parameters
beta0  <- 0.5
beta1  <- 1.2
alpha1 <- 0.95
alpha2 <- 0.35
theta1 <- 0.4
theta2 <- -0.2



#2. State-space representation
# State vector:
# s_t = [ y_t, y_{t-1}, w_t, w_{t-1} ]'
#
# s_t   = A s_{t-1} + B u_t + K * e_t   (here we just use w_t as noise)
# y_t   = C s_t + D u_t + w_t
#
# We encode the ARMA structure in A, B, C, D.

A <- matrix(c(
  alpha1, alpha2, theta1, theta2,
  1,       0,      0,      0,
  0,       0,      0,      0,
  0,       0,      1,      0
), nrow = 4, byrow = TRUE)

B <- matrix(c(beta1, 0, 0, 0), nrow = 4, byrow = TRUE)

C <- matrix(c(1, 0, 1, 0), nrow = 1)  # y_t depends on y_t and w_t in state
D <- beta0                            # intercept term

cat("A matrix:\n"); print(A)
cat("B matrix:\n"); print(B)
cat("C matrix:\n"); print(C)
cat("D (intercept):\n"); print(D)

set.seed(123)
n  <- 200
x  <- rnorm(n)
w  <- rnorm(n, sd = 1)

#3. Simulation from original ARMA
y <- numeric(n)

# Initialize with zeros
y[1:2] <- 0
w_lag1 <- 0
w_lag2 <- 0

for (t in 3:n) {
  y[t] <- beta0 +
    beta1 * x[t - 1] +
    alpha1 * y[t - 1] +
    alpha2 * y[t - 2] +
    w[t] +
    theta1 * w[t - 1] +
    theta2 * w[t - 2]
}

#4. Simulation in state-space form
s <- matrix(0, nrow = 4, ncol = n)
y_ss <- numeric(n)

for (t in 2:n) {
  noise_vec <- c(0, 0, w[t], 0)
  s[, t] <- A %*% s[, t - 1] + B * x[t - 1] + noise_vec
  y_ss[t] <- as.numeric(C %*% s[, t] + D + w[t])
}

#5. Stability for AR part
# AR polynomial: z^2 - alpha1 * z - alpha2 = 0
roots <- polyroot(c(1, -alpha1, -alpha2))
cat("\nAR(2) roots:\n"); print(roots)
cat("Moduli of roots:\n"); print(Mod(roots))

if (all(Mod(roots) > 1)) {
  cat("\nThe AR(2) part is STABLE\n")
} else {
  cat("\nThe AR(2) part is UNSTABLE\n")
}

#Plot
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))

plot(y, type = "l", col = "blue",
     main = "Simulated yt from ARMA",
     xlab = "t", ylab = "y_t")

plot(y_ss, type = "l", col = "red",
     main = "Simulated yt from state space form",
     xlab = "t", ylab = "y_t")

par(mfrow = c(1, 1))
