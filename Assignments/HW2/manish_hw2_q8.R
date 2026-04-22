# Difference Equation: Y_t = 0.9*Y_{t-1} - 0.5*Y_{t-2} + X_{t-1}

library(quantmod)

# Pull S&P 500 data
getSymbols("^GSPC", src = "yahoo", from = "1990-01-01", to = Sys.Date(), 
           auto.assign = TRUE, warnings = FALSE)

# Extract closing price and compute returns (to avoid scale issues)
X_price = as.numeric(Cl(GSPC))
X = diff(log(X_price)) * 100  # Log returns in percentage
X = c(0, X)  # Pad to match length
N = length(X)

# Problem 3 coefficients
a1 = 0.9
a2 = -0.5
b = 1

# Check eigenvalues
F = array(c(a1, a2, 1, 0), dim=c(2,2))
ev = eigen(F)
cat("Eigenvalues:", ev$values, "\n")
cat("Absolute values:", abs(ev$values), "\n")

# Simulate difference equation
Y = array(0, dim=c(N,1))
y0 = 0
y_1 = 0
x0 = 0

for (index in 1:N) {
  if (index == 1) {yt = a1*y0 + a2*y_1 + b*x0}
  if (index == 2) {yt = a1*Y[index-1] + a2*y0 + b*X[index-1]}
  if (index > 2) {yt = a1*Y[index-1] + a2*Y[index-2] + b*X[index-1]}
  Y[index] = yt
}

# Also compute using filter function
a = c(a1, a2)
YI = c(y0, y_1)
Z = filter(X, a, method=c("recursive"), sides=1, circular=FALSE, YI)

# Plot results
par(mfrow=c(3,1), mar=c(2,2,2,2))
plot(X, type="l", xlab="time", main="S&P 500 Log Returns (%)")
plot(Y, type="l", xlab="time", main="Difference Equation Output (Direct Simulation)")
plot(cbind(Y, Z), type="l", xlab="time", main="Direct Simulation vs filter()")
legend("topright", c("Direct", "filter()"), col=c("black", "red"), lty=c(1,1))

cat("\nSystem Stability Analysis:\n")
if (all(abs(ev$values) < 1)) {
  cat("System is STABLE (all eigenvalues < 1)\n")
  cat("Output should converge to steady state\n")
} else {
  cat("System is UNSTABLE (eigenvalue >= 1)\n")
}