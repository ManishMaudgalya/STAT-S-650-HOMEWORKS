# Task 9: State Space Model (Problem 3)
# Using S&P 500 Data
library(quantmod)

# Pull S&P 500 data
getSymbols("^GSPC", src = "yahoo", from = "1990-01-01", to = Sys.Date(), 
           auto.assign = TRUE, warnings = FALSE)

# Extract closing price and compute returns
X_price = as.numeric(Cl(GSPC))
X = diff(log(X_price)) * 100 
X = c(0, X)
N = length(X)

# Problem 3 coefficients (from difference equation)
a1 = 0.9
a2 = -0.5

# State space model
# F = State transition matrix
F = array(c(a1, a2, 1, 0), dim=c(2,2))
# G = Input matrix
G = array(c(1, 0), dim=c(2,1))

cat("State Transition Matrix F:\n")
print(F)
cat("Input Matrix G:\n")
print(G)

# eigenvalues
ev = eigen(F)
cat("\nEigenvalues of F:\n")
print(ev$values)
cat("Absolute values:\n")
print(abs(ev$values))

# Initialize state
Z = array(c(0, 0), dim=c(2,1))
ZM = array(0, dim=c(2,N))

ZM[,1] = F%*%Z + G*X[1]

for (time in 2:N) {
  LZ = ZM[,time-1]
  LZ = array(LZ, dim=c(2,1))
  Z = F%*%LZ + G*X[time]
  ZM[,time] = Z
}

ZM = t(ZM)
Ys = ZM[,1]

# Compare with DE
y0 = 0
y_1 = 0
Y = array(0, dim=c(N,1))

for (index in 1:N) {
  if (index == 1) {yt = a1*y0 + a2*y_1 + 0}
  if (index == 2) {yt = a1*Y[index-1] + a2*y0 + X[index-1]}
  if (index > 2) {yt = a1*Y[index-1] + a2*Y[index-2] + X[index-1]}
  Y[index] = yt
}

par(mfrow=c(4,1), mar=c(2,2,2,2))
plot(X, type="l", xlab="time", main="S&P 500 Log Returns (%)")
plot(Ys, type="l", xlab="time", main="State Space Model Output")
plot(Y, type="l", xlab="time", main="Difference Equation Output")
plot(cbind(Ys, Y), type="l", xlab="time", main="State Space vs Difference Equation", col=c("blue", "red"))
legend("topright", c("State Space", "Diff Eq"), col=c("blue", "red"), lty=c(1,1))

# Comparison
diff = Ys - as.numeric(Y)
cat("\nComparison Results:\n")
cat("Max absolute difference:", max(abs(diff), na.rm=TRUE), "\n")
cat("Mean squared error:", mean(diff^2, na.rm=TRUE), "\n")

if (max(abs(diff), na.rm=TRUE) < 1e-10) {
  cat("\nBoth methods produce equivalent results!\n")
}