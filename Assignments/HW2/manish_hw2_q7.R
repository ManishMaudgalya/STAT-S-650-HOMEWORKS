
# Using S&P 500 Data and smoothening using different window sizes

library(quantmod)
getSymbols("^GSPC", src = "yahoo", from = "1990-01-01", to = Sys.Date(), 
           auto.assign = TRUE, warnings = FALSE)

# Closing price
X5 = as.numeric(Cl(GSPC))
N = length(X5)

# First window size: ws1 = 10
ws1 = 10
w1 = ws1 + 1
Yf1 = filter(X5, rep(1/ws1, ws1), sides = 1)

# Second window size: ws2 = 30
ws2 = 30
w2 = ws2 + 1
Yf2 = filter(X5, rep(1/ws2, ws2), sides = 1)

# Third window size: ws3 = 60
ws3 = 60
w3 = ws3 + 1
Yf3 = filter(X5, rep(1/ws3, ws3), sides = 1)


par(mfrow=c(4,1), mar=c(2,2,2,2))
plot(X5, type="l", xlab="time", main="S&P 500 Original Data")
plot(Yf1, type="l", xlab="time", main=paste("Smoothed (window size =", ws1, ")"))
plot(Yf2, type="l", xlab="time", main=paste("Smoothed (window size =", ws2, ")"))
plot(Yf3, type="l", xlab="time", main=paste("Smoothed (window size =", ws3, ")"))