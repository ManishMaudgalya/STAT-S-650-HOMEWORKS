
# Impulse response example
N=1000
# X1 = array(1:N, dim=c(N,1))
X1 = 1:N
f = 2
X2 = sin(2*pi*(X1/N)*f)
 plot(X2,type="l",xlab="time")
 X3 = array( 0, dim=c(N,1))
X3[N/2]=1
 plot(X3,type="l",xlab="time")

SX3=array( 0, dim=c(N,1))
sum = 0
for (index in 1:N) {
  sum = sum + X3[index]
  SX3[index] = sum 
}
 plot(SX3,type="l",xlab="time") 

std = 1
X4 = std*rnorm(N)
# X4 = array(std*rnorm(N), dim=c(N,1))
X5 = X2+ X4
 plot(X5,type="l",xlab="time")

width = 25
w1 = width+1

w = array( (1/width), dim=c(width,1))
Y=array( 0, dim=c(N-(width-1),1))
nrow(Y)

for (index in (w1):(N+1)) {
  yt=0
  for ( Ind in 1:width )  {
       yt = yt + w[Ind]*X5[index-Ind]
      }
  Y[index-width] = yt
 # print(Y[index-width])
}

 ws = 25
 Yf = filter(X5, rep(1/ws, ws), sides = 1)
 # set sides = 1 for impulse response model. 
 # Yf[n]= 1/ws*X5[n] + 1/ws*X5[n-1] +...+ 1/ws*X5[n-(ws-1)]. 
 #  length(Yf) = 1000 
 # first ws-1 are NA, Yf[ws] is first real number, ends at N 
 # set sides = 2 for two sided filter model. 
 # Yf[n]= 1/ws*X5[n+ (ws-1)/2] +...+1/ws*X5[n]+...1/ws*X5[n-(ws-1)/2]
 # length(Yf) = 1000
 # first (ws-1)/2 are NA, Yf[ (ws-1)/2 +1] is first real number, Yf[N-(ws-1)/2 + 1] is last real number, 
 # last (ws-1)/2 elements are NA
 
 
 # windows(7,10)
 par(mar=c(4,4,4,4))
 layout( matrix( c(2,1), nrow=2, ncol=1, byrow="FALSE"))
 plot(Y,type="l",xlab="time") 
 plot(Yf,type="l",xlab="time") 
# plot(X5,type="l",xlab="time")

# par(mfrow=c(3,1))
#  acf(Y, 10, main="output")
#   acf(X5, 10, main="input")
#  ccf(Yf[25:1000], X5[25:1000], 10, main="Output vs Input", ylab="CCF")


