# Linear Difference Equation example
N=100
X1 = array(1:N, dim=c(100,1))
f = 2
X2 = sin(2*pi*(X1/N)*f)
# plot(X2,type="l",xlab="time")
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
X4 = array(std*rnorm(N), dim=c(N,1))
X5 = X2+X4
# plot(X5,type="1",xlab="time")

X = SX3;
Y=array( 0, dim=c(N,1))
y0 = 1 # -5   # 1
y_1 = 0
x0 = 0
# a1 = .5
 a1 = .9
# a2 = -.8
a2 = -.5

b = 1
s = .1
for (index in 1:N) {
  if (index == 1) {yt = a1*y0 +a2*y_1 + b*x0}
  if (index == 2) {yt = a1*Y[index-1] + a2*y0 + b*X[index-1] + 0*X4[index] }
  if (index >2 ){yt = a1*Y[index-1] + a2*Y[index-2] + b*X[index-1] + 0*X4[index] }
  Y[index] = yt
 
}
 B = c(0,b)
 X = filter(X,B) + 0*X4       # wgted lag inputs using B wgts
 a <- c(a1,a2)          # vector
 YI = c(y0,y_1)
Z = filter(X,a,method=c("recursive"),sides=1,circular=FALSE,YI)  # feedback model 

# windows(7,10)
# layout( matrix( c(1,2), nrow=2, ncol=1, byrow="FALSE"))
# plot(Y,type="l",xlab="time")
# plot(Z,type="l",xlab="time") 
# plot(SX3,type="l",xlab="time")
M = cbind(SX3,Y,Z)
plot(M,type="b",xlab="time",main="diff eq")
