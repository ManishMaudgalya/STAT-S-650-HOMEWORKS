# State Space model 
N=100
X1 = array(1:N, dim=c(N,1))
f = 2
X2 = sin(2*pi*(X1/N)*f)
# plot(X2,type="l",xlab="time")
X3 = array( 0, dim=c(N,1))
X3[N/2]=1
# plot(X3,type="l",xlab="time")

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

X = SX3;
plot(X,type="l",xlab="time")


X[1,1] = 1;



# F = c(.75,.25,-.5,.7)
 F  = c(.9,-.5,1,0)   # .9 -.5
F = array( F, dim = c(2,2))
F 
ev = eigen(F)
# Z = c(.5,.5)
ev$values
abs(ev$values)

Z = c(0,0);
Z = array( Z, dim = c(2,1))
ZM = array( 0, dim = c(2,N))
# G = c(0,1)
 G = c(1,0)
G = array( G, dim=c(2,1))

ZM[,1] = F%*%Z + G*X[1]

for (time in 2:N) {
  LZ = ZM[,time-1]
  LZ = array( LZ, dim = c(2,1))
  Z = F%*%LZ + G%*%X[time] 
  ZM[,time]=Z
 }
ZM = t(ZM)

Ys = ZM[,1]

par(mar=c(2,2,2,2))
 # windows(7,10)
 layout( matrix( c(1,2), nrow=2, ncol=1, byrow="FALSE"))
 plot(Ys,type="l",xlab="time") 
 plot(SX3,type="l",xlab="time")

 plot(ZM,type="l")

# windows(7,10)
# layout( matrix( c(1,2), nrow=2, ncol=1, byrow="FALSE"))
# plot(ZM[,1],type="l",xlab="time") 
# plot(ZM[,2],type="l",xlab="time")
 
 y0 = 1
 y_1 = 0
 YI = c(y0,y_1)
 b = c(.9,-.5)
 Yd = filter(SX3,b,method=c("recursive"),sides=1,circular=FALSE,YI)  # feedback model 
 plot(Yd,type='l',xlab='time')

 M = cbind(Ys,Yd)
 plot(M,type="l",xlab="time",main="state space")
 
 