#1a

#Theoretical value of B
Nsim=10^5
set.seed(123)
x <- runif(Nsim)
f=function(x){
  1/(1+x^2)
}
integrate(f, 0, Inf) #1.570796
act_int <- 1.570796
realb <- 1/act_int #0.6366199

#Simulation method of B
set.seed(123)
Nsim=10^5
u=runif(Nsim)
x=2*sqrt(1-u^2)
est_int=mean(x) #1.5719759
est_int
B <- 1/est_int #0.63614

#Measure of approximation accuracy(of our estimation): Standard Error + 95% CI
#diff between estimated and actual value of 1/B (x)
diff=(1/B)-(1/realb)
#SE of 1/B (x)
se=sqrt(mean(x^2)-((est_int)^2))/sqrt(Nsim) #sd=sqrt(Var[X])=sqrt(E[X^2]-{E[X]}^2})
#CI
LCI = (est_int)-1.96*se         #Lower Bound of CI
ULI = (est_int)+1.96*se         #Upper Bound of CI
compare_oneoverb <- rbind(1/realb, 1/B, diff, se, LCI, ULI)
row.names(compare_oneoverb) <- c("Actual value (1/B)", "Simulated Value (1/B)", 
                                 "Difference (1/B)", "Std Error", "Lower 95% CI", "Upper 95% CI" )
print(compare_oneoverb)

#diff between estimated and actual value of B
diff=B-realb
diff #-0.0004778628

###############################################################################################
###############################################################################################
#1b

#ITM
#estimated
set.seed(123)
u=runif(Nsim)
Y <- tan(pi*u/2)
hist(Y, probability=T, col="brown",  main="Histogram and PDF")
lines(density(Y), col="blue",lty=2, lwd=2, xlab="Y")
mean(Y) #7.338874
#like in Ex3Q2c, the theoretical graph doesn't show up by adding the next few lines anyway
pdf <- function(x){
  B/(1+(x^2))
}
limits=seq(0.01,15,0.05) #15 is an arbitrary value since x goes up to infinity
lines(limits,pdf(limits), type='p', cex=0.8, col="red")
legend(10, 0.01, legend=c("Estimated PDF"),
       col=c("blue", "red"), lty=c(2,NA), pch=c(NA,1), cex=0.8)
#In this case, the function goes to 0 quickly when x takes small or large values, so the smoothing technique used in R creates the roughness seen in the plot.

###############################################################################################
###############################################################################################
#1c

#Method 1: sub y=exp(-x)
set.seed(123)
Nsim=10^5
u=runif(Nsim)
x1=1/(1+(log(u))^2)
Estimated_Int1=mean(x1)
A1 <- 1/Estimated_Int1 #1.61099

#Method 2: sub y=1/(x+1)
set.seed(123)
Nsim=10^5
u=runif(Nsim)
x2=(1/u^2)*exp(1-(1/u))/(1+((1/u)-1)^2)
Estimated_Int2=mean(x2)
A2 <- 1/Estimated_Int2 #1.612766

#actual value
f=function(y){
  exp(-y)/(1+(y^2))
}
integrate(f, 0, Inf)
real_mean <- 0.6214496
actual_a <- as.numeric(1/0.6214496) #1.609141

#Comparison for 1/A
#Method 1
diff1=(1/actual_a)-(1/A1) 
se1=sqrt(mean((x1)^2)-(mean(Estimated_Int1)^2))/sqrt(Nsim)
#Method 2
diff2=(1/actual_a)-(1/A2)     
se2=sqrt(mean((x2)^2)-(mean(Estimated_Int2)^2))/sqrt(Nsim)
#compare
print(rbind(diff1, diff2, se1, se2))

#Comparison for A
#Method 1
diff11=A1-actual_a
#Method 2
diff22=A2-actual_a
print(rbind(diff11, diff22))

###############################################################################################
###############################################################################################
#1d

#envelop g with B
set.seed(123)
u=runif(Nsim)
Y=tan(pi*u/2) #found from ITM
V=runif(Nsim)
#A/R
i.y <- which(V<exp(-Y))
out <- Y[i.y]
est_accept1 <- length(out)/Nsim #estimated acceptance rate #0.3962
Actual_Acceptance_rate=1/(actual_a*pi/2) #1/c
Actual_Acceptance_rate #0.3956271

actual_a*pi/2 #c

#check the value of c indeeds gives the maximum
f=function(x){
  fg <- actual_a*pi*exp(-x)/2
  return(fg)
}
plot(f, 0, 100, main="Plot of f(x)/g(x)", ylab="f(x)/g(x)")
optimize(f, lower=0, upper=100, maximum=TRUE)
#maximum: 4.627768e-05
#max point indeed at 2.527516

#envelop g with expo fn
set.seed(123)
Nsim <- 10^5
u=runif(Nsim)
Y2=log(1/u)
V=runif(Nsim)
#A/R
i.y2 <- which(V<1/(1+(Y2)^2))
out2 <- Y2[i.y2] 
est_accept2 <- length(out2)/Nsim #prob of acceptance
est_accept2 #0.62019

Actual_Acceptance_rate2=1/actual_a #1/c
Actual_Acceptance_rate2 #0.6214496

actual_a #c=1.609141

#check the value of c indeeds gives the maximum
f=function(x){
  fg <- actual_a/(1+x^2)
  return(fg)
}
plot(f, 0, 100, main="Plot of f(x)/g(x)", ylab="f(x)/g(x)")
optimize(f, lower=0, upper=100, maximum=TRUE)
#maximum: 4.627768e-05
#max point indeed at 1.609141

#Compare
print(rbind(est_accept1,est_accept2, Actual_Acceptance_rate, Actual_Acceptance_rate2))

###############################################################################################
###############################################################################################
#1e

Nsim=10^5
set.seed(123)
xxx=function(lambda){
  u=runif(Nsim)
  integrand=1/(1+(log(u)/lambda)^2)
  return(1/mean(integrand))
}
Alambda <- Vectorize(xxx) #vectorise and curve in order to get a continuous curve
curve(Alambda, 0,100, xlab="lambda", ylab="A")

###############################################################################################
###############################################################################################
#2a

set.seed(123)
Nsim <- 10^5
quant1 <- quantile(out2, 0.9)   #quantile at 90% probability level
quant1
quant2 <- quantile(out2, 0.99)   #quantile at 99% probability level
quant2
prob2a <- quant2 - quant1
prob2a #1.371777

###############################################################################################
###############################################################################################
#2b

#simulation
set.seed(123)
u=runif(Nsim)
sum_claims <- function(x){
  n <- rpois(1,20)
  x <- sample(out2, n, replace=FALSE, prob=NULL) #out2 is from Q1d
  s <- sum(x)
  return(s)
}
x <- 1:10^5
ss <- sapply(x,sum_claims)
means <- mean(ss) #11.07481
vars <- var(ss) #12.27502
hist(ss, probability = T, main="Histogram of S against density", col="lightpink", ylim = c(0,0.13))
lines(density(ss), col="aquamarine4", lty=2, lwd=2, xlab="S")
legend(25,0.06, legend=c("Estimated PDF"), col=c("aquamarine4"), lty=c(2,NA), pch=c(NA,1), cex=0.8)
theta <- 0.9
premium_0 <- (1+theta)*mean(ss)
premium_0 #21.04215

#quantiles by sample approximation
alpha.vec <- c(0.5,0.8, 0.9, 0.975, 0.99, 0.999)
q.ss <- array(sapply(alpha.vec, FUN=quantile, x=ss)) 
q.ss
#quantiles by Normal approximation
mu <- means
sigma <- sqrt(vars)
q.musig <- mu + sigma*qnorm(alpha.vec)
# table with the quantiles comparison
quantiles <- cbind(q.musig, q.ss)
rownames(quantiles) <- c("0.5","0.8", "0.9", "0.975", "0.99", "0.999")
colnames(quantiles) <- c("Normal", "Sample")
quantiles

###############################################################################################
###############################################################################################
#2c

######Method 1########
ptm = proc.time()
set.seed(123)
Nsim <- 10^5
#d = 1.371777
sr <- rep(0, Nsim) #create vector of 0's
si <- rep(0, Nsim)
for(i in 1:Nsim){
  n <- rpois(1,20) #generate N
  x <- sample(out2, n, replace=FALSE, prob=NULL) #out2 is from Q1d #generate rv X of size N
  xr <- (x-1.371777)*(sign(x-1.371777)+1)/2 #build indicator fn
  xi <- x-xr
  si[i] <- sum(xi)
  sr[i] <- sum(xr)
} 
#Insurer
mean(si) #10.142
var(si) #8.521025

#Reinsurer
mean(sr) #0.9325627 
var(sr) #1.196506 

SimulationTime_Method1=proc.time() - ptm
SimulationTime_Method1[3] #27.338 seconds to process

#Check: mean(si)+mean(sr)=mean(ss)

hist(sr, probability=T, col = 'palevioletred4', main="Density of Sum for Reinsurer")
lines(density(sr), col="blue", lty=2, lwd=2, xlab="S")
legend(6,0.6, legend=c("Estimated PDF"), col=c("blue"), lty=c(2,NA), pch=c(NA,1), cex=0.8)

theta <- 0.9
#Premium for Insurer
premium_I2 <- (1+theta)*mean(si)
premium_I2 # 19.26981

#Premium for Reinsurer
premium_R2 <- (1+theta)*mean(sr)
premium_R2 # 1.771868

#######Method 2########
ptm = proc.time()
set.seed(123)
Nsim <- 10^5
myf <- function(ssize){
  out <- rep(0, ssize) # create a vector of 0’s for x with desired size
  bb <- 0 # number of accepted values
  while (bb < ssize){
    u <- runif(1)
    v <- runif(1)
    Y2=log(1/u)
    if(v<1/(1+(Y2)^2)){ #A/R decision
      out[bb+1] <- Y2
      bb <- bb+1
    }
  }
  return(out)
}
#Generate the function required for x in the scmp.ri function below
# simulate SI and SR
scmp.ri <- function(ssize){
  n <- rpois(ssize,20)
  x <- sapply(n, FUN = myf)
  ss.r <- rep(0, length(x))
  ss.i <- rep(0, length(x))
  
  for(k in 1:length(x)){
    x.l <- x[[k]]
    i.ni <- which(x.l-1.371777 <= 0) # d= 1.371777
    i.nr <- which(x.l-1.371777>0)
    ss.r[k] <- sum(x.l[i.nr]-1.371777)
    ss.i[k] <- sum(x.l[i.ni])
  }
  return(cbind(ss.r, ss.i))
}
srsi.mat<- scmp.ri(Nsim) 
sr <- srsi.mat[,1] 
si <- srsi.mat[,2] 
# means
colSums (srsi.mat)/Nsim #0.9221592   #7.9904721
# variances
var(sr) #1.166307
var(si) #5.561647

SimulationTime_Method1=proc.time() - ptm
SimulationTime_Method1[3] # 23.007 seconds

hist(sr, probability=T, col = 'palevioletred4', main="Density of Sum for Reinsurer")
lines(density(sr), col="blue", lty=2, lwd=2, xlab="S")
legend(6,0.6, legend=c("Estimated PDF"), col=c("blue"), lty=c(2,NA), pch=c(NA,1), cex=0.8)

theta <- 0.9
#Premium for Insurer
premium_I2 <- (1+theta)*mean(si)
premium_I2 # 15.1819

#Premium for Reinsurer
premium_R2 <- (1+theta)*mean(sr)
premium_R2 # 1.752102

######normal approximation######
#From the histogram of S_R, it is evident that the normality assumption is not suitable for the 
#total claim of the reinsurer as there is a positive skew and has lack of symmetry.


###############################################################################################
###############################################################################################
#2d

set.seed(123)
Nsim <- 10^5
d <- prob2a
sr <- rep(0, Nsim) #create vector of 0's
si <- rep(0, Nsim)

#####Poisson rate 10#####
for(i in 1:Nsim){
  n <- rpois(1,10) #generate N
  x <- sample(out2, n, replace=FALSE, prob=NULL) #out2 is from Q1d #generate rv X of size N
  xr <- (x-d)*(sign(x-d)+1)/2 #build indicator fn
  xi <- x-xr
  si[i] <- sum(xi)
  sr[i] <- sum(xr)
} 
#Premium
theta <- 0.9
#Insurer
premium_I1 <- (1+theta)*mean(si)
premium_I1 #9.639715
#Reinsurer
premium_R1 <- (1+theta)*mean(sr)
premium_R1 #0.8909659

#####Poisson rate 30#####
for(i in 1:Nsim){
  n <- rpois(1,30) #generate N
  x <- sample(out2, n, replace=FALSE, prob=NULL) #out2 is from Q1d #generate rv X of size N
  xr <- (x-d)*(sign(x-d)+1)/2 #build indicator fn
  xi <- x-xr
  si[i] <- sum(xi)
  sr[i] <- sum(xr)
} 
#Premium
theta <- 0.9
#Insurer
premium_I3 <- (1+theta)*mean(si)
premium_I3 #28.89977
#Reinsurer
premium_R3 <- (1+theta)*mean(sr)
premium_R3 #2.662707

#####Poisson rate 40#####
for(i in 1:Nsim){
  n <- rpois(1,40) #generate N
  x <- sample(out2, n, replace=FALSE, prob=NULL) #out2 is from Q1d #generate rv X of size N
  xr <- (x-d)*(sign(x-d)+1)/2 #build indicator fn
  xi <- x-xr
  si[i] <- sum(xi)
  sr[i] <- sum(xr)
} 
#Premium
theta <- 0.9
#Insurer
premium_I4 <- (1+theta)*mean(si)
premium_I4 #38.57927
#Reinsurer
premium_R4 <- (1+theta)*mean(sr)
premium_R4 #3.5502

#####Poisson rate 50#####
for(i in 1:Nsim){
  n <- rpois(1,50) #generate N
  x <- sample(out2, n, replace=FALSE, prob=NULL) #out2 is from Q1d #generate rv X of size N
  xr <- (x-d)*(sign(x-d)+1)/2 #build indicator fn
  xi <- x-xr
  si[i] <- sum(xi)
  sr[i] <- sum(xr)
} 
#Premium
theta <- 0.9
#Insurer
premium_I5 <- (1+theta)*mean(si)
premium_I5 #48.14499
#Reinsurer
premium_R5 <- (1+theta)*mean(sr)
premium_R5 #4.438042

#Altogether
prem_insur <- cbind(premium_I1, premium_I2, premium_I3, premium_I4, premium_I5)
prem_reins <- cbind(premium_R1, premium_R2, premium_R3, premium_R4, premium_R5)
rbind(prem_insur, prem_reins)

###############################################################################################
###############################################################################################
#3

library(scatterplot3d)
set.seed(123)
n <- 10^5
w1 <- runif(n,-1,1)
w2 <- runif(n,-1,1)
w3 <- runif(n,-1,1)
s <- w1^2 + w2^2 + w3^2
i.c <- which(s<=1) #accepts values that fit the condition
n1 <- length(i.c) #the effective size (the num of non-zero values)
x <- w1[i.c]
y <- w2[i.c]
z <- w3[i.c]

acceptance_rate=n1/n #0.52365

library(rgl)
myColorRamp <- function(colors, values) {
  v <- (values - min(values))/diff(range(values))
  x <- colorRamp(colors)(v)
  rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
}
cols <- myColorRamp(c("red", "blue"), z) 
scatterplot3d(x, y, z, color=cols, angle = 20)
#P(fall inside)=P(area of circle/area of square)=(4/3)*pi*(1^3)/(2*2*2)=pi/6

###############################################################################################
###############################################################################################





