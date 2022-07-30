###############################################################################################
#1c
#Most efficient alpha
Nsim=10^5
set.seed(123)
B <- function(alpha){
  U <- runif(Nsim)
  Y <- -log(U)/alpha #found from ITM
  V <- runif(Nsim)
  i.y <- which(V<exp((-(Y-alpha)^2/2))) #A/R decision
  out <- Y[i.y]
  est_accept <- length(out)/Nsim #acceptance rate
  return(est_accept)
}
alpha <-c(0.5,1,2,3,5)
ss <- sapply(alpha, B)
ss #can see that alpha=1 has the highest acceptance rate, 0.75874
#alpha_star is =1

#Actual Acceptance rate for alpha star
c_alpha1 <- (1/1)*sqrt(2/pi)*exp(1^2/2)
Actual_Acceptance_rate=1/(c_alpha1) #1/c
Actual_Acceptance_rate #0.7601735

################################################################################################

#1d
set.seed(123)
#Z= S|Z|

#Simulate X=|Z|
alpha=1
myf <- function(ssize){
  out <- rep(0, ssize) # create a vector of 0’s for x with desired size
  bb <- 0 # number of accepted values
  while (bb < ssize){
    u <- runif(1)
    v <- runif(1)
    Y2= -log(u)/1
    if(v<exp(-(Y2-1)^2)/2){ #A/R decision 
      out[bb+1] <- Y2
      bb <- bb+1
    }
  }
  return(out)
}
fixedn <- myf(10^5)

#Simulate S
#Step 1: Generate a rv U~U(0,1)
u <- runif(Nsim)
s <- rep(0, Nsim) #store rv S
#Step 2: Set Z=|Z| if U ≤ 0.5
s1 <- which(u<=1/2) #if uniform rv less than 0.5, it is S=-1
s[s1] <- 1
#Step 3: Set Z=-|Z| if U > 0.5
s2 <- which(u>1/2)
s[s2] <- -1
#Simulate Z
Z <- s*fixedn

################################################################################################

#1e

set.seed(123)
#X
hist(fixedn, col='mediumseagreen',probability=T, main="Histogram of X")
lines(density(fixedn),col="blue",lty=2, lwd=2, xlab="X")
legend(2, 0.4, legend=c("Estimated PDF") , col=c("aquamarine4") , lty=c(2,NA), pch=c(NA,1), cex=0.8)
#we can see it only takes the positive values of a standard normal distribution; it is the folded distribution

#Z
hist(Z, probability=T,col="lightblue", main="Histogram of Z")
lines(density(Z),col="blue",lty=2, lwd=2, xlab="Z")
legend(2, 0.2, legend=c("Estimated PDF") , col=c("aquamarine4") , lty=c(2,NA), pch=c(NA,1), cex=0.8)
#algorithm for Z cannot be improved by changing the rate of the exponential; that is, 
#by using an exponential density g(x) = λe−λx for some because λ = 1 minimizes
#the value of c obtained as proved by the calculus in 1c)).

################################################################################################

#1f

#acceptance rate for alpha star is 0.75874 (from 1c)

#enter explanation for polar method here

#Polar Method
n <- 10^5
w1 <- runif(n,-1,1)
w2 <- runif(n,-1,1)
s <- w1^2 + w2^2
i.c <- which(s<=1) #accepts values that fit the condition
n1 <- length(i.c) #the effective size (the num of non-zero values)
x <- w1[i.c]
y <- w2[i.c]
plot(x,y, col='royalblue4')
n1/n #0.78593
#acceptance rate for polar method 
#P(fall inside)=P(area of circle/area of square)=(pi(1^2))/4/(1*1)=pi/4~~ 0.785

#polar method is slightly better

################################################################################################

#2a

set.seed(123)
Nsim <- 10^6
#-log(U)/alpha generates a geometric rv
c <- sample(c(2,3),10^6, replace=T, prob=c(0.25,0.75))
c1 <- runif(Nsim)
x <- -log(c1)/c
hist(x, col='bisque',probability=T, main="Histogram of X")
lines(density(x),col="blue",lty=2, lwd=2, xlab="X")
legend(2, 0.4, legend=c("Estimated PDF") , col=c("aquamarine4") , lty=c(2,NA), pch=c(NA,1), cex=0.8)
mean(x) #0.3743057

################################################################################################

#2b

#Actual Value
integrand=function(x){
  ((0.5*exp(-2*x))+(2.25*exp(-3*x)))*x               
}
real_value=integrate(integrand, lower=0, upper=Inf)
real_value=real_value[[1]] #0.374999
###########################################################
#Sub y=e^(-x)
set.seed(123)
#Method 1: standard MC (without variance reduction)
Nsim=10^6
y=runif(Nsim)
x1_MC= ((0.5*y)+(2.25*y^2))*(-log(y))  #transformation to change integral limits to 0 and 1                 
Mu_MC=mean(x1_MC)
Mu_MC #0.3747709
SE_MC=sd(x1_MC) 
SE_MC #0.170942
var(x1_MC) #0.02922116

#Method 2 + 3: standard MC with SS (No AV)
#Step 1: Compute SS without AV
#Step 1.1: Generate a sample of 10^3 of U(0,1) rv
set.seed(123)
NSS <- 10^4 #smaller sample size than MC
u=runif(NSS) 

#Step 1.2: Generate j from 1:10^4
j=1:NSS

#Step 1.3: Generate x1 using the formula computed for SS
x1_SS <- ((0.5*((u+j-1)/NSS))+(2.25*((u+j-1)/NSS)^2))*(-log(((u+j-1)/NSS)))  

#Step 1.4: Find the mean
Mu_SS=mean(x1_SS) 
Mu_SS #0.3750001
SE_SS=sd(x1_SS)
SE_SS #0.1709063
var(x1_SS) #0.02920897

#compare these methods
compareSS <- cbind(real_value, Mu_MC, Mu_SS, SE_MC, SE_SS) 
compareSS
#Comments:
#All 2 methods means are close to the real value
#Note: the SS methods gets the value with SMALLER sample size

###########################################################

#try sub y=1/(x+1)
#Standard MC (without variance reduction)
Nsim=10^6  
set.seed(123)
y=runif(Nsim)
x2_MC = ((0.5*exp(-2*((1/y)-1)))+(2.25*exp(-3*((1/y)-1))))*((1-y)/(y^3))                    
Mu2_MC=mean(x2_MC)
Mu2_MC #0.3745244
SE2_MC=sd(x2_MC) #0.3010737
SE2_MC #0.3010737

#MC with SS
#Step 1.1: Generate a sample of 10^4 of U(0,1) rv
nss <- 10^4
u=runif(nss) #note the sample size is smaller

#Step 1.2: Generate j from 1:10^4
j=1:nss

#Step 1.3: Generate x1 using the formula computed for SS
x2_SS <- ((0.5*exp(-2*((1/((u+j-1)/nss))-1)))+(2.25*exp(-3*((1/((u+j-1)/nss))-1))))*((1-((u+j-1)/nss))/(((u+j-1)/nss)^3))  

#Step 1.4: Find the mean
Mu2_SS=mean(x2_SS) 
Mu2_SS #0.3750005
SE2_SS=sd(x2_SS)
SE2_SS #0.3010472

###########################################################
#What size should your sample be

#standard MC
set.seed(123)
n=100*(10^6)
y=runif(n)
x3_MC= ((0.5*y)+(2.25*y^2))*(-log(y))                
SE_MC=sd(x3_MC) 
SE_MC #0.1708998
#Variance under MC with sample size 100 times the original
var(x3_MC) #0.02920676

#SS
u=runif(10^4) 
#Generate j from 1:10^4
j=1:10^4

x3_SS <- ((0.5*((u+j-1)/10^4))+(2.25*((u+j-1)/10^4)^2))*(-log(((u+j-1)/10^4)))  

SE_SS=sd(x3_SS)
SE_SS #0.1709063
var(x3_SS) #0.02920897

################################################################################################

#2c

Nsim=10^6  
set.seed(123)
y=runif(Nsim)
#Verify function h is monotone over interval (0,1)
plot(y,((0.5*y)+(2.25*y^2))*(-log(y)), ylab="h", col='dodgerblue2', main="Plot of h")
#Function h is NOT monotone over the interval (0,1) so we cannot apply the method combined with antithetic variables.

################################################################################################

#3a

#T=20
set.seed(123)
rnf3=function(x){
  nn=rpois(1,20*90) #Step 1.1
  v=20*runif(nn) #Step 1.2
  vv=sort(v) #Step 1.3
  x <- mix.rv(nn) #Step 2
  xx=cumsum(x) 
  c = (1+ 2/3)*90*1*(3/8)
  ss=xx-(c*vv) #Step 3
  s=max(ss) #Step 4: max
  return(s)
}
#Step 5: Simulate a sample of 10^3
Zzz=sapply(1:10^6 ,rnf3)

#Step 6: Compute P(ruin) when u=0
P1=length(Zzz[Zzz>0])/10^6   # The probability of ruin for u=0
#Step 7: Compute P(ruin) when u=5
P2=length(Zzz[Zzz>5])/10^6  # The probability of ruin for u=5
#Step 8: Compute P(ruin) when u=10
P3=length(Zzz[Zzz>10])/10^6  # The probability of ruin for u=10

#probability of ruin for different initial capital
cbind(P1,P2,P3) #0.599833 0.003812 2.3e-05
#see that as the initial capital increases, the probability of ruin decreases

yr3=1:10^3
for (i in 1:10^3){
  yr3[i]=length(Zzz[Zzz>i/10])
}
xr3=1:10^3 
xr3=xr3/10 #initial capital, u
yr3=yr3/10^6 #ruin probability
plot(xr3,yr3,type="l", col="red", xlab="u", ylab="  Probability",
             main="Ruin Probability for T=20")

#T=50
rnf=function(x){
  nn=rpois(1,50*90) #Step 1.1
  v=50*runif(nn) #Step 1.2
  vv=sort(v) #Step 1.3
  x <- mix.rv(nn) #Step 2
  xx=cumsum(x) 
  c = (1+ 2/3)*90*1*(3/8)
  ss=xx-(c*vv) #Step 3
  s=max(ss) #Step 4: max
  return(s)
}

#Step 5: Simulate a sample of 10^6
Z=sapply(1:10^6 ,rnf)

#Step 6: Compute P(ruin) when u=0
P1=length(Z[Z>0])/10^6   # The probability of ruin for u=0
#Step 7: Compute P(ruin) when u=5
P2=length(Z[Z>5])/10^6  # The probability of ruin for u=5
#Step 8: Compute P(ruin) when u=10
P3=length(Z[Z>10])/10^6  # The probability of ruin for u=10

#probability of ruin for different initial capital
cbind(P1,P2,P3) #0.599573 0.003871 2.6e-05
#see that as the initial capital increases, the probability of ruin decreases

yr=1:10^3
for (i in 1:10^3){
  yr[i]=length(Z[Z>i/10])
}
xr=1:10^3 
xr=xr/10 #initial capital, u
yr=yr/10^6 #ruin probability
plot(xr,yr,type="l", col="red", xlab="u", ylab="  Probability",
     main="Ruin Probability for T=50")

#T=100
rnf2=function(x){
  nn=rpois(1,100*90) #Step 1.1
  v=100*runif(nn) #Step 1.2
  vv=sort(v) #Step 1.3
  x <- mix.rv(nn) #Step 2
  xx=cumsum(x) 
  c = (1+ 2/3)*90*1*(3/8)
  ss=xx-(c*vv) #Step 3
  s=max(ss) #Step 4: max
  return(s)
}
#Step 5: Simulate a sample of 10^3
Zz=sapply(1:10^6 ,rnf2)

#Step 6: Compute P(ruin) when u=0
P1=length(Zz[Zz>0])/10^6   # The probability of ruin for u=0
#Step 7: Compute P(ruin) when u=5
P2=length(Zz[Zz>5])/10^6  # The probability of ruin for u=5
#Step 8: Compute P(ruin) when u=10
P3=length(Zz[Zz>10])/10^6  # The probability of ruin for u=10

#probability of ruin for different initial capital
cbind(P1,P2,P3) #0.600586 0.003823 2.4e-05
#see that as the initial capital increases, the probability of ruin decreases

yr2=1:10^3
for (i in 1:10^3){
  yr2[i]=length(Zz[Zz>i/10])
}
xr2=1:10^3 
xr2=xr2/10 #initial capital, u
yr2=yr2/10^6 #ruin probability
plot(xr2,yr2,type="l", col="red", xlab="u", ylab="  Probability",
              main="Ruin Probability for T=100")

################################################################################################

#3b

#function r(u)
rufn <- function(u){
  r_u <- ((1/35)*exp(-(12/5)*u))+((4/7)*exp(-u))
}
u <- 0:10
ru <- sapply(u, rufn)
plot(u, ru, type="l",xlab="u",ylab="r(u)")

################################################################################################

#3c

#from T=100
#choose u=0
set.seed(156)
#Surplus process
nn=rpois(1,100*90) #T=time horizon , lambda=90
v=100*runif(nn) #arrival times
vv=sort(v) #t
x <- mix.rv(nn)
xx=cumsum(x) #S(t)
c = (1+2/3)*90*(3/8)
u_t <- (c*vv)-xx #u=0
# ss=xx-(225/4)*vv #c = (1+2/3)*lambda*t*E(X) = 225/4
# s=max(ss)
plot(vv, u_t, type="s", xlab="t", ylab = "Surplus Process, U(t)", xlim = c(0, 0.2), ylim=c(-1,3), main = "Surplus process when initial capital is 0")

################################################################################################

#4b

#####Find the values of k for the 3 quantiles####
set.seed(123)
n <- 10^6
quant_S_sim <- function(size){ #simulate S(2)
  N <- rpois(size, 2*30) #t=2, lambda = 30
  X <- sapply(N, FUN=rexpwunif, beta=15)
  S <- sapply(X, FUN=sum)
  return(S)
}
quant_S <- quant_S_sim(n)
k <- quantile(quant_S, probs=c(0.9, 0.95, 0.99)) #quantiles of S(2)
k #4.955222 5.254340 5.844909

#######Estimation without Importance Sampling#########
set.seed(123)
#S1: Build function to simulate from an exponential with the inverse method
rexpwunif <- function (ns, beta){
  unifs <- runif(ns)
  aux <- -(1/beta)*log(unifs) #from ITM
  aux
}
#S2: Build function for S(2)
scmp <- function(ssize){                #S(2)
  N <- rpois(ssize,60)                    #N~Pois(30*2)
  x <- sapply(N,FUN = rexpwunif, beta=15) #X~Exp(beta=15)
  ss <- sapply(x, FUN = sum)              #sum of X
  return(ss)
}
#S3: Get estimation for S(2)
ss <- scmp(n)
mean(ss) # 3.999528 #should be 60*(1/15)=4
var(ss) # 0.533375

#S4: Build indicator (S(2)−k)+
#90th percentile, k=4.955222
i.ss_90 <- which(ss-4.955222>=0)       
ssk_90 <- rep(0, n) #n=10^6
ssk_90[i.ss_90] <- (ss[i.ss_90] - 4.955222)
mean(ssk_90) #0.04005104
var(ssk_90) # 0.02710562

#95th percentile, k=5.254340
i.ss_95 <- which(ss-5.254340>=0)       
ssk_95 <- rep(0, n) #n=10^6
ssk_95[i.ss_95] <- (ss[i.ss_95] - 5.254340)
mean(ssk_95) #0.01819638
var(ssk_95) # 0.01169732

#99th percentile, k=5.844909
i.ss_99 <- which(ss-5.844909>=0)       
ssk_99 <- rep(0, n) #n=10^6
ssk_99[i.ss_99] <- (ss[i.ss_99] - 5.844909)
mean(ssk_99) #0.003053425
var(ssk_99) # 0.001757462

######Importance Sampling Estimation##########

set.seed(123)
n <- 10^6 ; beta <- 15 ; lambda <- 30 ; T <- 2
tt <- function(k){
  t1 <- 15+ sqrt((T*15*30)/k)
  t2 <- 15- sqrt((T*15*30)/k)
  return(cbind(t1,t2))
}
t <- sapply(k, tt)
t <- t[2,] #t2 is consistent with calculations
t #t[1]=1.523110, t[2]= 1.912335, t[3]= 2.591126

#S1: Simulate S∗ with t_i = beta -(30/sqrt(k_i))
impsamp <- rep(0, 3)
muu <- rep(0,3)
varr <- rep(0,3)
for(i in 1:3){
  scmp.t <- function(ssize){
    N <- rpois(ssize, (2*15*30)/(15-t[i]))            #N`~Pois(T*beta*lambda/(beta-t_i))
    x <- sapply(N,FUN = rexpwunif , beta = 15-t[i])   #X`~Exp(beta - t_i)
    ss <- sapply(x, FUN = sum)
    return(ss)
  }
  #S2: Simulate S∗_(2)
  ss.t <- scmp.t(n)
  
  #S3: Build indicator (S∗_(2) − k)+
  i.sst <- which(ss.t-k[i]>=0)       
  ssk.t <- rep(0, n)
  ssk.t[i.sst] <- ss.t[i.sst] - k[i]
  
  #IS estimator
  impsamp <- exp(-(t[i])*ss.t) * ssk.t * exp(2*30*((15/(15- t[i]))-1))
  muu[i] <- mean(impsamp)
  varr[i] <- var(impsamp)
}

is_stats <- cbind(muu, varr)
is_stats                    
#it obtained similar values to the MC estimator with significantly less variance

################################################################################################


