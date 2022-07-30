#CLEANING DATA
#Change variable class
dat <- read.csv("Market_data2019.csv", header=T)
head(dat)
names(dat)
summary(dat)
nrow(dat)
class(dat$ROE)  # ROE is wrongly viewed as a factor

# Converting variables which are wrongly recognised factors that should be numeric:
# ROE, EPS_Growth, Cost_of_Equity, CEO_holding, Institutional_holding
#strip out % and convert factor to numeric overwriting existing variable
dat$ROE <-  as.numeric(gsub("[\\%,]", "", dat$ROE))
dat$EPS_Growth <-  as.numeric(gsub("[\\%,]", "", dat$EPS_Growth))
dat$Cost_of_Equity <-  as.numeric(gsub("[\\%,]", "", dat$Cost_of_Equity))
dat$CEO_holding <-  as.numeric(gsub("[\\%,]", "", dat$CEO_holding))
dat$Institutional_holding <-  as.numeric(gsub("[\\%,]", "", dat$Institutional_holding))

class(dat$Region) #Region is factor
class(dat$Industry) #Industry is factor
class(dat$Number_of_firms) #Integer

write.csv(dat,'Market_data_clean.csv')

#Missing values
summary(dat)
#PBV has a missing value
which(is.na(dat$PBV)) #row 69:restaurants in PBV has NA value
dat <- na.omit(dat)

#take out Industry
dat$Industry<- NULL

#Initial Regression of all variables
reg1 <- lm(PE~.-PE, data=dat)
summary(reg1) #many unsignificant variables, adjusted R squared: 0.1288 without Industry


#VARIABLE SELECTION

#Check Multicollinearity
attach(dat)
#Take out Region and Industry because they must be numeric in order to see their CORRELATION
dat2 <- data.frame(cbind(PE, Number_of_firms, ROE, EPS_Growth, PBV, PS, Beta, Cost_of_Equity, CEO_holding, Institutional_holding))
detach(dat)
dat2cor <- cor(dat2) #there is high correlation among some predictors(EPS_GrowthvsROE) but does not necessarily imply there is severe multicollinearity.
corrplot(dat2cor, type="upper", order="hclust")

#TRY VIF
model_vif <- lm(log(PE)~.-Industry,data=dat)
car::vif(model_vif) #VIF for ROE, EPS_Growth and Cost_of_Equity exceed the cutoff of 10. This points to severe multicollinearity and makes sense because the calculations of EPS Growth and ROE are related, as well as the calcualtions for Beta vs Cost of Equity. Beta just narrowly missed the cutoff of VIF=10.
#To handle multicollinearity, we omit either Beta OR Cost of Equity, and EPS Growth OR ROE. See which one has higher adjusted R squared to see which one to omit.
reg_noroe <- lm(log(PE)~.-PE-Industry-ROE, data=dat)
summary(reg_noroe) #adjusted R squared is 0.1722 with Industry, 0.1307 without Industry
car::vif(reg_noroe)
reg_noeps <- lm(log(PE)~.-PE-Industry-EPS_Growth, data=dat)
summary(reg_noeps) #adjusted R squared is  0.1555 with Industry, 0.1171 without
car::vif(reg_noeps)
#No ROE has higher adjusted R squared than no EPS Growth, choose to omit ROE

reg_nocoe <- lm(log(PE)~.-PE-Industry-Cost_of_Equity, data=dat)
summary(reg_nocoe) #adjusted R squared is 0.172 with Industry, 0.1259 without
car::vif(reg_nocoe)
reg_nobeta <- lm(log(PE)~.-PE-Industry-Beta, data=dat)
summary(reg_nobeta) #adjusted R squared is 0.1675 with industry, 0.126 without
car::vif(reg_nobeta)
#No COE has higher adjusted R squared than no beta, choose to omit COE WITH industry, choose to omit Beta WITHOUT Industry

#New regression after handling multicollinearity
dat$ROE <- dat$Beta <- NULL
reg2 <- lm(log(PE)~.-PE, data=dat) #the variables are regionx2, no of firms, eps growth, pbv, ps, cost of equity, ceo holding, insti holding
summary(reg2) #adjusted R squared is 0.128 without Industry, which is slightly lower than reg1


#Best Subset Selection based on AIC criterion
model_bs<-leaps::regsubsets(log(PE)~.-PE,nvmax=11,data=dat)
plot(model_bs, scale="Cp") #best model: RegionEUR, EPS_Growth, PBV, PS, Institutional_holding

#Best Subset Selection based on BIC criterion
plot(model_bs, scale="bic") #best model: EPS_Growth, PBV

#Best Subset Selection based on R squared criterion
plot(model_bs, scale="adjr2") #best model: RegionEUR, EPS_Growth, PBV, PS, Institutional_holding

#R squared criterion fits with aic, bic took too little variables(penalised too much)

#Try Backward Elimination
model_back <-lm(log(PE)~.-PE,data=dat)
model_backb<-step(model_back, trace=0)
summary(model_backb) #model: EPS_Growth, PBV, PS, Institutional holding <& intercept RegionEMG??>
#significant variables: EPS_Growth, PBV, PS. adjusted R squared is 0.129

#NOTE: THESE ARE NOT NECESSARILY THE BEST MODELS, CHECK THE RESIDUALS THAT THE ASSUMPTIONS ARE SATISFIED AND CHECK FOR MULTICOLLINEARITY PROBLEMS IN THE COEFFICIENTS

#Try Forward Elimination, for this im not sure if the code is right, he didn't teach properly, dunno if should minus industry
null<-lm(log(PE)~1, data=dat)
full<-lm(log(PE)~.-PE,data=dat)
model_f<-step(null, scope=list(lower=null, upper=full),direction="forward", trace=0)
summary(model_f) #model: EPS Growth, PBV, PS, RegionEUR, RegionUS,
#significant variables are: EPS_Growth, PBV
# adjusted R squared is 0.1333

#stepwise regression, start from nothing
model_step1 <- step(lm(log(PE)~Region+EPS_Growth+PBV+PS+Institutional_holding #variables from AIC 
        ,data=dat),direction="both",k=2, trace=0)
summary(model_step1)

#stepwise regression, start from everything
model_step2 <- step(lm(log(PE)~1,data=dat), direction="both", k=2,scope=~Region+ #variables from AIC
       EPS_Growth+PBV+PS+Institutional_holding, trace=0)
summary(model_step2)

#New regression after variable selection, choose the variables from AIC
reg2 <- lm(log(PE)~Region+EPS_Growth+PBV+PS+Institutional_holding, data=dat) #the variables are regionx2, no of firms, eps growth, pbv, ps, cost of equity, ceo holding, insti holding
summary(reg2) #adjusted R squared is 0.1334 without Industry, which is higher than reg1


#TRANSFORMATION


#Visual scatterplot, this is to try to see if there is a different regression for the colour variable(which would be the dummy variable(which must be 0 or 1), Lab 3)
attach(dat)
par(mfrow=c(5,2))
plot(Industry, log(PE), xlab="Industry", ylab="PE", main="Industry boxplot")
plot(Region, log(PE), col ="red", xlab="Industry", ylab="PE", main="Industry boxplot")
plot(Number_of_firms, log(PE), xlab="Number of firms", ylab="PE", main="Number of firms scatterplot")
abline(lm(log(PE)~Number_of_firms),col="red")
plot(ROE, log(PE), xlab="ROE", ylab="PE", main="ROE scatterplot")
abline(lm(log(PE)~ROE),col="red")
plot(EPS_Growth, log(PE), xlab="EPS Growth", ylab="PE", main="EPS Growth scatterplot")
abline(lm(log(PE)~EPS_Growth),col="red")
plot(PBV, log(PE), xlab="PBV", ylab="PE", main="PBV scatterplot")
abline(lm(log(PE)~PBV),col="red")
plot(PS, log(PE), xlab="PS", ylab="PE", main="PS scatterplot")
abline(lm(log(PE)~PS),col="red")
plot(Beta, log(PE), xlab="Beta", ylab="PE", main="Beta scatterplot")
abline(lm(log(PE)~Beta),col="red")
plot(Cost_of_Equity, log(PE), xlab="Cost of Equity", ylab="PE", main="Cost of Equity scatterplot")
abline(lm(log(PE)~Cost_of_Equity),col="red")
plot(CEO_holding, log(PE), xlab="CEO holding", ylab="PE", main="CEO holding scatterplot")
abline(lm(log(PE)~CEO_holding),col="red")
plot(Institutional_holding, log(PE), xlab="Institutional holding", ylab="PE", main="Institutional holding scatterplot")
abline(lm(log(PE)~Institutional_holding),col="red")

# library(MASS)
# boxcox()


#Factor variables: Region (3 levels), Industry(not suitable)
levels(dat$Region) #Check reference group: EMG


#REGRESSION DIAGNOSTIC
#Partial residual plot
reg1 <- lm(PE~.-PE, data=dat)
car::crPlots(reg1, line=F, smooth=T) #detect var and covar problem using partial residual plot
reg2 <- lm(log(PE)~Region+EPS_Growth+PBV+PS+Institutional_holding, data=dat)
car::crPlots(reg2, line=F, smooth=T) #handle by transformation on y

#Standardised Residual Values
#Method 1: from Lab
# library(ggplot2)
# #plot of y vs x
# p1 <- ggplot(mSim,aes(x=x,y=y))+
#   geom_point()+
#   geom_abline(intercept=-19, slope=2.5, colour="red")
# p2<-ggplot(reg2, aes(x=x,y=.stdresid))+ #be careful here a .stdresid is from the re gression object (not the data)
#   geom_point()+
#   geom_hline(yintercept=2, col="red", linetype="dashed")+
#   geom_hline(yintercept=-2, col="red", linetype="dashed")+
#   labs(y="standardized residuals")
# gridExtra::grid.arrange(p1,p2,ncol=2)

#Method 2: from senior
#Leverage scatterplot
library(MASS)
sres <- stdres(reg2)
plot(sres, ylab="Standardised residuals", main="Standardised residuals scatterplot")
abline(h=3,col="red")
abline(h=-3,col="red")
#the residuals fall in a band(constant variance), there are a few outliers which could be influential points, they fall around zero and has no pattern
sres[sres>3]
sres[sres<-3]
#potential outliers: 3,81

#Hat Values
lev=hatvalues(reg2)
lev[lev>2*6/281] #2p/n
plot(lev, ylab="hat value", main="leverage scatterplot")
abline(h=2*6/281,col="red")


#But the influence of a point depends on the balance between the size of the outlier, residual size and its leverage. Use Cook's distance as a single measure of influence of a data point to the regression model.

#Influential points
#Cook's D
plot(reg2, which=4)
abline(h=4/281,col="blue")
abline(h=0.03,col="pink")
legend(0, 1, legend=c("4/281", "0.03"),
       col=c("blue", "pink"), lty=1, cex=0.8)

cd=cooks.distance(reg2)
cd[cd>4/281] #21 influential points, too many
cd[cd>0.03] #6 influential points

#remove outliers
reg2
summary(reg2) #adjusted r squared: 0.1933
reg3 <- lm(log(PE)~Region+EPS_Growth+PBV+PS+Institutional_holding, data=dat[-89,]) #61,71,81,89,157,263
reg3
summary(reg3) #adjusted r squared: 0.2115
plot(reg3,which=4)
#not very similar estimated coeffs for reg2(all points) and reg3(points minus outliers), so take out influential points

#how does this affect the adjusted R squared value?

#AVP
car::avPlots(reg2) #doesnt really tell much
#AVP can identify groups of influential points, and not just single points (Phil Chan)

#Normal QQ plot
library(ggplot2)
ggplot(reg1) +
  stat_qq(aes(sample = .stdresid)) +
  geom_abline()
#there is a departure from normality in the tails
reg3 <- lm(log(PE)~.-PE, data=dat) #transformed response variable shows the residuals more closely normally distributed
ggplot(reg3) +
  stat_qq(aes(sample = .stdresid)) +
  geom_abline()


#Linear Regression Assumptions
#Our multiple regression model has satisfy the NICE assumptions: Normality, Independence of Residuals, Constant Variance and zero Expectation of residual vector, and we have to verify that the model is linear.
plot(reg1,add.smooth=T)
#In the Residuals vs Fitted plot,the fitted(red) line is mostly horizontal, which data points randomly spaced around zero. We can further infer the constant variance assumption is satisfied and that the data has homoscedascity.
#<It currently shows a pattern, must fix this> Further, this plot does not reveal any patterns, so independence is satisfied.






