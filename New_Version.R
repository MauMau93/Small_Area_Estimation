##################################################################################

#Integrating Census and Survey Information for Small Area Estimation

##################################################################################



rm(list=ls())
#setwd("C:/Isa/AdvisorWB/AdvisorBootstrapPaulACorral/Simulations/NewBootstrap_Modified")

library(nlme)
library(sae)

#########################################
# Generation of a population
#########################################

# Setting of sizes and initialization of vectors

N<-20000        # Number of population units
D<-80           # Number of areas
Nd<-rep(250,D)      # Number of population units in each area
p<-3            # Number of auxiliary variables (intercept+2 dummies)
p1<-rep(0,D)       # Probability of the first category of X1
x1<-NULL       # X1 values
p2<-rep(0,D)       # Probability of the second category of X2
x2<-NULL       # X2 values
#Generation of auxiliary variables (2 dummies with probability of 1 given by p1 and p2)

for (d in 1:D){
  p1[d]<-0.3+0.5*d/D
  p2[d]<-0.2*d/D
  x1<-c(x1,rbinom(Nd[d],size=1,prob=p1[d]))
  x2<-c(x2,rbinom(Nd[d],size=1,prob=p2[d]))
}

write.table(data.frame(x1,x2),"AuxVarD80.txt")
x<-read.table("AuxVarD80.txt",header=TRUE)
x1<-x[,1]
x2<-x[,2]

#########################################
# Extraction of a stratified SRS
#########################################

# Sample sizes

nd<-Nd/5
n<-sum(nd)

# Extraction of sample indexes

samp<-NULL
sumNd<-0


for (d in 1:D){
  samp<-c(samp,sample((sumNd+1):(sumNd+Nd[d]),size=nd[d],replace=FALSE))
  sumNd<-sumNd+Nd[d]
}

write.table(samp,"SampleD80.txt")
sample<-read.table("SampleD80.txt",header=TRUE)
samp<-sample[,1]

#Sample elements of aux. variables
x1s<-x1[samp]
x2s<-x2[samp]
Xs<-cbind(rep(1,n),x1s,x2s)

# Non-Sample elements of aux. variables
x1r<-x1[-samp]
x2r<-x2[-samp]
Xr<-cbind(rep(1,N-n),x1r,x2r)

# Non-sample sizes of the areas
rd<-Nd-nd

# Area/province indicator

area<-NULL
for(d in 1:D){
  area<-c(area,rep(d,nd[d]))
  areaP<-NULL #*** New
  areaP<-c(areaP,rep(d,Nd[d]))
}

weights <- Nd/nd
out_of_samp <- (Nd/nd) -1

# Augmented Auxiliary Variables (aav)
domain <- rep(c(1:D), each=Nd[1]) 
x1_aug <- rep(x1s, each= weights)
x2_aug <- rep(x2s, each= weights)
X_aug<-cbind(domain,x1_aug,x2_aug)

write.table(data.frame(X_aug),"AuxVarAug.txt")
X_aug<-read.table("AuxVarAug.txt",header=TRUE)
X_aug$domain <- as.factor(X_aug$domain) #*** I think X_aug is not used any more

# Out-of sample Augmented Auxiliary Variables 
domain_out <- rep(c(1:D), each=rd[1]) 
x1_o_aug <- rep(x1s, each= out_of_samp)
x2_o_aug <- rep(x2s, each= out_of_samp)
X_o_aug<-cbind(domain_out,x1_o_aug,x2_o_aug)

write.table(data.frame(X_o_aug),"AuxVarAug_out.txt")
X_o_aug<-read.table("AuxVarAug_out.txt",header=TRUE)
X_o_aug$domain_out <- as.factor(X_o_aug$domain_out)

#*** For unit-context CEB (UC-CEB)

#x1.mean<-as.vector(tapply(x1,areaP,mean))
#x2.mean<-as.vector(tapply(x1,areaP,mean))

#X1s.mean<-NULL
#X2s.mean<-NULL
#for(d in 1:D){
#  X1s.mean<-c(X1s.mean,rep(x1.mean[d],nd[d]))
#  X2s.mean<-c(X2s.mean,rep(x2.mean[d],nd[d]))
#}

# Function for indicators: Poverty Incidence and Poverty Gap

# Alpha = 0
povertyincidence <- function(y) {
  result <- mean(y < z) 
  return (result)}

# Alpha = 1
povertygap <- function(y) {
  result <- mean((y<z) * (z-y) / z)
  return (result)}

#########################################
# Monte Carlo simulations
#########################################

# Parameter values

beta<-c(3,0.03,-0.04)
sigmae2<-(0.5)^2
sigmau2<-(0.15)^2
sigmau2/(sigmau2+sigmae2/nd[1])
z<-12    # Poverty line: It is approximately 0.6 Median(E) for the whole population

# gamma=sigmau2/(sigmau2+sigmae2/nd[1])=0.81
# We should try with different values for gamma, such as gamma=0.3,0.6

# Monte Carlo iterations

#nMC<-500
nMC<-100

# Initialization of matrices that will contain poverty proportions and gaps
# Each column contains the D area quantities for each Monte Carlo simulation

propMC<-matrix(0,nr=D,nc=nMC)
gapMC<-matrix(0,nr=D,nc=nMC)
propMC.samp <- matrix(0,nr=D,nc=nMC)
gapMC.samp <- matrix(0,nr=D,nc=nMC)

propfinMC<-matrix(0,nr=D,nc=nMC)
gapfinMC<-matrix(0,nr=D,nc=nMC)
propfinMC.ELL<-matrix(0,nr=D,nc=nMC)
gapfinMC.ELL<-matrix(0,nr=D,nc=nMC)
propfinMC.S<-matrix(0,nr=D,nc=nMC)
gapfinMC.S<-matrix(0,nr=D,nc=nMC)
propfinMC.uc <- matrix(0,nr=D,nc=nMC)
gapfinMC.uc <- matrix(0,nr=D,nc=nMC)
propfinMC.New <- matrix(0,nr=D,nc=nMC)
gapfinMC.New <- matrix(0,nr=D,nc=nMC)


MSEpropMC.B<-matrix(0,nr=D,nc=nMC)
MSEgapMC.B<-matrix(0,nr=D,nc=nMC)
MSEpropMC.ELL<-matrix(0,nr=D,nc=nMC)
MSEgapMC.ELL<-matrix(0,nr=D,nc=nMC)
MSEpropMC.S<-matrix(0,nr=D,nc=nMC)
MSEgapMC.S<-matrix(0,nr=D,nc=nMC)
MSEpropMC.New <- matrix(0,nr=D,nc=nMC)
MSEgapMC.New <- matrix(0,nr=D,nc=nMC)

time1<-Sys.time()
time1

for (i in 1:nMC){  # Start of Monte Carlo cycle: Generation of populations and calculation of EBPs
  
  #i<-1
  #print(i)
  #print(Sys.time())
  
  sumNd<-0
  sumnd<-0
  y<-NULL
  ys<-NULL
  prop<-rep(0,D)
  gap<-rep(0,D)
  prop.samp <- rep(0,D)
  gap.samp <- rep(0,D)
  
  # Generation of population values of the response, area by area
  # Calculation of true area poverty proportions and gaps
  
  for (d in 1:D){
    
    ed<-rnorm(Nd[d],0,sqrt(sigmae2))
    ud<-rnorm(1,0,sqrt(sigmau2))
    Xd<-cbind(rep(1,Nd[d]),x1[(sumNd+1):(sumNd+Nd[d])],x2[(sumNd+1):(sumNd+Nd[d])])
    
    mud<-Xd%*%beta
    yd<-mud+ud+ed
    y<-c(y,yd)
    
    Ed<-exp(yd)
    prop[d]<-mean(Ed<z)
    gap[d]<-mean((Ed<z)*(z-Ed)/z)
    
    # I save them in the columns of a matrix
    propMC[d,i]<-prop[d]
    gapMC[d,i]<-gap[d]
    
    # Generate sample values
    
    Xsd<-Xs[(sumnd+1):(sumnd+nd[d]),]
    musd<-Xsd%*%beta
    esd<-rnorm(nd[d],0,sqrt(sigmae2))
    ysd<-musd+ud+esd
    ys<-c(ys,ysd)
    
    sumNd<-sumNd+Nd[d]
    sumnd<-sumnd+nd[d]
    
  }
  
  E<-exp(y)       # People earnings
  
  # Sample data vectors and matrices
  
  ys<-y[samp]
  Es<-E[samp]
  Es<-exp(ys)
  
  prop.samp[d]<-mean(Es<z)
  gap.samp[d]<-mean((Es<z)*(z-Es)/z)
  
  # I save them in the columns of a matrix
  
  propMC.samp[d,i]<-prop.samp[d]
  gapMC.samp[d,i]<-gapMC.samp[d]
  
  # Fitting of nested-error model to sample data
  
  fit<-lme(ys~x1s+x2s,random=~1|area,method="REML")
  betaest<-fixed.effects(fit)
  upred<-random.effects(fit)
  sigmae2est<-fit$sigma^2 # Varianza residual
  sigmau2est<-as.numeric(VarCorr(fit)[1,1]) # Matriz de covarianzas de las componentes aleatorias del modelo
  Xs<-model.matrix(fit)
  
  #####################################################################    
  # Calculating Census EB, Survey EB, ELL With Analytical Formulas
  #####################################################################  
  
  propfin<-rep(0,D)
  gapfin<-rep(0,D)
  propfin.ELL<-rep(0,D)
  gapfin.ELL<-rep(0,D)
  propfin.S<-rep(0,D)
  gapfin.S<-rep(0,D)
  propfin.uc<-rep(0,D)
  gapfin.uc<-rep(0,D)
  propfin.New <- rep(0,D)
  gapfin.New <- rep(0,D)
  
  sumNd<-0
  sumnd<-0
  
  for (d in 1:D){
    Xd<-cbind(rep(1,Nd[d]),x1[(sumNd+1):(sumNd+Nd[d])],x2[(sumNd+1):(sumNd+Nd[d])])
    
    # Census EB
    mudcond<-Xd%*%matrix(betaest,nr=p,nc=1)+upred[d,1]
    gammad<-sigmau2est/(sigmau2est+sigmae2est/nd[d])
    sigma2cond<-sigmau2est*(1-gammad)+sigmae2est
    alphad<-(log(z)-mudcond)/sigma2cond  ## alphad is a vector of size Nd[d] because Xd is a matrix
    propfin[d]<-mean(pnorm(alphad))
    gapfin[d]<-mean( pnorm(alphad)*(1-(exp(mudcond+sigma2cond/2)*pnorm(alphad-sqrt(sigma2cond))/pnorm(alphad) )/z)  )
    
    # ELL
    mud.ELL<-Xd%*%matrix(betaest,nr=p,nc=1)
    sigma2.ELL<-sigmau2est+sigmae2est
    alphad.ELL<-(log(z)-mud.ELL)/sigma2.ELL  ## alphad is a vector of size Nd[d] because Xd is a matrix
    propfin.ELL[d]<-mean(pnorm(alphad.ELL))
    gapfin.ELL[d]<-mean( pnorm(alphad.ELL)*(1-(exp(mud.ELL+sigma2.ELL/2)*pnorm(alphad.ELL-sqrt(sigma2.ELL))/pnorm(alphad.ELL) )/z)  )
    
    
    
    #*** Survey EB (Census EB, based only on the sample values)
    Xsd.S<-Xs[(sumnd+1):(sumnd+nd[d]),]
    mucond.S<-Xsd.S%*%matrix(betaest,nr=p,nc=1)+upred[d,1]
    sigma2cond.S <- sigma2cond 
    alphad.S<-(log(z)-mucond.S)/sigma2cond.S   # Now alphad.S is a vector of size nd[d]
    propfin.S[d]<-mean(pnorm(alphad.S))  # Now is the sample mean (equal to the weighted sample mean with weights Nd[d]/nd[d])
    gapfin.S[d]<-mean( pnorm(alphad.S)*(1-(exp(mucond.S+sigma2cond/2)*pnorm(alphad.S-sqrt(sigma2cond))/pnorm(alphad.S) )/z)  )
    
    # Census EB
    mudcond<-Xd%*%matrix(betaest,nr=p,nc=1)+upred[d,1]
    gammad<-sigmau2est/(sigmau2est+sigmae2est/nd[d])
    sigma2cond<-sigmau2est*(1-gammad)+sigmae2est
    alphad<-(log(z)-mudcond)/sigma2cond  ## alphad is a vector of size Nd[d] because Xd is a matrix
    propfin[d]<-mean(pnorm(alphad))
    gapfin[d]<-mean( pnorm(alphad)*(1-(exp(mudcond+sigma2cond/2)*pnorm(alphad-sqrt(sigma2cond))/pnorm(alphad) )/z)  )
    
    
    
    
    #*** UC-CEB
    
    #Xsd.uc<-matrix(c(1,x1.mean[d],x2.mean[d]),nr=1,nc=p)
    #mucond.uc<-Xsd.uc%*%matrix(betaest,nr=p,nc=1)+upred.uc[d,1]
    #gammad.uc<-sigmau2est.uc/(sigmau2est.uc+sigmae2est.uc/nd[d])
    #sigma2cond.uc<-sigmau2est.uc*(1-gammad.uc)+sigmae2est.uc
    #alphad.uc<-(log(z)-mudcond.uc)/sigma2cond.uc   # Now alphad.uc is a scalar
    #propfin.uc[d]<-pnorm(alphad.uc)  # Now all the elements in the sum for the area are equal
    #gapfin.uc[d]<-pnorm(alphad.uc)*(1-(exp(mudcond.uc+sigma2cond.uc/2)*pnorm(alphad.uc-sqrt(sigma2cond.uc))/pnorm(alphad.uc) )/z)  
    
    sumNd<-sumNd+Nd[d]
    sumnd<-sumnd+nd[d]
  }

  
  ############################################################
  # Parametric Bootstrap
  # for the calculation of the MSE of the Survey EB
  ############################################################
  
  #B<-500
  B<-50
  
  MSEpropsum.B<-rep(0,D)
  MSEgapsum.B<-rep(0,D)
  
  for (b in 1:B){   ### Start of bootstrap cycle
    
    print(cbind("Bootstrap iteration",b,"Monte Carlo iteration",i))
    
    y.B<-NULL
    ys.B<-NULL
    ys.B2<-NULL
    prop.B<-rep(0,D)
    gap.B<-rep(0,D)
    
    sumnd<-0
    sumrd<-0
    sumNd<-0
    
    for (d in 1:D){
      
      # Generation of population values of y from the fitted model
      # Parametric bootstrap for MSE estimation
      
      ed.B<-rnorm(Nd[d],0,sqrt(sigmae2est))
      ud.B<-rnorm(1,0,sqrt(sigmau2est))
      
      Xd<-cbind(rep(1,Nd[d]),x1[(sumNd+1):(sumNd+Nd[d])],x2[(sumNd+1):(sumNd+Nd[d])])
      mud.B<-Xd%*%matrix(betaest,nr=p,nc=1)
      yd.B<-mud.B+ud.B+ed.B
      y.B<-c(y.B,yd.B)
      
      # True poverty measures for the bootstrap population
      
      Ed.B<-exp(yd.B)
      prop.B[d]<-mean(Ed.B<z)
      gap.B[d]<-mean((Ed.B<z)*(z-Ed.B)/z)
      
      # Generate bootstrap sample values
      
      Xsd<-Xs[(sumnd+1):(sumnd+nd[d]),]
      musd.B<-Xsd%*%matrix(betaest,nr=p,nc=1)
      esd.B<-rnorm(nd[d],0,sqrt(sigmae2est))
      ysd.B<-musd.B+ud.B+esd.B
      ys.B<-c(ys.B,ysd.B)
      
      sumnd<-sumnd+nd[d]
      sumrd<-sumrd+rd[d]
      sumNd<-sumNd+Nd[d]
      
    }
    
    E.B<-exp(y.B)
    
    # Bootstrap sample data
    
    Es.B<-exp(ys.B)
    
    # Fitting of nested-error model to bootstrap sample data
    
    fit.B<-lme(ys.B~x1s+x2s,random=~1|area,method="REML")
    betaest.B<-fixed.effects(fit.B)
    upred.B<-random.effects(fit.B)
    sigmae2est.B<-fit.B$sigma^2    # Varianza residual
    sigmau2est.B<-as.numeric(VarCorr(fit.B)[1,1]) # Matriz de covarianzas de las componentes aleatorias del modelo
    
    propfin.B<-rep(0,D)
    gapfin.B<-rep(0,D)
    
    sumNd<-0
    sumnd<-0
    
    for (d in 1:D){
      #Xd<-cbind(rep(1,Nd[d]),x1[(sumNd+1):(sumNd+Nd[d])],x2[(sumNd+1):(sumNd+Nd[d])])
      
      Xsd<-Xs[(sumnd+1):(sumnd+nd[d]),]
      mucond.B<-Xsd%*%matrix(betaest,nr=p,nc=1)+upred[d,1]
      sigma2cond.B <- sigma2cond.S
      alphad.B<-(log(z)-mucond.B)/sigma2cond.B   # Now alphad.S is a vector of size nd[d]
      propfin.B[d]<-mean(pnorm(alphad.B))  # Now is the sample mean (equal to the weighted sample mean with weights Nd[d]/nd[d])
      gapfin.B[d]<-mean( pnorm(alphad.B)*(1-(exp(mucond.B+sigma2cond.B/2)*pnorm(alphad.B-sqrt(sigma2cond.B))/pnorm(alphad.B) )/z)  )
      
      sumNd<-sumNd+Nd[d]
    }
    
    #*** Until here the alternative code using analytical formulas
    
    MSEpropsum.B<-MSEpropsum.B+(propfin.B-prop.B)^2
    MSEgapsum.B<-MSEgapsum.B+(gapfin.B-gap.B)^2
    
  }  # End of the bootstrap cycle
  
  # I save predictor values for each Monte Carlo simulation
  # in the columns of a matrix
  
  MSEpropMC.B[,i]<-MSEpropsum.B/B
  MSEgapMC.B[,i]<-MSEgapsum.B/B
  

  #New Fay-Harrio Model
  
  
  New_EB_inc <- eblupFH(formula = propfin.S ~ propfin  , vardir = MSEpropMC.B[,i], method="ML")
  New_EBmse_inc <- mseFH(formula = propfin.S ~ propfin  , vardir = MSEpropMC.B[,i], method="REML",MAXITER = 100, PRECISION = 0.0001)
  
  New_EB_gap <- eblupFH(formula = gapfin.S ~ gapfin , vardir = MSEgapMC.B[,i], method="ML")
  New_EBmse_gap <- mseFH(formula = gapfin.S ~ gapfin  , vardir =  MSEgapMC.B[,i], method="REML")
  
  
  # I save predictor values for each Monte Carlo simulation
  # in the columns of a matrix
  
  propfinMC[,i]<-propfin
  gapfinMC[,i]<-gapfin
  propfinMC.ELL[,i]<-propfin.ELL 
  gapfinMC.ELL[,i]<-gapfin.ELL 
  propfinMC.S[,i] <- propfin.S
  gapfinMC.S[,i] <- gapfin.S
  propfinMC.uc[,i] <- propfin.uc
  gapfinMC.uc[,i]<- gapfin.uc
  propfinMC.New[,i] <- New_EB_inc$eblup
  gapfinMC.New[,i] <- New_EB_gap$eblup
  
}   # End of Monte Carlo simulations


# Calculating the Mean of Bootstrap MSE

MSEpropMC.Br<-rep(0,D)
MSEgapMC.Br<-rep(0,D)

for (d in 1:D){
  MSEpropMC.Br[d]<-mean(MSEpropMC.B[d,])
  MSEgapMC.Br[d]<-mean(MSEgapMC.B[d,])
}

EmpiricalMSEpropPre<-rep(0,D)
EmpiricalMSEgapPre<-rep(0,D)
EmpiricalMSEpropPre.ELL<-rep(0,D)
EmpiricalMSEgapPre.ELL<-rep(0,D)
EmpiricalMSEpropPre.New <- rep(0,D)
EmpiricalMSEgapPre.New <- rep(0,D)
MSEpropMCr.ELL<-rep(0,D)
MSEgapMCr.ELL<-rep(0,D)

for (d in 1:D){
  EmpiricalMSEpropPre[d]<-mean((propfinMC[d,]-propMC[d,])^2)
  EmpiricalMSEgapPre[d]<-mean((gapfinMC[d,]-gapMC[d,])^2)
  EmpiricalMSEpropPre.ELL[d]<-mean((propfinMC.ELL[d,]-propMC[d,])^2)
  EmpiricalMSEgapPre.ELL[d]<-mean((gapfinMC.ELL[d,]-gapMC[d,])^2)
  EmpiricalMSEpropPre.New[d] <- mean((propfinMC.New[d,]-propMC[d,])^2)
  EmpiricalMSEgapPre.New[d] <- mean((gapfinMC.New[d,]-gapMC[d,])^2)
  #MSEpropMCr.ELL[d]<-mean(MSEpropMC.ELL[d,])
  #MSEgapMCr.ELL[d]<-mean(MSEgapMC.ELL[d,])
}

EmpiricalMSEPre<-data.frame(EmpiricalMSEpropPre,EmpiricalMSEgapPre,EmpiricalMSEpropPre.ELL,EmpiricalMSEgapPre.ELL,EmpiricalMSEpropPre.New,EmpiricalMSEgapPre.New)
write.table(EmpiricalMSEPre,"EmpiricalD80MC10000_Modified.txt",row.names=F)


# Plots
differenceMSE <- EmpiricalMSEgapPre*10000-EmpiricalMSEgapPre.New*10000
table <- data.frame(EmpiricalMSEgapPre*10000,EmpiricalMSEgapPre.New*10000,differenceMSE)
table

# Difference between MSEs (CEB and NEW EB)

plot(differenceMSE)

a<-min(EmpiricalMSEgapPre*10000,EmpiricalMSEgapPre.New*10000)
b<-max(EmpiricalMSEgapPre*10000,EmpiricalMSEgapPre.New*10000)
plot(1:D,EmpiricalMSEgapPre*10000,type="n",xlab="Area",ylab="MSE Poverty incidence (x10000)",ylim=c(0,b))
points(1:D,EmpiricalMSEgapPre*10000,type="b",col="blue",pch=1,lty=1,lwd=2)
points(1:D,EmpiricalMSEgapPre.New*10000,type="b",col="red",pch=4,lty=4,lwd=2)
legend("topright",legend=c("CEB","New"),pch=c(1,4),col=c("blue","red"),cex=0.5)
#savePlot(file="ELLVarpovgapD80B300MC100",type="pdf")
#savePlot(file="ELLVarpovgapD80B300MC100",type="eps")
#savePlot(file="ELLVarpovgapD80B300MC100",type="jpg")

############################################################
############################################################

# Calculate EBLUP under FH: incidence and gap

############################################################
############################################################

# Alpha = 0
povertyincidence <- function(y) {
  result <- mean(y < z) 
  return (result)}

# Alpha = 1
povertygap <- function(y) {
  result <- mean((y<z) * (z-y) / z)
  return (result)}

# Direct Estimators

poor <- numeric(n)
poor [Es < z] <- 1

PropDir <- direct(y=poor,dom=area,sweight=weight,domsize=Nd)

# For povery incidence

