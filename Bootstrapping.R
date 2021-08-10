# Bootrstraping

rm(list=ls())
#setwd("C:/Isa/AdvisorWB/AdvisorBootstrapPaulACorral/Simulations/NewBootstrap_Modified")

library(nlme)

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
#
 Generation of auxiliary variables (2 dummies with probability of 1 given by p1 and p2)

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
}


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
nMC<-20

# Initialization of matrices that will contain poverty proportions and gaps
# Each column contains the D area quantities for each Monte Carlo simulation

propMC<-matrix(0,nr=D,nc=nMC)
gapMC<-matrix(0,nr=D,nc=nMC)

MSEpropMC.B<-matrix(0,nr=D,nc=nMC)
MSEgapMC.B<-matrix(0,nr=D,nc=nMC)

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
  
  # Fitting of nested-error model to sample data
  
  fit<-lme(ys~x1s+x2s,random=~1|area,method="REML")
  betaest<-fixed.effects(fit)
  upred<-random.effects(fit)
  sigmae2est<-fit$sigma^2 # Varianza residual
  sigmau2est<-as.numeric(VarCorr(fit)[1,1]) # Matriz de covarianzas de las componentes aleatorias del modelo
  Xs<-model.matrix(fit)
  
  ############################################################
  # Parametric Bootstrap
  # for the calculation of the MSE of the EBP
  ############################################################
  
  #B<-500
  B<-50
  
  MSEpropsum.B<-rep(0,D)
  MSEgapsum.B<-rep(0,D)
  
  for (b in 1:B){   ### Start of bootstrap cycle
    
    #print(cbind("Bootstrap iteration",b,"Monte Carlo iteration",i))
    
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
  
    
    # Alternative code using analytical formulas for Census EB predictors of poverty incidence and poverty gap
    
    propfin.B<-rep(0,D)
    gapfin.B<-rep(0,D)
    
    sumNd<-0
    
    for (d in 1:D){
      Xd<-cbind(rep(1,Nd[d]),x1[(sumNd+1):(sumNd+Nd[d])],x2[(sumNd+1):(sumNd+Nd[d])])
      mudcond.B<-Xd%*%matrix(betaest.B,nr=p,nc=1)+upred.B[d,1]
      gammad.B<-sigmau2est.B/(sigmau2est.B+sigmae2est.B/nd[d])
      sigma2cond.B<-sigmau2est.B*(1-gammad.B)+sigmae2est.B
      alphad.B<-(log(z)-mudcond.B)/sigma2cond.B  ## alphad is a vector because Xd is a matrix
      propfin.B[d]<-mean(pnorm(alphad.B))
      gapfin.B[d]<-mean( pnorm(alphad.B)*(1-(exp(mudcond.B+sigma2cond.B/2)*pnorm(alphad.B-sqrt(sigma2cond.B))/pnorm(alphad.B) )/z)  )
      
      sumNd<-sumNd+Nd[d]
    }
    
    # Calculate the sums over each bootstrap to calculate the average between repetitions
    
    
    MSEpropsum.B<-MSEpropsum.B+(propfin.B-prop.B)^2
    MSEgapsum.B<-MSEgapsum.B+(gapfin.B-gap.B)^2
    
  }  # End of the bootstrap cycle
  
  # I save predictor values for each Monte Carlo simulation
  # in the columns of a matrix
  
  MSEpropMC.B[,i]<-MSEpropsum.B/B
  MSEgapMC.B[,i]<-MSEgapsum.B/B
  
}   # End of Monte Carlo simulations

time2 <- Sys.time()
print(difftime(time2,time1,units="mins"))

# Means over Monte Carlo simulations of true values, sample values,
# EB and ELL predictors

MSEpropMC.Br<-rep(0,D)
MSEgapMC.Br<-rep(0,D)

for (d in 1:D){
  MSEpropMC.Br[d]<-mean(MSEpropMC.B[d,])
  MSEgapMC.Br[d]<-mean(MSEgapMC.B[d,])
}

# True MSE (approximated by Monte Carlo)

TrueMSEPre<-read.table("EmpiricalD80MC10000_Modified.txt",header=TRUE)
TrueMSEpropPre<-TrueMSEPre[,1]
TrueMSEgapPre<-TrueMSEPre[,2]

a<-min(TrueMSEpropPre*10000,MSEpropMC.Br*10000)
b<-max(TrueMSEpropPre*10000,MSEpropMC.Br*10000)
plot(1:D,TrueMSEpropPre*10000,type="n",xlab="Area",ylab="MSE Poverty incidence (x10000)",ylim=c(0,b))
points(1:D,TrueMSEpropPre*10000,type="b",col="blue",pch=1,lty=1,lwd=2)
points(1:D,MSEpropMC.Br*10000,type="b",col="red",pch=4,lty=4,lwd=2)
#points(1:D,CMSEpropMC.Br*10000,type="b",col=2,pch=2,lty=2,lwd=2)
legend(0,5,legend=c("True MSE(Census EB)","Bootstrap MSE(Census EB)"),col=c(1,4),pch=c(1,4),lwd=rep(2,2),lty=c(1,4))
#legend(0,5,legend=c("True MSE","Bootstrap MSE","C Bootstrap MSE"),col=c(1,4,2),pch=c(1,4,2),lwd=rep(2,3),lty=c(1,4,2))

#savePlot(file="BootstrapCensusEBpovincD80B50MC20_Modif",type="pdf")
#savePlot(file="BootstrapCensusEBpovincD80B50MC20_Modif",type="eps")
#savePlot(file="BootstrapCensusEBpovincD80B50MC20_Modif",type="jpg")


a<-min(TrueMSEgapPre*10000,MSEgapMC.Br*10000)
b<-max(TrueMSEgapPre*10000,MSEgapMC.Br*10000)
plot(1:D,TrueMSEgapPre*10000,type="n",xlab="Area",ylab="MSE, Poverty gap (x10000)",ylim=c(0,b))
points(1:D,TrueMSEgapPre*10000,type="p",col="blue",pch=1,lty=1,lwd=2)
points(1:D,MSEgapMC.Br*10000,type="p",col="red",pch=4,lty=4,lwd=2)
#points(1:D,MSEgapMC.B2r*10000,type="b",col=2,pch=2,lty=2,lwd=2)
legend(0,0.2,legend=c("True MSE (Census EB)","Bootstrap MSE (Census EB)"),col=c(1,4),pch=c(1,4),lwd=rep(2,2),lty=c(1,4))
#legend(0,2,legend=c("True MSE","Bootstrap MSE","C Bootstrap MSE"),col=c(1,4,2),pch=c(1,4,2),lwd=rep(2,3),lty=c(1,4,2))

#savePlot(file="BootstrapCensusEBpovgapD80B50MC20_Modif",type="pdf")
#savePlot(file="BootstrapCensusEBpovgapD80B50MC20",type="eps")
#savePlot(file="BootstrapCensusEBpovgapD80B50MC20",type="jpg")

MeanRes.MSE<-data.frame(TrueMSEpropPre,MSEpropMC.Br,TrueMSEgapPre,MSEgapMC.Br)
write.table(MeanRes.MSE,"ResultsBootMC100B300D80.txt")

MeanResMSE<-read.table("ResultsBootMC100B300D80.txt",header=T)


TrueMSEpropPre.ELL<-TrueMSEPre[,3]
TrueMSEgapPre.ELL<-TrueMSEPre[,4]

MSEpropMCr.ELL<-TrueMSEPre[,5]
MSEgapMCr.ELL<-TrueMSEPre[,6]

a<-min(TrueMSEpropPre.ELL*10000,MSEpropMCr.ELL*10000)
b<-max(TrueMSEpropPre.ELL*10000,MSEpropMCr.ELL*10000)
plot(1:D,TrueMSEpropPre.ELL*10000,type="n",xlab="Area",ylab="MSE Poverty incidence (x10000)",ylim=c(0,b))
points(1:D,TrueMSEpropPre.ELL*10000,type="b",col=1,pch=1,lty=1,lwd=2)
points(1:D,MSEpropMCr.ELL*10000,type="b",col=4,pch=4,lty=4,lwd=2)
legend(0,10,legend=c("True MSE","ELL Bootstrap Var"),col=c(1,4),pch=c(1,4),lwd=rep(2,2),lty=1:2)
#savePlot(file="ELLVarpovincD80B300MC100",type="pdf")
#savePlot(file="ELLVarpovincD80B300MC100",type="eps")
#savePlot(file="ELLVarpovincD80B300MC100",type="jpg")

a<-min(TrueMSEgapPre.ELL*10000,MSEgapMCr.ELL*10000)
b<-max(TrueMSEgapPre.ELL*10000,MSEgapMCr.ELL*10000)
plot(1:D,TrueMSEgapPre.ELL*10000,type="n",xlab="Area",ylab="MSE Poverty incidence (x10000)",ylim=c(0,b))
points(1:D,TrueMSEgapPre.ELL*10000,type="b",col=1,pch=1,lty=1,lwd=2)
points(1:D,MSEgapMCr.ELL*10000,type="b",col=4,pch=4,lty=4,lwd=2)
legend(0,1,legend=c("True MSE","ELL Bootstrap Var"),col=c(1,4),pch=c(1,4),lwd=rep(2,2),lty=1:2)
#savePlot(file="ELLVarpovgapD80B300MC100",type="pdf")
#savePlot(file="ELLVarpovgapD80B300MC100",type="eps")
#savePlot(file="ELLVarpovgapD80B300MC100",type="jpg")
