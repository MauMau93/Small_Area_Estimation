##################################################################################

#Integrating Census and Survey Information for Small Area Estimation

##################################################################################

rm(list=ls())
#setwd("C:/Isa/Docencia/TFM_MauricioMarcos_2020-21/Simulations")

# Libraries

library(nlme)
library(sae)

# Setting seed
set.seed(826)

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
  #p2[d]<-0.2*d/D
  p2[d]<- 0.2
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
areaP<-NULL #*** New 26/08/2021
for(d in 1:D){
  area<-c(area,rep(d,nd[d]))
  areaP<-c(areaP,rep(d,Nd[d]))
}

weights <- Nd/nd
out_of_samp <- (Nd/nd) -1

# Augmented Auxiliary Variables (aav)
domain <- rep(c(1:D), each=Nd[1]) #***26/08/2021 Es lo mismo que areaP, porque en estas simulaciones, todos los elementos de Nd son iguales
x1_aug <- rep(x1s, each= weights) #***26/08/2021 x1s es de tamaño n y weights es de tamaño D. No es correcto, pero sale bien porque todos los elementos de weights son iguales y toma solo el primero.
x2_aug <- rep(x2s, each= weights) #***26/08/2021 x2s es de tamaño n y weights es de tamaño D. No es correcto, pero sale bien porque todos los elementos de weights son iguales.
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

# Design of Matrixes for Direct Estimation
weightsaug <- rep(weights,50)
Nd.df <- matrix(c(1:D,Nd),ncol=2)
#area_size <- matrix(c(area,weightsaug),ncol=2)#*** 26/08/2021 Not needed

# Design matrix for FH Model #*** 26/08/2021 The X values are fixed in the MC simulations, so much faster to manipulate them out of the MC cycle

#Xs.area <- data.frame(area,Xs[,2:3]) #*** 26/08/2021 FH model uses the population means of the X's, not the sample means
#Xs.means <- aggregate(Xs.area[,2:3], list(Xs.area$area), mean)
#Xs.means <- as.matrix(Xs.means)

x1.mean<-as.vector(tapply(x1,areaP,mean))
x2.mean<-as.vector(tapply(x2,areaP,mean))

#*** For unit-context CEB (UC-CEB)

X1s.mean<-NULL
X2s.mean<-NULL
for(d in 1:D){
  X1s.mean<-c(X1s.mean,rep(x1.mean[d],nd[d]))
  X2s.mean<-c(X2s.mean,rep(x2.mean[d],nd[d]))
}

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

# Census Parameter values

beta<-c(3,0.03,-0.04)
sigmae2<-(0.5)^2
sigmau2<-(0.15)^2
sigmau2/(sigmau2+sigmae2/nd[1])
z<-12    # Poverty line: It is approximately 0.6 Median(E) for the whole population
# gamma=sigmau2/(sigmau2+sigmae2/nd[1])=0.81


# Sample Parameter Values

# C and S produced with same parameters
sigmau2samp <- sigmau2
sigmae2samp <- sigmae2
betasamp <- beta

# C and S produced with different variance
#sigmau2samp <- (0.07)^2
#sigmae2samp <- (0.7)^2
#gamma =0.33
#betasamp <- beta

# C and S produced with different coeficients
#sigmau2samp <- sigmau2
#sigmae2samp <- sigmae2
#gamma= 0.81
#betasamp<- c(2,0.01,-0.02)


# Monte Carlo iterations

#nMC<-500
nMC<-30

# Initialization of matrices that will contain poverty proportions and gaps
# Each column contains the D area quantities for each Monte Carlo simulation

propMC<-matrix(0,nr=D,nc=nMC) #Real Population Poverty Proportion
gapMC<-matrix(0,nr=D,nc=nMC) #Real Population Poverty Gap
propMC.samp <- matrix(0,nr=D,nc=nMC) #Real Sample Poverty Proportion
gapMC.samp <- matrix(0,nr=D,nc=nMC) #Real Sample Poverty Gap

propfinMC<-matrix(0,nr=D,nc=nMC) #Census EB for Poverty Proportion
gapfinMC<-matrix(0,nr=D,nc=nMC)  #Census EB for Poverty Gap
propfinMC.ELL<-matrix(0,nr=D,nc=nMC) #ELL for Poverty Proportion
gapfinMC.ELL<-matrix(0,nr=D,nc=nMC) #ELL for Poverty Gap
propfinMC.S<-matrix(0,nr=D,nc=nMC) # Survey EB for Proportion
gapfinMC.S<-matrix(0,nr=D,nc=nMC) # Survey EB for Gap
#propfinMC.uc <- matrix(0,nr=D,nc=nMC)  #*** 16/08/2021: Not used
#gapfinMC.uc <- matrix(0,nr=D,nc=nMC)  #*** 16/08/2021: Not used
propfinMC.Dir <-matrix(0,nr=D,nc=nMC) # Direct Estimator for Proportion
gapfinMC.Dir <-matrix(0,nr=D,nc=nMC) # Direct Estimator for Gap
propfinMC.FH <-matrix(0,nr=D,nc=nMC) # FH EBLUP for Proportion
gapfinMC.FH <-matrix(0,nr=D,nc=nMC) # FH EBLUP for Gap
propfinMC.New <- matrix(0,nr=D,nc=nMC) # New EB for Proportion
gapfinMC.New <- matrix(0,nr=D,nc=nMC) # New EB for Gap
#propsMC<- matrix(0,nr=D,nc=nMC)
#gapsMC<- matrix(0,nr=D,nc=nMC)


MSEpropMC.B<-matrix(0,nr=D,nc=nMC)
MSEgapMC.B<-matrix(0,nr=D,nc=nMC)
MSEpropMC.ELL<-matrix(0,nr=D,nc=nMC)
MSEgapMC.ELL<-matrix(0,nr=D,nc=nMC)
MSEpropMC.S<-matrix(0,nr=D,nc=nMC)
MSEgapMC.S<-matrix(0,nr=D,nc=nMC)
MSEprop.Dir <-matrix(0,nr=D,nc=nMC)
MSEgap.Dir<-matrix(0,nr=D,nc=nMC)
MSEpropMC.FH<-matrix(0,nr=D,nc=nMC)
MSEgapMC.FH<-matrix(0,nr=D,nc=nMC)
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
  y<-c()
  ys<-c()
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
    usd<-rnorm(1,0,sqrt(sigmau2samp))
    musd<-Xsd%*%betasamp
    esd<-rnorm(nd[d],0,sqrt(sigmae2samp))
    ysd<-musd+usd+esd
    ys<-c(ys,ysd)
    
    sumNd<-sumNd+Nd[d]
    sumnd<-sumnd+nd[d]
    
  }
  
  E<-exp(y)       # People earnings
  
  # Sample data vectors and matrices
  
  #ys<-y[samp]
  #Es<-E[samp]
  Es<-exp(ys)
  
  sumnd<-0
  
  for (d in 1:D){
    
    Esd<-Es[(sumnd+1):(sumnd+nd[d])]
    
    # I save them in the columns of a matrix
    propMC.samp[d,i]<-mean(Esd<z)
    gapMC.samp[d,i]<-mean((Esd<z)*(z-Esd)/z)
    
    sumnd<-sumnd+nd[d]
    
  }
  
  
  #prop.samp[d]<-mean(Es<z)
  #gap.samp[d]<-mean((Es<z)*(z-Es)/z)
  
  # I save them in the columns of a matrix
  
  #propMC.samp[d,i]<-prop.samp[d]
  #gapMC.samp[d,i]<-gapMC.samp[d]
  
  # Fitting of nested-error model to sample data
  
  fit<-lme(ys~x1s+x2s,random=~1|area,method="REML")
  betaest<-fixed.effects(fit)
  upred<-random.effects(fit)
  sigmae2est<-fit$sigma^2 # Varianza residual
  sigmau2est<-as.numeric(VarCorr(fit)[1,1]) # Matriz de covarianzas de las componentes aleatorias del modelo
  Xs<-model.matrix(fit)
  
  ##################################################################
  # Calculating Direct, Eblup FH, Census EB, Survey EB, ELL
  ##################################################################
  

  #Direct Estimators
  
  poor <- numeric(n)
  poor [Es < z] <- 1
  gap <- ((Es<z) * (z-Es) / z)
  
  PropDir <- direct(y=poor,dom=area,sweight=weightsaug,domsize=Nd.df)
  GapDir <- direct(y=gap, dom=area,sweight=weightsaug,domsize = Nd.df)
  
  propfin.Dir <- PropDir$Direct
  propfinVar.Dir <- (PropDir$SD)^2
  gapfin.Dir <- GapDir$Direct
  gapfinVar.Dir <- (GapDir$SD)^2
  
  # Calculation of EBLUP under FH
  #eblupFH(formula, vardir, method = "REML", MAXITER = 100, PRECISION = 0.0001,B = 0, data)
  EblupFHprop <- eblupFH(propfin.Dir ~ x1.mean+x2.mean,vardir=propfinVar.Dir)
  propfin.FH <- EblupFHprop$eblup
  EblupFHgap <- eblupFH(gapfin.Dir ~  x1.mean+x2.mean,vardir=gapfinVar.Dir)
  gapfin.FH <- EblupFHgap$eblup
  
  #*** 26/08/2021 Hasta aquí el código estaba dentro del bucle para cada área
  
  #propfin.Dir<-rep(0,D)
  #gapfin.Dir<-rep(0,D)
  #propfinVar.Dir<-rep(0,D)
  #gapfinVar.Dir<-rep(0,D)
  #propfin.FH<-rep(0,D)
  #gapfin.FH<-rep(0,D)
  propfin<-rep(0,D)
  gapfin<-rep(0,D)
  propfin.ELL<-rep(0,D)
  gapfin.ELL<-rep(0,D)
  propfin.S<-rep(0,D)
  gapfin.S<-rep(0,D)
  #propfin.uc<-rep(0,D) #*** 16/08/2021: Not used
  #gapfin.uc<-rep(0,D) #*** 16/08/2021: Not used
  propfin.New <- rep(0,D)
  gapfin.New <- rep(0,D)
  
  
  sumNd<-0
  sumnd<-0
  
  for (d in 1:D){
    Xd<-cbind(rep(1,Nd[d]),x1[(sumNd+1):(sumNd+Nd[d])],x2[(sumNd+1):(sumNd+Nd[d])])
    
    #Direct Estimators
    
    # Functions for poverty proportion and gap
    
    #poor <- numeric(n) 
    #poor [Es < z] <- 1
    #incidence <- ((Es<z) * (z-Es) / z)
    
    #PropDir <- direct(y=poor,dom=area,sweight=weightsaug,domsize=sampsizes)
    #GapDir <- direct(y=incidence, dom=area,sweight=weightsaug,domsize = sampsizes)
    
    #propfin.Dir[d] <- PropDir$Direct
    #propfinVar.Dir[d] <- (PropDir$SD)^2
    #gapfin.Dir[d] <- GapDir$Direct
    #gapfinVar.Dir[d] <- (GapDir$SD)^2
    
    # Matrix Design for FH Model
    
    #Xs.area <- data.frame(area,Xs[,2:3])
    #Xs.means <- aggregate(Xs.area[,2:3], list(Xs.area$area), mean)
    #Xs.means <- as.matrix(Xs.means)
    
    
    # Calculation of EBLUP under FH
    
    #EblupFHprop <- eblupFH(propfin.Dir ~ Xs.means,vardir=PropDir$SD^2)
    #propfin.FH <- EblupFHprop$eblup
    #EblupFHgap <- eblupFH(gapfin.Dir ~ Xs.means,vardir=GapDir$SD^2)
    #gapfin.FH <- EblupFHgap$eblup
    #*** 26/08/2921 Hasta aquí el código que no debería estar dentro del bucle para cada área
    
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
  
  B<-400
  #B<-200
  
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
    
    #sumNd<-0   #*** 16/08/2021: Not used for Survey EB
    sumnd<-0
    
    for (d in 1:D){
      #Xd<-cbind(rep(1,Nd[d]),x1[(sumNd+1):(sumNd+Nd[d])],x2[(sumNd+1):(sumNd+Nd[d])])
      
      Xsd<-Xs[(sumnd+1):(sumnd+nd[d]),]
      mucond.B<-Xsd%*%matrix(betaest.B,nr=p,nc=1)+upred.B[d,1]
      gammad.B<-sigmau2est.B/(sigmau2est.B+sigmae2est.B/nd[d])
      sigma2cond.B<-sigmau2est.B*(1-gammad.B)+sigmae2est.B
      alphad.B<-(log(z)-mucond.B)/sigma2cond.B   # Now alphad.B is a vector of size nd[d]
      propfin.B[d]<-mean(pnorm(alphad.B))  # Now is the sample mean (equal to the weighted sample mean with weights Nd[d]/nd[d])
      gapfin.B[d]<-mean( pnorm(alphad.B)*(1-(exp(mucond.B+sigma2cond.B/2)*pnorm(alphad.B-sqrt(sigma2cond.B))/pnorm(alphad.B) )/z)  )

      #sumNd<-sumNd+Nd[d]  #*** 16/08/2021: Not used for Survey EB
      sumnd<-sumnd+nd[d]  #*** 16/08/2021: New
    }
    
    #*** Until here the alternative code using analytical formulas
    
    MSEpropsum.B<-MSEpropsum.B+(propfin.B-prop.B)^2
    MSEgapsum.B<-MSEgapsum.B+(gapfin.B-gap.B)^2
    
  }  # End of the bootstrap cycle
  
  # I save predictor values for each Monte Carlo simulation
  # in the columns of a matrix
  
  MSEpropMC.B[,i]<-MSEpropsum.B/B
  MSEgapMC.B[,i]<-MSEgapsum.B/B
  
  
  #New Fay-Harriot Model
  
  
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
  #propfinMC.uc[,i] <- propfin.uc #*** 16/08/2021: Not used
  #gapfinMC.uc[,i]<- gapfin.uc #*** 16/08/2021: Not used
  propfinMC.Dir[,i] <- propfin.Dir
  gapfinMC.Dir[,i] <- gapfin.Dir
  propfinMC.FH[,i] <- propfin.FH
  gapfinMC.FH[,i] <- gapfin.FH
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
EmpiricalMSEpropPre.Dir<-rep(0,D)
EmpiricalMSEgapPre.Dir<-rep(0,D)
EmpiricalMSEpropPre.FH<-rep(0,D)
EmpiricalMSEgapPre.FH<-rep(0,D)
EmpiricalMSEpropPre.New <- rep(0,D)
EmpiricalMSEgapPre.New <- rep(0,D)
#EmpiricalMSEpropPre.uc <- rep(0,D)
#EmpiricalMSEgapPre.uc <- rep(0,D)
MSEpropMCr.ELL<-rep(0,D)
MSEgapMCr.ELL<-rep(0,D)

#x1.mean<-as.vector(tapply(x1,areaP,mean))

# Value of estimators for each area
propRate <- apply(propMC,1,mean)
#gapRate <- apply(gapMC,1,mean)
#sampPropRate <- apply(propMC.samp,1,mean)
#gapPropRate <- apply(gapMC.samp,1,mean)
#propPred.ELL
#gapPred.ELL
#propPred.Dir
#gapPred.Dir
#propPred.FH
#gapPred.FH
#propPred.New
#gapPred.New
#propPred.SEB
#gapPred.SEB
#propPred.CEB
#gapPred.CEB


for (d in 1:D){
  EmpiricalMSEpropPre[d]<-mean((propfinMC[d,]-propMC[d,])^2)
  EmpiricalMSEgapPre[d]<-mean((gapfinMC[d,]-gapMC[d,])^2)
  EmpiricalMSEpropPre.ELL[d]<-mean((propfinMC.ELL[d,]-propMC[d,])^2)
  EmpiricalMSEgapPre.ELL[d]<-mean((gapfinMC.ELL[d,]-gapMC[d,])^2)
  EmpiricalMSEpropPre.Dir[d] <- mean((propfinMC.Dir[d,]-propMC[d,])^2)
  EmpiricalMSEgapPre.Dir[d] <- mean((gapfinMC.Dir[d,]-gapMC[d,])^2)
  EmpiricalMSEpropPre.FH[d] <- mean((propfinMC.FH[d,]-propMC[d,])^2)
  EmpiricalMSEgapPre.FH[d] <- mean((gapfinMC.FH[d,]-gapMC[d,])^2)
  EmpiricalMSEpropPre.New[d] <- mean((propfinMC.New[d,]-propMC[d,])^2)
  EmpiricalMSEgapPre.New[d] <- mean((gapfinMC.New[d,]-gapMC[d,])^2)
  #EmpiricalMSEpropPre.uc[d] <- mean((propfinMC.uc[d,]-propMC[d,])^2)
  #EmpiricalMSEgapPre.uc[d] <- mean((gapfinMC.uc[d,]-gapMC[d,])^2)
  #MSEpropMCr.ELL[d]<-mean(MSEpropMC.ELL[d,])
  #MSEgapMCr.ELL[d]<-mean(MSEgapMC.ELL[d,])
}

EmpiricalMSEPre<-data.frame(EmpiricalMSEpropPre,EmpiricalMSEgapPre,EmpiricalMSEpropPre.ELL,EmpiricalMSEgapPre.ELL,EmpiricalMSEpropPre.New,EmpiricalMSEgapPre.New,EmpiricalMSEpropPre.Dir,EmpiricalMSEgapPre.Dir,EmpiricalMSEpropPre.FH,EmpiricalMSEgapPre.FH)
write.table(EmpiricalMSEPre,"EmpiricalD80MC10000_Modified.txt",row.names=F)
EmpiricalMSEPre_prop<- data.frame(EmpiricalMSEpropPre,EmpiricalMSEpropPre.ELL,EmpiricalMSEpropPre.New,EmpiricalMSEpropPre.Dir,EmpiricalMSEpropPre.FH)
EmpiricalMSEPre_gap<- data.frame(EmpiricalMSEgapPre,EmpiricalMSEgapPre.ELL,EmpiricalMSEgapPre.New,EmpiricalMSEgapPre.Dir,EmpiricalMSEgapPre.FH)


# Plot
differenceMSE <- (EmpiricalMSEgapPre)-(EmpiricalMSEgapPre.New)
print(differenceMSE)
dataframe <- data.frame(EmpiricalMSEgapPre*10000,EmpiricalMSEgapPre.New*10000)

#table <- data.frame(EmpiricalMSEgapPre*10000,EmpiricalMSEgapPre.New*10000,differenceMSE)
#table

# Difference between MSEs (CEB and NEW EB)

#plot(differenceMSE)

Dataframe <- data.frame(EmpiricalMSEpropPre,EmpiricalMSEpropPre.Dir)

a<-min(EmpiricalMSEgapPre.ELL*10000,EmpiricalMSEgapPre*10000)
b<-max(EmpiricalMSEgapPre.ELL*10000,EmpiricalMSEgapPre*10000)
plot(1:D,EmpiricalMSEgapPre*10000,type="n",xlab="Area",ylab="MSE Poverty incidence (x10000)",ylim=c(0,b))
points(1:D,EmpiricalMSEgapPre.ELL*10000,type="b",col="blue",pch=1,lty=1,lwd=2)
points(1:D,EmpiricalMSEgapPre.New*10000,type="b",col="red",pch=4,lty=4,lwd=2)
points(1:D,EmpiricalMSEgapPre*10000,type="b",col="yellow",pch=2,lty=2,lwd=2)
legend("top",legend=c("ELL","CEB","NEW"),pch=c(1,4,2),col=c("blue","red","yellow"),cex=0.5)

#savePlot(file="ELLVarpovgapD80B300MC100",type="pdf")
#savePlot(file="ELLVarpovgapD80B300MC100",type="eps")
#savePlot(file="ELLVarpovgapD80B300MC100",type="jpg")

############################################################
############################################################

