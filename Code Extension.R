##################################################################################
# Simulation experiment to compare the performance of Census EB and ELL estimators
# of poverty indicators, when the survey is not a subset of the census.
# 
# Author: Isabel Molina Peralta
# Date: 30/08/2020
##################################################################################


#########################################
# Generation of a population
#########################################

rm(list=ls())
#setwd("C:/Isa/AdvisorWB/AdvisorBootstrapPaulACorral/Simulations/NewModelBasedSimulations")

# Setting of sizes and initialization of vectors

library(nlme)
N<-20000        # Number of population units
D<-80           # Number of areas
Nd<-rep(250,D)      # Number of population units in each area
p<-3            # Number of auxiliary variables (intercept+2 dummies)
p1<-rep(0,D)       # Probability of the first category of X1
x1<-NULL       # X1 values
p2<-rep(0,D)       # Probability of the second category of X2
x2<-NULL       # X2 values
#
## Generation of auxiliary variables (2 dummies with probability of 1 given by p1 and p2)
#
for (d in 1:D){
    p1[d]<-0.3+0.5*d/D
    p2[d]<-0.2
    x1<-c(x1,rbinom(Nd[d],size=1,prob=p1[d]))
    x2<-c(x2,rbinom(Nd[d],size=1,prob=p2[d]))
}
#
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

# Sample elements of aux. variables
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

# gamma=sigmau2/(sigmau2+sigmae2/nd[1])=0.71
# We should try with different values for gamma, such as gamma=0.3,0.6

# Monte Carlo iterations

#nMC<-10000
nMC<-100

# Initialization of matrices that will contain poverty proportions and gaps
# Each column is for the D area quantities for each Monte Carlo simulation

propMC<-matrix(0,nr=D,nc=nMC)
gapMC<-matrix(0,nr=D,nc=nMC)
propsMC<-matrix(0,nr=D,nc=nMC)
gapsMC<-matrix(0,nr=D,nc=nMC)
propfinMC<-matrix(0,nr=D,nc=nMC)
gapfinMC<-matrix(0,nr=D,nc=nMC)
propfinMC.ELL<-matrix(0,nr=D,nc=nMC)
gapfinMC.ELL<-matrix(0,nr=D,nc=nMC)
MSEpropMC.ELL<-matrix(0,nr=D,nc=nMC)
MSEgapMC.ELL<-matrix(0,nr=D,nc=nMC)

time1<-Sys.time()
time1

for (i in 1:nMC){  # Start of Monte Carlo cycle: Generation of populations and calculation of EBPs
  
  print(i)
  print(Sys.time())
  
  sumNd<-0
  sumnd<-0
  y<-NULL
  ys<-NULL
  
  # Generation of population values of the response, area by area
  # Calculation of true area poverty proportions and gaps
  
  for (d in 1:D){
    
    ed<-rnorm(Nd[d],0,sqrt(sigmae2))
    ud<-rnorm(1,0,sqrt(sigmau2))
    Xd<-cbind(rep(1,Nd[d]),x1[(sumNd+1):(sumNd+Nd[d])],x2[(sumNd+1):(sumNd+Nd[d])])
    
    mud<-Xd%*%beta
    yd<-mud+ud+ed
    Ed<-exp(yd)
    y<-c(y,yd)
    
    # I save them in the columns of a matrix
    propMC[d,i]<-mean(Ed<z)
    gapMC[d,i]<-mean((Ed<z)*(z-Ed)/z)
    
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
  
  # National poverty line
  # 0.6*median(E) # Poverty line
  # Resultado: 20.256
  
  # 100*mean(E<z)     # Population poverty incidence: Debería estar alrededor de un 15%
  # Resultado: 15.585
  # 100*mean((E<z)*(z-E)/z)   # Population poverty gap
  # Resultado: 3.437
  # data.frame(prop*100,gap*100)
  
  # Sample data 
  
  Es<-exp(ys)
  
  # Non-Sample data vectors and matrices
  
  yr<-y[-samp]
  
  # Sample poverty proportions and gaps at national level
  #100*mean(Es<z)
  #100*mean((Es<z)*(z-Es[1:n])/z)
  
  # Sample poverty proportions and gaps for the areas
  
  sumnd<-0
  
  for (d in 1:D){
    
    Esd<-Es[(sumnd+1):(sumnd+nd[d])]
    
    # I save them in the columns of a matrix
    propsMC[d,i]<-mean(Esd<z)
    gapsMC[d,i]<-mean((Esd<z)*(z-Esd)/z)
    
    sumnd<-sumnd+nd[d]
    
  }
  
  #data.frame(prop*100,props*100,gap*100,gaps*100)
  
  # Fitting of nested-error model to sample data
  
  fit<-lme(ys~x1s+x2s,random=~1|area,method="REML")
  Xs<-model.matrix(fit)
  betaest<-fixed.effects(fit)
  upred<-random.effects(fit)
  sigmae2est<-fit$sigma^2    # Varianza residual
  sigmau2est<-as.numeric(VarCorr(fit)[1,1]) # Matriz de covarianzas de las componentes aleatorias del modelo
  
  
  ############################################################
  # Generation of census values and calculation of Census 
  # EB and ELL poverty proportions and gaps for the areas
  ############################################################
  
  L<-50
  propsum<-rep(0,D)
  gapsum<-rep(0,D)
  propsum.ELL<-rep(0,D)
  gapsum.ELL<-rep(0,D)
  sumSqprop<-rep(0,D)
  sumSqgap<-rep(0,D)
  
  for (ell in 1:L){   ### Start of generations
    
    sumNd<-0
    
    for (d in 1:D){
      
      # Generation of census values from the distrib. given the sample data (Census EB prediction)
      
      Xd<-cbind(rep(1,Nd[d]),x1[(sumNd+1):(sumNd+Nd[d])],x2[(sumNd+1):(sumNd+Nd[d])])
      mudpred<-Xd%*%matrix(betaest,nr=p,nc=1)+upred[d,1]
      
      gammad<-sigmau2est/(sigmau2est+sigmae2est/nd[d])
      sigmav2<-sigmau2est*(1-gammad)
      vd<-rnorm(1,0,sqrt(sigmav2))
      ed<-rnorm(Nd[d],0,sqrt(sigmae2est))
      
      ydnew<-mudpred+vd+ed
      Ednew<-exp(ydnew)
      
      # EB predictors of poverty measures for each generation
      
      propsum[d]<-propsum[d]+mean(Ednew<z)
      gapsum[d]<-gapsum[d]+mean((Ednew<z)*(z-Ednew)/z)
      
      # Generation of population values of y from the fitted model
      # (ELL method with clusters = areas), using parametric bootstrap
      # and without heteroscedasticity
      
      ed.ELL<-ed
      ud.ELL<-rnorm(1,0,sqrt(sigmau2est))
      
      mud.ELL<-Xd%*%matrix(betaest,nr=p,nc=1)
      yd.ELL<-mud.ELL+ud.ELL+ed.ELL
      Ed.ELL<-exp(yd.ELL)
      
      # ELL predictors of poverty measures for each generation
      
      propsum.ELL[d]<-propsum.ELL[d]+mean(Ed.ELL<z)
      gapsum.ELL[d]<-gapsum.ELL[d]+mean((Ed.ELL<z)*(z-Ed.ELL)/z)
      
      # For the calculation of ELL estimator of MSE
      
      sumSqprop[d]<-sumSqprop[d]+(mean(Ed.ELL<z))^2
      sumSqgap[d]<-sumSqgap[d]+(mean((Ed.ELL<z)*(z-Ed.ELL)/z))^2
      
      sumNd<-sumNd+Nd[d]
      
    }
    
  }  # End of generations
  
  # EB predictors of poverty measures
  
  propfin<-propsum/L
  gapfin<-gapsum/L
  
  # ELL predictors of poverty measures and ELL mean squared error
  
  propfin.ELL<-propsum.ELL/L
  gapfin.ELL<-gapsum.ELL/L
  
  MSEpropMC.ELL[,i]<-sumSqprop/L-propfin.ELL^2
  MSEgapMC.ELL[,i]<-sumSqgap/L-gapfin.ELL^2
  
  data.frame(prop,propfin,propfin.ELL,gap,gapfin,gapfin.ELL)
  
  # I save predictor values for each Monte Carlo simulation
  # in the columns of a matrix
  
  propfinMC[,i]<-propfin
  gapfinMC[,i]<-gapfin
  propfinMC.ELL[,i]<-propfin.ELL
  gapfinMC.ELL[,i]<-gapfin.ELL
  
}   # End of Monte Carlo simulations

time2 <- Sys.time()
print(difftime(time2,time1,units="mins"))


#####################################################################

EmpiricalMSEpropPre<-rep(0,D)
EmpiricalMSEgapPre<-rep(0,D)
EmpiricalMSEpropPre.ELL<-rep(0,D)
EmpiricalMSEgapPre.ELL<-rep(0,D)
MSEpropMCr.ELL<-rep(0,D)
MSEgapMCr.ELL<-rep(0,D)

for (d in 1:D){
  EmpiricalMSEpropPre[d]<-mean((propfinMC[d,]-propMC[d,])^2)
  EmpiricalMSEgapPre[d]<-mean((gapfinMC[d,]-gapMC[d,])^2)
  EmpiricalMSEpropPre.ELL[d]<-mean((propfinMC.ELL[d,]-propMC[d,])^2)
  EmpiricalMSEgapPre.ELL[d]<-mean((gapfinMC.ELL[d,]-gapMC[d,])^2)
  MSEpropMCr.ELL[d]<-mean(MSEpropMC.ELL[d,])
  MSEgapMCr.ELL[d]<-mean(MSEgapMC.ELL[d,])
}

EmpiricalMSEPre<-data.frame(EmpiricalMSEpropPre,EmpiricalMSEgapPre,EmpiricalMSEpropPre.ELL,EmpiricalMSEgapPre.ELL,MSEpropMCr.ELL,MSEgapMCr.ELL)
write.table(EmpiricalMSEPre,"EmpiricalD80MC10000_Modified.txt",row.names=F)

#############################################################################################
#From now on, my own code.
# First, I need to replicate the survey artificially to create a Census. I consider that all of the observations have the same weight in the sample: ni/Ni

weights <- Nd/nd
out_of_samp <- (Nd/nd) -1
# Augmented Auxiliary Variables (aav)
domain <- rep(c(1:80), each=250)
x1_aug <- rep(x1s, each= weights)
x2_aug <- rep(x2s, each= weights)
X_aug<-cbind(domain,x1_aug,x2_aug)

write.table(data.frame(X_aug),"AuxVarAug.txt")
X_aug<-read.table("AuxVarAug.txt",header=TRUE)
X_aug$domain <- as.factor(X_aug$domain)

# Out-of sample Augmented Auxiliary Variables 
domain_out <- rep(c(1:80), each=200)
x1_o_aug <- rep(x1s, each= out_of_samp)
x2_o_aug <- rep(x2s, each= out_of_samp)
X_o_aug<-cbind(domain_out,x1_o_aug,x2_o_aug)

write.table(data.frame(X_o_aug),"AuxVarAug_out.txt")
X_o_aug<-read.table("AuxVarAug_out.txt",header=TRUE)
X_o_aug$domain_out <- as.factor(X_o_aug$domain_out)

# Calculating  Artifitial Census EB with sae package
package(sae)
area_codes <- unique(area)
area_codes <- as.factor(area_codes)
povertyincidence <- function(y) {
  result <- mean(y < 20.256) 
  return (result)}


ArtCen_EB <- ebBHF(Es ~ x1s + x2s, dom = area,
             selectdom = area_codes, Xnonsample = X_o_aug, MC = 50,
             constant = 0, indicator = povertyincidence)

ArtCen_EB$fit$summary

plot(ArtCen_EB$fit$residuals, xlab = "Index", ylab = "Residuals", cex.axis = 1.5,
    cex.lab = 1.5, ylim = c(-2, 2), col = 4)
abline(h = 0)

hist(ArtCen_EB$fit$residuals, prob = TRUE, xlab = "Residuals", ylab = "", main = "",
       cex.axis = 1.5, cex.lab = 1.5, xlim = c(-2, 2), ylim = c(0, 1))

# Calculating the MSE

pbmse.EB <- pbmseebBHF(Es ~ x1s + x2s, dom = area,
                       selectdom = area_codes, Xnonsample = X_o_aug,
                        B = 200, MC = 50, constant = 3500,
                        indicator = povertyincidence)

# Calculating the CVs

pbcv.EB <- 100 * sqrt(pbmse.EB$mse$mse) / abs(pbmse.EB$est$eb$eb)

# Results with this new ArtCen_EB

results.EB <- data.frame(Domain = pbmse.EB$est$eb$domain,
                          SampleSize = pbmse.EB$est$eb$sampsize,
                          EB = pbmse.EB$est$eb$eb, cv.EB = pbcv.EB, MSE.EB = pbmse.EB$mse$mse)


# Building the New Fay-Harriot Model

New_EB <- eblupFH(formula = ArtCen_EB$eb$eb ~ propfin  , vardir = pbmse.EB$mse$mse, method="REML",MAXITER = 100, PRECISION = 0.0001)
New_EBmse <- mseFH(formula = ArtCen_EB$eb$eb ~ propfin  , vardir = pbmse.EB$mse$mse, method="REML",MAXITER = 100, PRECISION = 0.0001)
cv.FH <- 100 * sqrt(New_EBmse$mse) / New_EBmse$est$eblup

results <- data.frame(Area = area_codes, DIR = ArtCen_EB$eb$eb ,
                       cv.DIR = pbcv.EB, eblup.FH = New_EBmse$est$eblup, cv.FH)

# Ploting

# Plot of the estimations: direct, new eblupFH
 plot(results$DIR, type = "n", ylab = "Estimate", ylim = c(-0.4, 1.6),
          xlab = "area (sorted by decreasing sample size)", cex.axis = 1.5,
          cex.lab = 1.5)
 points(results$DIR, type = "b", col = 1, lwd = 2, pch = 1, lty = 1)
 points(results$eblup.FH, type = "b", col = 4, lwd = 2, pch = 4, lty = 2)
 legend("top", legend = c("Direct", "EBLUP FH"), ncol = 2, col = c(1, 4), lwd = 2,
          pch = c(1, 4), lty = c(1, 2), cex = 1.3)
 plot(results$cv.DIR, type = "n", ylab = "CV", ylim = c(1, 40),
        xlab = "area (sorted by decreasing sample size)", cex.axis = 1.5,
        cex.lab = 1.5)
 points(results$cv.DIR, type = "b", col = 1, lwd = 2, pch = 1, lty = 1)
 points(results$cv.FH, type = "b", col = 4, lwd = 2, pch = 4, lty = 2)
 legend("top", legend = c("Direct", "EBLUP FH"), ncol = 2, col = c(1, 4), lwd = 2,
          pch = c(1, 4), lty = c(1, 2), cex = 1.3)
 
 # Plot for the mse: Census EB, new eblupFH
 New_EB_mse_augmented <- New_EBmse$mse*10000
 CensusEB_mse_augmented<- EmpiricalMSEpropPre*10000
 plot(New_EBmse$mse, type = "n", ylab = "CV", ylim = c(0, 30),
      xlab = "area", cex.axis = 1.5,
      cex.lab = 1.5)
 #points(ArtCen_EB$eb$eb, type = "b", col = 1, lwd = 2, pch = 1, lty = 1)
 points(New_EB_mse_augmented, type = "b", col = 4, lwd = 2, pch = 4, lty = 2)
 points(CensusEB_mse_augmented, type = "b", col = 2, lwd = 2, pch = 2, lty = 3)
 legend("top", legend = c( "EBLUPFH_mse","Census EB"), ncol = 2, col = c(1, 4), lwd = 2,
        pch = c(1, 4), lty = c(1, 2), cex = 1.3)

