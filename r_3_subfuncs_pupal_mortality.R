# This file fits pupal temperature-dependent mortality function
# using laboratory data available in:
#Phelps 1973 The effect of temperature on fat consumption during the puparial stages of
# Glossina morsitans morsitans Westw. (Dipt., Clossinidae) under laboratory conditions,
# and its implication in the field. Bull. Ent. Res. 62, 423-438
# It then plots Figure 2

library("binom")                                            # required pacakges
require(bbmle)

# #**************************Data*************************************
# phelps.temp <- c(16,17,18,20,22,24,25,26,27,28,29,30,31,32)                                 # Vector of temperatures used by Phelps
# 
# phelps.mortality <- c(52.2,21.9,15.3,10.3,4.9,3.6,5,3.5,10.3,3.8,4,7.2,20.5,45.8)/100         # Proportion of pupae that did not survive to emergence sex unknown
# 
# phelps.m.emerge <- c(25,53,26,27,25,27,25,34,24,14,26,30,48,26)                             # Number of males that emerged
# 
# phelps.f.emerge <- c(28,37,24,25,32,26,32,20,28,12,22,22,49,39)                             # Number of females that emerged


phelps.temp <- c(17,18,20,22,24,25,26,27,28,29,30,31,32)                                 # Vector of temperatures used by Phelps

phelps.mortality <- c(21.9,15.3,10.3,4.9,3.6,5,3.5,10.3,3.8,4,7.2,20.5,45.8)/100         # Proportion of pupae that did not survive to emergence sex unknown

phelps.m.emerge <- c(53,26,27,25,27,25,34,24,14,26,30,48,26)                             # Number of males that emerged

phelps.f.emerge <- c(37,24,25,32,26,32,20,28,12,22,22,49,39)                             # Number of females that emerged



phelps.emerged <- phelps.m.emerge + phelps.f.emerge                                      # Total number that emerged

phelps.number.died <- round(phelps.mortality*phelps.emerged / (1-phelps.mortality),0)    # Estimate number that died

phelps.total <- phelps.emerged + phelps.number.died                                      # Total

phelps.b <- binom.confint(phelps.number.died,phelps.total,methods="exact")               # Confidence intervals
#***************************************************************



# #****************Plot the data*******************************
par(mar=c(4,4,2,2),cex=1,mfcol=c(1,1))

plot(phelps.temp
     ,phelps.b$mean
     ,bty="n"
     ,xlab=expression(paste("Temperature (",degree,"C)"))
     ,ylab="Proportion of pupae failing to emerge"
     ,pch=19
     ,col="blue"
     ,ylim=c(0,1)
     ,xlim=c(16,36))
segments(phelps.temp,phelps.b$lower,phelps.temp,phelps.b$upper)
#**************************************************************



#**************Convert probabilities to instantaneous mortality********
inst.mort <- -log(1-phelps.b$mean)
#********************************************************************



# #**************Mortality function*******************************
mort_func <- function(k1,k2,k3,k4,temp){                                 # This function has 4 parameters to be fitted - k1 to k4
  inst.mort <- k1*exp(-(k2*temp))+k3*exp(k4*(temp-25))
  return(inst.mort)
}
#*******************************************************************



#*************Likelihood function**************************
nll <- function(log_k1
                ,log_k2
                ,log_k3
                ,log_k4
                ,dat=cbind.data.frame(inst.mort,phelps.temp)){
   k1 <- exp(log_k1)                                                       # Parameters to be fitted
   k2 <- exp(log_k2)
   k3 <- exp(log_k3)
   k4 <- exp(log_k4)

   mod.mort <- mort_func(k1=k1,k2=k2,k3=k3,k4=k4,temp=dat$phelps.temp)    # Model output
   obs.mort <- dat$inst.mort                                              # Observed data
   ll <- sum(dnorm(mod.mort,mean=obs.mort,log=T))                         # Likelihood assuming normal distribution
   return(-ll)
}
#***********************************************************



#************initial parameters***************************
init.pars <- c(log_k1=log(24)
               ,log_k2=log(0.27)
               ,log_k3=log(0.0002)
               ,log_k4=log(1))
#************************************************************



#***********Fitting******************
 objFXN <- function(fit.params
                    ,dat=cbind.data.frame(inst.mort,phelps.temp)) {
   parms <- fit.params
   nll(fit.params[1],fit.params[2],fit.params[3],fit.params[4], dat = dat) ## then call likelihood
 }
 
 trace <- 3
 
 optim.vals <- optim(par = init.pars
                     , objFXN
                     , dat = cbind.data.frame(inst.mort,phelps.temp)
                     , control = list(trace = trace, maxit = 200)
                     , method = "SANN")
 exp(optim.vals$par) #
 
 
 optim.vals <- optim(par = optim.vals$par
                     , objFXN
                     , dat = cbind.data.frame(inst.mort,phelps.temp)
                     , control = list(trace = trace, maxit = 200)
                     , method = "SANN")
 exp(optim.vals$par) #

 
 optim.vals <- optim(par = optim.vals$par
                     , objFXN
                     , dat = cbind.data.frame(inst.mort,phelps.temp)
                     , control = list(trace = trace, maxit = 1000, reltol = 10^-7)
                     , method = "Nelder-Mead" #
                     , hessian = T)
 optim.vals # convergence 0 means algorithm converged

 MLEfits <- optim.vals$par
 exp(MLEfits)
#**********************************************************
 
 
 
#********************Confidence intervals*********************************** 
fisherInfMatrix <- solve(optim.vals$hessian) ## invert the Hessian, to estimate the covar-var matrix of parameter estimates
fisherInfMatrix

# Finds the critical z value
conf.level <- 0.95
crit <- qnorm((1 + conf.level)/2)

ci <- optim.vals$par[1] + c(-1, 1) * crit * sqrt(abs(fisherInfMatrix[1, 1]))
exp(ci)

ci <- optim.vals$par[2] + c(-1, 1) * crit * sqrt(abs(fisherInfMatrix[2, 2]))
exp(ci)

ci <- optim.vals$par[2] + c(-1, 1) * crit * sqrt(abs(fisherInfMatrix[3, 3]))
exp(ci)

ci <- optim.vals$par[2] + c(-1, 1) * crit * sqrt(abs(fisherInfMatrix[4, 4]))
exp(ci)
#********************************************************************


#***********************Plot*******************************************
tiff("Fig_2_pupalmort.tiff", height = 3, width =5, units = 'in', compression="lzw", res=400)
par(mar=c(3,3,1,1),mgp=c(2,1,0),cex=0.5)



plot(phelps.temp
     ,inst.mort
     ,bty="n"
     ,xlab=expression(paste("Temperature (",degree,"C)"))
     ,ylab="Proportion of pupae failing to emerge"
     ,pch=19
     ,col="blue"
     ,ylim=c(0,1.2)
     ,xlim=c(15,35))
segments(phelps.temp,phelps.b$lower,phelps.temp,phelps.b$upper)


fit.inst.mort <- mort_func(exp(MLEfits[1])
                           ,exp(MLEfits[2])
                           ,exp(MLEfits[3])
                           ,exp(MLEfits[4])
                           ,seq(15,35,0.1))
lines(seq(15,35,0.1),fit.inst.mort,col="grey",xlim=c(15,35),ylim=c(0,1.2))

fit.inst.mort <- mort_func(1.078462e+01                       #from r_4
                           ,1.899046e-01
                           ,5.712021e-06
                           ,1.673068e+00
                           ,seq(15,35,0.1))
lines(seq(15,35,0.1),fit.inst.mort,col="black",xlim=c(15,35),ylim=c(0,1))

legend("topleft",legend=c("Model fit to laboratory data"
                          ,"Model fit to population data")
       ,lty=c(1,1),col=c("grey","black"),bty="n")

dev.off()
#**********************************************************************

pk1 <- exp(MLEfits[1])            # starting params
pk2 <- exp(MLEfits[2])
pk3 <- exp(MLEfits[3])
pk4 <- exp(MLEfits[4])