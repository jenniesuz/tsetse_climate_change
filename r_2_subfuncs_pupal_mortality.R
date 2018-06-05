# This file fits pupal temperature-dependent mortality function
# using laboratory data available in:
#Phelps 1973 The effect of temperature on fat consumption during the puparial stages of
# Glossina morsitans morsitans Westw. (Dipt., Clossinidae) under laboratory conditions,
# and its implication in the field. Bull. Ent. Res. 62, 423-438
# It then plots Figure 2

library("binom")                                            # required pacakges
require(bbmle)

# #**************************Data*************************************
phelps.temp <- c(17,18,20,22,24,25,26,27,28,29,30,31,32)                                 # Vector of temperatures used by Phelps

phelps.mortality <- c(21.9,15.3,10.3,4.9,3.6,5,3.5,10.3,3.8,4,7.2,20.5,45.8)/100         # Proportion of pupae that did not survive to emergence sex unknown

phelps.m.emerge <- c(53,26,27,25,27,25,34,24,14,26,30,48,26)                             # Number of males that emerged

phelps.f.emerge <- c(37,24,25,32,26,32,20,28,12,22,22,49,39)                             # Number of females that emerged

phelps.emerged <- phelps.m.emerge + phelps.f.emerge                                      # Total number that emerged

phelps.number.died <- round(phelps.mortality*phelps.emerged / (1-phelps.mortality),0)    # Estimate number that died

phelps.total <- phelps.emerged + phelps.number.died                                      # Total

phelps.b <- binom.confint(phelps.number.died,phelps.total,methods="exact")               # Confidence intervals


#**************Convert probabilities to instantaneous mortality********
inst.mort <- -log(1-phelps.b$mean)
inst.lower <- -log(1-phelps.b$lower)
inst.upper <- -log(1-phelps.b$upper)
phelps <- cbind.data.frame(temp=phelps.temp,mean=inst.mort,lower=inst.lower,upper=inst.upper)
#***************************************************************

# #****************Set up plot of the data*******************************

# plot
p <- ggplot(phelps, aes(x=temp,y=mean)) +
 xlim(low=16, high=35) +
  geom_pointrange(aes(ymin=lower,ymax=upper,x=temp)) +
  geom_point(size=2) +
   ylim(low=0, high=1.2) +
  labs( y="Proportion of pupae failing to emerge"
        , x=bquote("Temperature"~"(°C)")) + 
  theme_set(theme_bw()) +
  theme( panel.border = element_blank()
         ,axis.line = element_line(color = 'black')
         ,text=element_text(size=12)
         ,plot.margin=unit(c(0.2,0.2,0.2,0.1), "cm")
         ,axis.text=element_text(size=12)
  )
#******************************************************

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
 pk1 <- exp(MLEfits[1])            # starting params
 pk2 <- exp(MLEfits[2])
 pk3 <- exp(MLEfits[3])
 pk4 <- exp(MLEfits[4])
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

#********************Fitted values as a function of temperature***************
mort <- mort_func(pk1,pk2,pk3,pk4,seq(16,35,0.1))

modmort <- cbind.data.frame(temp = seq(16,35,0.1), mort = mort)
#*************************************************************************


#***********************Plot*******************************************
p +
  geom_line(data=modmort
            ,mapping=aes(x=temp,y=mort))



#tiff("Fig_2_pupalmort.tiff", height = 3, width =5, units = 'in', compression="lzw", res=400)

#dev.off()
#**********************************************************************
