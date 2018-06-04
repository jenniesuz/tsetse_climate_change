require(bbmle)  # required packages
require(deSolve)
library(zoo)

#****************Required R scripts*************
source("r_1_data_bioassay.R")                   # count and climate data      

source("r_3_subfuncs_pupal_mortality.R")        # initial estimates for pupal temperature-dependent mortality

source("r_2_fit_funcs.R")                       # likelihood function
#***********************************************

#********************Model parameters and conditions*************************************
temps <- b.temps$MeanC # temperature data for each day
length(temps) 

times <- seq(from = 1, to = length(temps),by=1)        # time steps to solve model at (in days)

tsetse_params <- function(   larv = 1/10          ## rate larvae are produced (days)
                             , pup.a = 17.94          ## rate pupae emerge as adults (days) - if not temperature dependent
                             , pup.b = 82.3
                             , pup.c = -0.253
                             , pup.t1 = 16
                             , nboxes=3
                             , mud.p = 0.0000     ## pupal density-dependent mortality
                             , mu.p.k1 = exp(2.03)     ## parameters for pupal temperature-dependent mortality
                             , mu.p.k2 = pk2
                             , mu.p.k3 = exp(-7.83)
                             , mu.p.k4 = pk4
                             , mu.a.k1 = 0      ## parameter 1 adult temperature-dependent mortality
                             , mu.a.k2 = 0        ## parameter 2 adult temperature-dependent mortality
                             , adults.zero = 100 ## numbers of adults at start
                             , pupae.zero = 100 ## numbers of pupae at start
)
  return(as.list(environment()))

# Initial conditions. A vector containing the numbers of starting pupae P and adults A.
initial <- c( P=tsetse_params()$pupae.zero/tsetse_params()$nboxes 
              ,rep(tsetse_params()$pupae.zero/tsetse_params()$nboxes,tsetse_params()$nboxes-1)
              ,A=tsetse_params()$adults.zero
) 

#******************************************************************************



#********************Population dynamics model***********************************
tsetse_mod <- function(tt,yy,parms) with(c(parms,as.list(yy)), {
  
  temp <-  temps[tt]                                          # the temperature data is read in by another script. This takes the temperature for the time point in the model (tt)
  
  pup <- (1/(pup.a + pup.b*exp(pup.c*(temp-pup.t1))))/nboxes  # pupal development rate for day tt
  
  # calculate how far back to go to average temperature for proportion of pupae surviving
  store.days <- 1                                             # store the total pupal duration
  day <- tt                                                   # time t in the model
  cum.rate <- pup.a + pup.b*exp(pup.c*(temp-pup.t1))          # start off with the rate of pupal development for that day
  temperature <- numeric(0)                                   # empty vector to store temperatures
  temperature[1] <- temp                                      # first value is temperature at day tt
  while (cum.rate < 1) {
    store.days <- store.days + 1                              # update store.days
    day <- day - 1                                            # look at previous day
    if ( day <= 0 ) {                                         # if before model start assign to 25 - running model for a year before the actual data means this doesn't matter
      day.temp <- 25
    }else{
      day.temp <- temps[day]                                  # otherwise set next temp to previous days temp
      }                    
    temperature[store.days] <- day.temp                       # store in vector of temperatures
    cum.rate <- cum.rate + (pup.a + pup.b*exp(pup.c*(day.temp-pup.t1)))   # increase development
  }
  mean.temp <- mean(temperature)                              # now calculate mean for that period for input into pupal survival as function of temperature

  inst.mort.pup <- mort_func(k1=mu.p.k1                       # calculate pupal mortality rate at the given temperature
                             ,k2=mu.p.k2
                             ,k3=mu.p.k3
                             ,k4=mu.p.k4
                             ,temp=mean.temp)                 # note using the above mean.temp here
  
  prob.surv.pup <- exp(-inst.mort.pup)                        # convert to a probability
  
  if (temp > 24){                                             # calculate adult temperature-dependent mortality
    mu.at <- mu.a.k1*exp(mu.a.k2*(temp-24))
  }else{
    mu.at <- mu.a.k1
  }

  # ODEs
  deriv <- rep(NA,length(yy)) # creating an empty vector to store ODEs in here
  
  deriv[-(nboxes+1)] <- - (mud.p*sum(yy[-(nboxes+1)]))*yy[-(nboxes+1)] - pup*yy[-(nboxes+1)] # all compartments rate at which they leave nboxes subcompartments of the pupal stage
  
  deriv[1] <-  deriv[1] + larv*A  # rate at which pupae enter first subcompartment
  
  deriv[2:nboxes] <- deriv[2:nboxes] + pup*yy[-c(nboxes,nboxes+1)] # rate at which pupae enter other compartments
  
  deriv[nboxes+1] <-  yy[length(yy)-1]*pup/2*(prob.surv.pup) - mu.at*A  # change in adults 
  
  return(list(deriv))
  
  
})
#*******************************************************

#**************SIMULATE***************************************
simPop <- function(init=initial, tseq = times, modFunction=tsetse_mod      # function which runs the model
                   , parms = tsetse_params()) { 
  simDat <- as.data.frame(lsoda(init, tseq, modFunction, parms=parms))
  return(simDat)
}
########################################################

# # #************initial parameters**********************
# #**********************************************
# init.pars <- c(
#   log_mu.p.k2=log(2.212627e-01)
#   ,log_mu.p.k4=log(1.091462e+00)
#   ,log_mu.a.k1=log(3.968048e-02)
#   ,log_mu.a.k2=log(6.405333e-02)
#   ,log_mud.p=log(6.726220e-06)
# )
# 
# #**************Optimise*******************************************
# trace <- 3
# optim.vals <- optim(par = init.pars
#                     , objFXN
#                     , fixed.params = tsetse_params()
#                     , dat = temps.count
#                     , control = list(trace = trace, maxit = 200)
#                     , method = "SANN")
# 
# 
# 
# exp(optim.vals$par) # 
# optim.vals <- optim(par = optim.vals$par
#                     , objFXN
#                     , fixed.params = tsetse_params()
#                     , dat = temps.count
#                     , control = list(trace = trace, maxit = 200)
#                     , method = "SANN")
# exp(optim.vals$par) # 
# 
# 
# # want to then feed output from SANN to Nelder-Mead
# optim.vals <- optim(par = init.pars #optim.vals$par#optim.vals$par
#                     , objFXN
#                     , fixed.params = tsetse_params()
#                     , dat = temps.count
#                     , control = list(trace = trace, maxit = 1000, reltol = 10^-7)
#                     , method = "Nelder-Mead" # 
#                     , hessian = T)
# optim.vals # convergence 0 means algorithm converged
# 
# MLEfits <- optim.vals$par 
# exp(MLEfits)
# #******************************************************************
# 
# 
# #*****************confidence intervals for parameter estimates*********
# fisherInfMatrix <- solve(optim.vals$hessian) ## invert the Hessian, to estimate the covar-var matrix of parameter estimates
# fisherInfMatrix
# 
# # Finds the critical z value
# conf.level <- 0.95
# crit <- qnorm((1 + conf.level)/2)
# 
# ci <- optim.vals$par[1] + c(-1, 1) * crit * sqrt(abs(fisherInfMatrix[1, 1]))
# exp(ci)
# ci <- optim.vals$par[2] + c(-1, 1) * crit * sqrt(abs(fisherInfMatrix[2, 2]))
# exp(ci)
# ci <- optim.vals$par[3] + c(-1, 1) * crit * sqrt(abs(fisherInfMatrix[3,3]))
# exp(ci)
# ci <- optim.vals$par[4] + c(-1, 1) * crit * sqrt(abs(fisherInfMatrix[4, 4]))
# exp(ci)
# 
# ci <- optim.vals$par[5] + c(-1, 1) * crit * sqrt(abs(fisherInfMatrix[5, 5]))
# exp(ci)
# 
# #**********************************************************************
# 
# 
# # #********************************************************************
# test <- simPop(parms=tsetse_params(
#                                    mu.p.k2=exp(MLEfits[1])
#                                    ,mu.p.k4=exp(MLEfits[2])
#                                    ,mu.a.k1=exp(MLEfits[3])
#                                    ,mu.a.k2=exp(MLEfits[4])
#                                    ,mud.p=exp(MLEfits[5])
# ))
# 
# 
# # #********************************************************************


#*****TEST**************
test <- simPop(parms=tsetse_params(
  mu.p.k2=pk2
  ,mu.p.k4=pk4
  ,mu.a.k1=0.03
  ,mu.a.k2=0.05
  ,mud.p=0.0000003
))
#*************

plot(b.temps$Date
  ,test$A
  ,bty="n"
  ,type="l"
  ,xlab="Date"
  ,ylab="Numbers of adult tsetse")

#****************Plots************************************************
#tiff("Fig_2.tiff", height = 3, width = 5, units = 'in', compression="lzw", res=400)
# par(mar=c(3,3,1,1),mgp=c(2,1,0),cex=0.5)
# plot( 
#       temps.count$time[65:length(temps.count$temp)]
#       ,temps.count$count[65:length(temps.count$temp)]
#      ,col="black"
#      ,cex.lab=1.2
#      ,type="n"
#      ,bty="n"
#      ,xlab="Date"
#      ,xlim=c(as.yearqtr('1959 Q4'),
#              as.yearqtr('2017 Q4'))
#     # ,xlim=c(as.yearqtr('1990 Q4'),
#      #        as.yearqtr('2017 Q4'))
#     ,log="y"
#      ,pch=1
#      ,ylab="Numbers of tsetse"
#      ,ylim=c(0.1,100)#,ylim=c(0.1,150)
# )
# points(temps.count$time[65:length(temps.count$temp)]
#        ,temps.count$count[65:length(temps.count$temp)]
#        ,cex=0.6
#        ,col="red"
#        ,pch=19)
# par(new=T)
# plot(temps.count$time[65:length(test$A)]
#      #,model.output
#      ,test$A[65:length(test$A)]
#      ,cex.lab=1.2
#      ,type="l"
#      ,col="black"
#      ,bty="n"
#      ,xaxt="n"
#      ,yaxt="n"
#       ,xlim=c(as.yearqtr('1959 Q4'),
#               as.yearqtr('2017 Q4'))
#     # ,xlim=c(as.yearqtr('1990 Q4'),
#    #          as.yearqtr('2017 Q4'))
#     ,log="y"
#      ,xlab=" "
#      ,ylab=" "
#      ,lty=1
#     ,ylim=c(0.1,100)
#    # ,ylim=c(0.1,150)
# )
# 
# 
# legend("topright",title="Numbers of tsetse",legend=c("Observed","Model")
#        ,lty=c(NA,1,3),pch=c(19,NA),col=c("red","black"),bty="n")
#dev.off()
#********************************************************************************

# 
# 
# temps <- rep(25,length(1:length(temps.count$temp)))         # temperature data for each month
# 
# #***************run model with fitted parameters*********************
# test <- simPop(parms=tsetse_params(mu.p.k1=exp(MLEfits[1])
#                                    ,mu.p.k2=exp(MLEfits[2])
#                                    ,mu.p.k3=exp(MLEfits[3])
#                                    ,mu.p.k4=exp(MLEfits[4])
#                                    ,mu.a.k1=exp(MLEfits[5])
#                                    ,mu.a.k2=exp(MLEfits[6])
#                                    ,mud.p=exp(MLEfits[7])
# ))
# 
# plot(temps.count$time[1:length(temps.count$temp)]
#      #,model.output
#      ,test$A
#      #,test.smooth
#      ,type="l"
#      ,col="black"
#      ,bty="n"
#      ,xaxt="n"
#      ,yaxt="n"
#      ,xlim=c(as.yearqtr('1950 q1'),
#              as.yearqtr('2017 Q4'))
#      #,log="y"
#      ,xlab=" "
#      ,ylab=" "
#      ,lty=3
#      ,ylim=c(0.1,500)
# )
# 
# 
# 
