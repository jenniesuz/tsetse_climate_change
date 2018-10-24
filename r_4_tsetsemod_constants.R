require(bbmle)  # required packages
require(deSolve)
library(zoo)
library(beepr)
library(MASS)
library(reshape)
library(ggplot2)
library(gridExtra)
library(grid)

#****************Required R scripts*************
source("r_1_data_bioassay.R")                   # count and climate data      

source("r_2_subfuncs_pupal_duration.R")

source("r_2_subfuncs_larviposition.R")

source("r_2_subfuncs_pupal_mortality.R")        # initial estimates for pupal temperature-dependent mortality

source("r_2_subfuncs_adult_mortality.R")

source("r_3_fit_funcs.R")                       # likelihood function
#***********************************************

#********************Model parameters and conditions*************************************
temps <- temps.count$temp[1:length(temps.count$temp)]         # temperature data for each month
 
times <- seq(from = 0, to = length(temps)*30-1, by = 30)        # time steps to solve model at (in days)

tsetse_params <- function(mu.a.a1=a1      # adult mortality parameters
                          ,mu.a.a2=a2
                          ,mu.p.b1=b1     # pupal mortality parameters
                          ,mu.p.b2=b2
                          ,mu.p.b3=b3
                          ,mu.p.b4=b4
                          ,mu.p.b5=b5
                          ,pupemerg=0.03/2    # pupal emergence parameters
                          ,larvipos=0.1        # larviposition  
                          , mud.p = 0.00005     ## pupal density-dependent mortality
                          , adults.zero = 100   ## numbers of adults at start
                          , juv.zero=25
                          , pupae.zero = 100    ## numbers of pupae at start
                          , tempdep="p"
)
return(as.list(environment()))

# Initial conditions. A vector containing the numbers of starting pupae P and adults A.
initial <- c( P=tsetse_params()$pupae.zero   
              ,J=tsetse_params()$juv.zero
              ,A=tsetse_params()$adults.zero
) 
#******************************************************************************

#********************Population dynamics model***********************************
tsetse_mod <- function(tt,yy,parms) with(c(parms,as.list(yy)), {
  
  temp <-  temps[(tt/30)+1]                     # the temperature data is read in by another script. This takes the temperature for the time point in the model (tt)
  
  pup <-   pupemerg                              # pupal emergence rate           
  
  if(tempdep=="p"){
    pmort <- pmFunc(mu.p.b1,mu.p.b2,mu.p.b3,mu.p.b4,mu.p.b5,temp=temp)  # pupal mortality rate
    amort <- mu.a.a1   # adult mortality rate
    
  }else{
    if(tempdep=="a"){
  pmort <- mu.p.b1  # pupal mortality rate
  amort <- amFunc(mu.a.a1,mu.a.a2,temp=temp)    # adult mortality rate
    }else{
      if(tempdep=="n"){
      pmort <- mu.p.b1  # pupal mortality rate
      amort <- mu.a.a1    # adult mortality rate
      }
    }
  }
  
  larva <- larvipos   # larviposition rate
  larvt <- larvipos
  
  prob.surv <- exp(-pmort)^(1/pup) #exp(-pmort*(1/pup)) #exp(-pmort)^(1/pup)
  # ODEs
  deriv <- rep(NA,2) # There are two ODEs so creating an empty vector to store them in here
  
  deriv[1] <-  larva*A + larvt*J  - pup*P - mud.p*P*P #- pmort*P          # change in pupae
  
  deriv[2] <-  pup*P/2*prob.surv - amort*J - J*(larvt)                             # change in nulliparous
  
  deriv[3] <-  J*(larvt) - amort*A                                         # change in adults 
  
  return(list(deriv))
})
#*******************************************************

#**************SIMULATE***************************************
simPop <- function(init=initial, tseq = times, modFunction=tsetse_mod      # function which runs the model
                   , parms = tsetse_params(tempdep=="p")) { 
  simDat <- as.data.frame(lsoda(init, tseq, modFunction, parms=parms))
  return(simDat)
}
########################################################

test <- simPop(parms=tsetse_params(mud.p=0.0003, mu.p.b1=0.006,mu.a.a1=0.01
,tempdep="n") )

temps.count$A <- test$A+test$J
plot.dat <- temps.count[65:length(temps.count$time),] # remove starting runs

ggplot(plot.dat, aes(x=time,y=count)) +
  geom_line(mapping=aes(x=time,y=A),size=0.2,col="black") +
  geom_point(size=0.4,col="red") +
  labs( y= "Numbers of tsetse"
        , x="Date") + 
  scale_x_yearmon(
    limits = as.yearmon(c('1960-02','2017-12'))) +
  scale_y_continuous(trans='log10',limits=c(0.1, high=200)) +
  theme_set(theme_bw()) +
  theme( panel.border = element_blank()
         ,axis.line = element_line(color = 'black')
         ,text=element_text(size=7)
         ,plot.margin=unit(c(0.3,0.7,0.3,0.2), "cm")
         ,axis.text=element_text(size=6)
  )
#****************************************************************



#****************************************************************************************
#*************fit 1 all constant********************************
init.pars.none <- c(
  log_mud.p=log(0.0003)
  ,log_mu.a.a1=log(0.01)
  ,log_mu.p.b1=log(0.006)     # pupal mortality parameters
)
#**************Optimise*******************************************
trace <- 3
optim.vals.none  <- optim(par = init.pars.none 
                       , objFXN
                       , fixed.params = tsetse_params(tempdep="n")
                       , dat = temps.count
                       , control = list(trace = trace, maxit = 200)
                       , method = "SANN")

exp(optim.vals.none$par) #
beep()
optim.vals.none  <- optim(par = optim.vals.none$par
                       , objFXN
                       , fixed.params = tsetse_params(tempdep="n")
                       , dat = temps.count
                       , control = list(trace = trace, maxit = 200)
                       , method = "SANN")
exp(optim.vals.none$par) #
beep()

optim.vals.none <- optim(par = optim.vals.none$par
                       , objFXN
                       , fixed.params =  tsetse_params(tempdep="n")
                       , dat = temps.count
                       , control = list(trace = trace, maxit = 2000, reltol = 10^-7)
                       , method = "Nelder-Mead" # 
                       , hessian = T)
optim.vals.none # convergence 0 means algorithm converged
beep()

MLEfits.none <- optim.vals.none$par 
exp(MLEfits.none)
beep()

-2*-optim.vals.none$value+ 2*3  #6762.134

#*****************confidence intervals for parameter estimates*********
 fisherInfMatrix <-solve(optim.vals.none$hessian) ## invert the Hessian, to estimate the covar-var matrix of parameter estimates
fisherInfMatrix

# Finds the critical z value
conf.level <- 0.95
crit <- qnorm((1 + conf.level)/2)

ci.n1 <- optim.vals.none$par[1] + c(-1, 1) * crit * sqrt(abs(fisherInfMatrix[1, 1]))
exp(ci.n1)
ci.n2 <- optim.vals.none$par[2] + c(-1, 1) * crit * sqrt(abs(fisherInfMatrix[2, 2]))
exp(ci.n2)
ci.n3 <- optim.vals.none$par[3] + c(-1, 1) * crit * sqrt(abs(fisherInfMatrix[3,3]))
exp(ci.n3)
#**********************************************************


# #*********************Fitted model*******************************
sim <- simPop(parms=tsetse_params(
  mud.p=exp(MLEfits.none)[1]
  ,pupemerg = exp(MLEfits.none)[2]
  ,larvipos = exp(MLEfits.none)[3]
  , mu.a.a1=exp(MLEfits.none)[4]
  , mu.p.b1=exp(MLEfits.none)[5]
  ,tempdep="n"
  
))

temps.count$A <- sim$A+sim$J
plot.dat <- temps.count[65:length(temps.count$time),] # remove starting runs

#****************Plots************************************************
tiff("Fig_SXa.tiff", height =4, width = 5, units = 'in', compression="lzw", res=400)
ggplot(plot.dat, aes(x=time,y=count)) +
  geom_line(mapping=aes(x=time,y=A),size=0.2,col="black") +
  geom_point(size=0.4,col="red") +
  labs( y= "Numbers of tsetse"
        , x="Date") + 
  scale_x_yearmon(
    limits = as.yearmon(c('1960-02','2017-12'))) +
  scale_y_continuous(trans='log10',limits=c(0.1, high=200)) +
  theme_set(theme_bw()) +
  theme( panel.border = element_blank()
         ,axis.line = element_line(color = 'black')
         ,text=element_text(size=7)
         ,plot.margin=unit(c(0.3,0.7,0.3,0.2), "cm")
         ,axis.text=element_text(size=6)
  )
dev.off()
#********************************************************************************
#****************************************************************************************





#*************fit 2 pupal temp.dep. const. adult********************************
init.pars.pm <- c(
  log_mud.p=log(0.00002)
  ,log_mu.a.a1=log(0.01)
  ,log_mu.p.b1=log(as.numeric(b1))     # pupal mortality parameters
  ,log_mu.p.b3=log(as.numeric(b3)) 
  ,log_mu.p.b5=log(as.numeric(b5)) 
)
#**************Optimise*******************************************
trace <- 3
optim.vals.pm <- optim(par = init.pars.pm
                       , objFXN
                       , fixed.params = tsetse_params(tempdep="p")
                       , dat = temps.count
                       , control = list(trace = trace, maxit = 200)
                       , method = "SANN")

exp(optim.vals.pm$par) #
beep()
optim.vals.pm <- optim(par = optim.vals.pm$par
                       , objFXN
                       , fixed.params = tsetse_params(tempdep="p")
                       , dat = temps.count
                       , control = list(trace = trace, maxit = 200)
                       , method = "SANN")
exp(optim.vals.pm$par) #
beep()

optim.vals.pm <- optim(par = optim.vals.pm$par
                       , objFXN
                       , fixed.params = tsetse_params(tempdep="p")
                       , dat = temps.count
                       , control = list(trace = trace, maxit = 2000, reltol = 10^-7)
                       , method = "Nelder-Mead" # 
                       , hessian = T)
optim.vals.pm # convergence 0 means algorithm converged
beep()

MLEfits.pm <- optim.vals.pm$par 
exp(MLEfits.pm)
beep()

-2*-optim.vals.pm$value+ 2*3  # 2518

#*****************confidence intervals for parameter estimates*********
fisherInfMatrix <- solve(optim.vals.pm$hessian) ## invert the Hessian, to estimate the covar-var matrix of parameter estimates
fisherInfMatrix

# Finds the critical z value
conf.level <- 0.95
crit <- qnorm((1 + conf.level)/2)

ci.pm1 <- optim.vals.pm$par[1] + c(-1, 1) * crit * sqrt(abs(fisherInfMatrix[1, 1]))
exp(ci.pm1)
ci.pm2 <- optim.vals.pm$par[2] + c(-1, 1) * crit * sqrt(abs(fisherInfMatrix[2, 2]))
exp(ci.pm2)
ci.pm3 <- optim.vals.pm$par[3] + c(-1, 1) * crit * sqrt(abs(fisherInfMatrix[3,3]))
exp(ci.pm3)
ci.pm4 <- optim.vals.pm$par[4] + c(-1,1)*crit * sqrt(abs(fisherInfMatrix[4,4]))
exp(ci.pm4)
ci.pm5 <- optim.vals.pm$par[5] + c(-1,1)*crit *sqrt(abs(fisherInfMatrix[5,5]))
exp(ci.pm5)
#**********************************************************************


#*********************fit 3 vary adult mortality rate function********************************
init.pars.am <- c(
  log_mud.p=log(0.00001)
  ,log_mu.a.a1=log(as.numeric(a1))
  ,log_mu.a.a2=log(as.numeric(a2))
  ,log_mu.p.b1=log(as.numeric(b1))
)
#**************Optimise*******************************************
trace <- 3
optim.vals.am <- optim(par = init.pars.am
                    , objFXN
                    , fixed.params = tsetse_params(tempdep="a")
                    , dat = temps.count
                    , control = list(trace = trace, maxit = 200)
                    , method = "SANN")

exp(optim.vals.am$par) #
beep()
optim.vals.am <- optim(par = optim.vals.am$par
                    , objFXN
                    , fixed.params = tsetse_params(tempdep="a")
                    , dat = temps.count
                    , control = list(trace = trace, maxit = 200)
                    , method = "SANN")
exp(optim.vals.am$par) #
beep()

optim.vals.am <- optim(par = optim.vals.am$par
                    , objFXN
                    , fixed.params = tsetse_params(tempdep="a")
                    , dat = temps.count
                    , control = list(trace = trace, maxit = 2000, reltol = 10^-7)
                    , method = "Nelder-Mead" # 
                    , hessian = T)
optim.vals.am # convergence 0 means algorithm converged
beep()

MLEfits.am <- optim.vals.am$par 
exp(MLEfits.am)
beep()

-2*-optim.vals.am$value+ 2*3 # 2523

#*****************confidence intervals for parameter estimates*********
fisherInfMatrix <- solve(optim.vals.am$hessian) ## invert the Hessian, to estimate the covar-var matrix of parameter estimates
fisherInfMatrix

# Finds the critical z value
conf.level <- 0.95
crit <- qnorm((1 + conf.level)/2)

ci <- optim.vals.am$par[1] + c(-1, 1) * crit * sqrt(abs(fisherInfMatrix[1, 1]))
exp(ci)
ci <- optim.vals.am$par[2] + c(-1, 1) * crit * sqrt(abs(fisherInfMatrix[2, 2]))
exp(ci)
ci <- optim.vals.am$par[3] + c(-1, 1) * crit * sqrt(abs(fisherInfMatrix[3,3]))
exp(ci)
ci <- optim.vals.am$par[4] + c(-1, 1) * crit * sqrt(abs(fisherInfMatrix[4,4]))
exp(ci)


#**********************************************************************

# #*********************Fitted model*******************************
sim <- simPop(parms=tsetse_params(
                                   mud.p=exp(MLEfits.am)[1]
                                  , mu.a.a1=exp(MLEfits.am)[2]
                                  , mu.a.a2=exp(MLEfits.am)[3]
                                  ,mu.p.b1=exp(MLEfits.am)[4]
                                  ,tempdep="a"
                                 
))



temps.count$A <- sim$A+sim$J
plot.dat <- temps.count[65:length(temps.count$time),] # remove starting runs

#****************Plots************************************************
tiff("Fig_SXc.tiff", height =4, width = 5, units = 'in', compression="lzw", res=400)

ggplot(plot.dat, aes(x=time,y=count)) +
  geom_line(mapping=aes(x=time,y=A),size=0.2,col="black") +
  geom_point(size=0.4,col="red") +
  labs( y= "Numbers of tsetse"
        , x="Date") + 
  scale_x_yearmon(
               limits = as.yearmon(c('1960-02','2017-12'))) +
  scale_y_continuous(trans='log10',limits=c(0.1, high=200)) +
  theme_set(theme_bw()) +
  theme( panel.border = element_blank()
         ,axis.line = element_line(color = 'black')
         ,text=element_text(size=7)
         ,plot.margin=unit(c(0.3,0.7,0.3,0.2), "cm")
         ,axis.text=element_text(size=6)
  )

dev.off()
#********************************************************************************

