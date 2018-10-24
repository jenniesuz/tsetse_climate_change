require(bbmle)                 # required packages
require(deSolve)
library(zoo)
library(MASS)
library(reshape)
library(ggplot2)
library(gridExtra)
library(grid)

#****************Required R scripts*************
source("r_1_data_bioassay.R")                   # tsetse count and temperature data      

source("r_2_subfuncs_pupal_duration.R")         # pupal duration as a function of temperature

source("r_2_subfuncs_larviposition.R")          # larviposition as a function of temperature

source("r_2_subfuncs_pupal_mortality.R")        # pupal mortality as a function of temperature

source("r_2_subfuncs_adult_mortality.R")        # adult mortality as a function of temperature

source("r_3_fit_funcs.R")                       # likelihood function and function to vary which parameters to fit
#***********************************************

#**************Model parameters and conditions***********************
temps <- temps.count$temp[1:length(temps.count$temp)]           # mean temperature for each month - see r_1_data_bioassay.R

times <- seq(from = 0, to = length(temps)*30-1, by = 30)        # time steps to solve model at (in days)

tsetse_params <- function(mu.a.a1=a1      # adult mortality parameters - see r_2_subfuncs
                          ,mu.a.a2=a2
                          ,mu.p.b1=b1     # pupal mortality parameters
                          ,mu.p.b2=b2
                          ,mu.p.b3=b3
                          ,mu.p.b4=b4
                          ,mu.p.b5=b5
                          ,p.dur.c1=c1    # pupal emergence parameters
                          ,p.dur.c2=c2
                          ,p.dur.c3=c3
                          ,lrp.d1=d1      # larviposition parameters
                          ,lrp.d2=d2
                          ,lrn.d3=d3
                          ,lrn.d4=d4
                          , mud.p = 0.00005     # pupal density-dependent mortality
                          , adults.zero = 100   # numbers of adults (parous) at start
                          , juv.zero = 25       # numbers of adults (nulliparous) at start
                          , pupae.zero = 100    # numbers of pupae at start
)
return(as.list(environment()))

initial <- c( P=tsetse_params()$pupae.zero      # initial conditions
              ,J=tsetse_params()$juv.zero
              ,A=tsetse_params()$adults.zero
) 
#**********************************************************************

#********************Population dynamics model*************************
tsetse_mod <- function(tt,yy,parms) with(c(parms,as.list(yy)), {
  
  temp <-  temps[(tt/30)+1]   # take temperature for the time point in the model (tt)
  
  pup <-   pdFunc(p.dur.c1,p.dur.c2,p.dur.c3,temp=temp)                       # pupal emergence rate - see r_2_subfuncs       
  pmort <- pmFunc(mu.p.b1,mu.p.b2,mu.p.b3,mu.p.b4,mu.p.b5,temp=temp)          # pupal mortality rate
  amort <- amFunc(mu.a.a1,mu.a.a2,temp=temp)                                  # adult mortality rate
  larvp <- lrFunc(lrp.d1,lrp.d2,temp=temp)                                    # larviposition rate - parous
  larvn <- lrFunc(lrn.d3,lrn.d4,temp=temp)                        # larviposition rate - nulliparous
  
  # ODEs
  deriv <- rep(NA,2) # empty vector to store ODEs
  
  deriv[1] <-  larvp*A + larvn*J  - pup*P - mud.p*P*P - pmort*P   # change in pupae
  
  deriv[2] <-  pup*P/2 - amort*J - J*(larvn)                      # change in nulliparous adults
  
  deriv[3] <-  J*(larvn) - amort*A                                # change in parous adults 
  
  return(list(deriv))
})
#**********************************************************************

#**************Function to simulate model******************************
simPop <- function(init=initial, tseq = times, modFunction=tsetse_mod      # function which runs the model
                   , parms = tsetse_params()) { 
  simDat <- as.data.frame(lsoda(init, tseq, modFunction, parms=parms))
  return(simDat)
}
#*********************************************************************
test <- simPop(parms=tsetse_params(mud.p = 0.00006,mu.p.b3=1.481-0.681*1.96))                      # test run with default parameter estimates from available data

temps.count$A <- test$A+test$J                                            # add to dataframe
plot.dat <- temps.count[65:length(temps.count$time),]                     # remove starting runs

ggplot(plot.dat, aes(x=time,y=count)) +                                   # plot
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
#*********************************************************************


#***********fitting function********************

mle.func <- function(param,param.val,dd.coef=0.0001){
  parms<-tsetse_params()
  parms[param]<-param.val
  
init.pars <- c(
  log_mud.p=log(dd.coef)
  ,log_mu.a.a1=log(as.numeric(a1))
  ,log_mu.a.a2=log(as.numeric(a2))
)
trace <- 3
optim.vals <- optim(par = init.pars
                    , objFXN
                    , fixed.params = parms
                    , dat = temps.count
                    , control = list(trace = trace, maxit = 200)
                    , method = "SANN")

exp(optim.vals$par) #
optim.vals <- optim(par = optim.vals$par
                    , objFXN
                    , fixed.params = parms
                    , dat = temps.count
                    , control = list(trace = trace, maxit = 200)
                    , method = "SANN")
exp(optim.vals$par) #
optim.vals <- optim(par = optim.vals$par
                       , objFXN
                       , fixed.params = parms
                       , dat = temps.count
                       , control = list(trace = trace, maxit = 1000, reltol = 10^-7)
                       , method = "Nelder-Mead" # 
                       , hessian = T)
MLEfits <- optim.vals$par 
AIC <- -2*-optim.vals$value+ 2*3       

#*****************confidence intervals for parameter estimates*********
fisherInfMatrix <- solve(optim.vals$hessian) ## invert the Hessian, to estimate the covar-var matrix of parameter estimates
fisherInfMatrix

# Finds the critical z value
conf.level <- 0.95
crit <- qnorm((1 + conf.level)/2)

ci1 <- optim.vals$par[1] + c(-1, 1) * crit * sqrt(abs(fisherInfMatrix[1, 1]))
ci2 <- optim.vals$par[2] + c(-1, 1) * crit * sqrt(abs(fisherInfMatrix[2, 2]))
ci3 <- optim.vals$par[3] + c(-1, 1) * crit * sqrt(abs(fisherInfMatrix[3,3]))
return(list(MLEfits,AIC,ci1,ci2,ci3))

}

b1min <- mle.func(3,0.0019-0.0004*1.96)
b1max <- mle.func(3,0.0019+0.0004*1.96,0.00008)

b3min <- mle.func(5,1.481-0.681*1.96,0.00006)
b3max <- mle.func(5,1.481+0.681*1.96,0.00008)

b5min <- mle.func(7,1.211-0.117*1.96,0.00008)
b5max<- mle.func(7,1.211+0.117*1.96,0.00008)


