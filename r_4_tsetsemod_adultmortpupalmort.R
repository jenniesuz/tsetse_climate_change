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
test <- simPop(parms=tsetse_params(mud.p = 0.00006))                      # test run with default parameter estimates from available data

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

#****fit dd - only vary density-dependent mortality coefficient******
#************initial parameter estimate**********************
init.pars.dd <- c(
  log_mud.p=log(0.00006)
)
#**************Optimise*******************************************
 trace <- 3
# optim.vals.dd <- optim(par = init.pars.dd
#                     , objFXN
#                     , fixed.params = tsetse_params()
#                     , dat = temps.count
#                     , control = list(trace = trace, maxit = 200)
#                     , method = "SANN")
# exp(optim.vals.dd$par) # 
# optim.vals.dd <- optim(par = optim.vals.dd$par
#                     , objFXN
#                     , fixed.params = tsetse_params()
#                     , dat = temps.count
#                     , control = list(trace = trace, maxit = 200)
#                     , method = "SANN")
# exp(optim.vals.dd$par) 
init.pars.dd <- c(
  log_mud.p=log(8.453963e-05)          # from 2 SANN runs above
)
optim.vals.dd <- optim(par = init.pars.dd                  # fit using Nelder-Mead
                    , objFXN
                    , fixed.params = tsetse_params()
                    , dat = temps.count
                    , control = list(trace = trace, maxit = 2000, reltol = 10^-7)
                    , method = "Nelder-Mead" # 
                    , hessian = T)
optim.vals.dd           
MLEfits.dd <- optim.vals.dd$par 
exp(MLEfits.dd)
- 2 * -optim.vals.dd$value + 2 * 1 #2867
#*********fit 2 vary adult mortality rate function**********
init.pars.am <- c(
  log_mud.p=log(0.00006)
  ,log_mu.a.a1=log(as.numeric(a1))
  ,log_mu.a.a2=log(as.numeric(a2))
)
#**************Optimise*******************************************
# trace <- 3
# optim.vals.am <- optim(par = init.pars.am
#                     , objFXN
#                     , fixed.params = tsetse_params()
#                     , dat = temps.count
#                     , control = list(trace = trace, maxit = 200)
#                     , method = "SANN")
# 
# exp(optim.vals.am$par) # 
# beep()
# optim.vals.am <- optim(par = optim.vals.am$par
#                     , objFXN
#                     , fixed.params = tsetse_params()
#                     , dat = temps.count
#                     , control = list(trace = trace, maxit = 200)
#                     , method = "SANN")
# exp(optim.vals.am$par) # 
# beep()
init.pars.am <- c(
  log_mud.p=log(2.030794e-05)
  ,log_mu.a.a1=log(3.367511e-02)
  ,log_mu.a.a2=log(1.158975e-01)
)


init.pars.am <- c(
  log_mud.p=log(2.009573e-05)
  ,log_mu.a.a1=log(3.368781e-02)
  ,log_mu.a.a2=log(1.167421e-01)
)
optim.vals.am <- optim(par = init.pars.am
                    , objFXN
                    , fixed.params = tsetse_params()
                    , dat = temps.count
                    , control = list(trace = trace, maxit = 1000, reltol = 10^-7)
                    , method = "Nelder-Mead" # 
                    , hessian = T)
optim.vals.am # convergence 0 means algorithm converged
beep()

MLEfits.am <- optim.vals.am$par 
exp(MLEfits.am)
beep()
-2*-optim.vals.am$value+ 2*3 # 1606
#********************************************************************


#*********fit 2 vary pupal mortality rate function******************
init.pars.pm <- c(
  log_mud.p=log(0.00006)
  ,log_mu.p.b1=log(as.numeric(b1))
  ,log_mu.p.b3=log(as.numeric(b3))
  ,log_mu.p.b5=log(as.numeric(b5))
)
#**************Optimise*******************************************
# trace <- 3
# optim.vals.pm <- optim(par = init.pars.pm
#                     , objFXN
#                     , fixed.params = tsetse_params()
#                     , dat = temps.count
#                     , control = list(trace = trace, maxit = 200)
#                     , method = "SANN")
# 
# exp(optim.vals.pm$par) #
# optim.vals.pm <- optim(par = optim.vals.pm$par
#                     , objFXN
#                     , fixed.params = tsetse_params()
#                     , dat = temps.count
#                     , control = list(trace = trace, maxit = 200)
#                     , method = "SANN")
# exp(optim.vals.pm$par) #
init.pars.pm <- c(
  log_mud.p=log(0.000028951)           # from 2 SANN runs
  ,log_mu.p.b1=log(0.008642859)
  ,log_mu.p.b3=log(2.116735401)
  ,log_mu.p.b5=log(2.883192124)
)
optim.vals.pm <- optim(par = init.pars.pm
                       , objFXN
                       , fixed.params = tsetse_params()
                       , dat = temps.count
                       , control = list(trace = trace, maxit = 2000, reltol = 10^-7)
                       , method = "Nelder-Mead" # 
                       , hessian = T)
optim.vals.pm 
MLEfits.pm <- optim.vals.pm$par 
exp(MLEfits.pm)
-2*-optim.vals.pm$value+ 2*4 # 1789
#**********************************************************************
#*********fit 2 vary both mortality rate functions******************
init.pars.pam <- c(
  log_mud.p=log(0.00006)
  ,log_mu.a.a1=log(as.numeric(a1))
  ,log_mu.a.a2=log(as.numeric(a2))
  ,log_mu.p.b1=log(as.numeric(b1))
  ,log_mu.p.b3=log(as.numeric(b3))
  ,log_mu.p.b5=log(as.numeric(b5))
)
#**************Optimise*******************************************
trace <- 3
optim.vals.pam <- optim(par = init.pars.pam
                    , objFXN
                    , fixed.params = tsetse_params()
                    , dat = temps.count
                    , control = list(trace = trace, maxit = 200)
                    , method = "SANN")

exp(optim.vals.pam$par) #
optim.vals.pam <- optim(par = optim.vals.pam$par
                    , objFXN
                    , fixed.params = tsetse_params()
                    , dat = temps.count
                    , control = list(trace = trace, maxit = 200)
                    , method = "SANN")
exp(optim.vals.pam$par) #
init.pars.pam <- c(
  log_mud.p=log(3.542124e-05)           # from 2 SANN runs
  ,log_mu.a.a1=log(3.116092e-02)
  ,log_mu.a.a2=log(1.374914e-01)
  ,log_mu.p.b1=log(2.645885e-03)
  ,log_mu.p.b3=log(2.106021e+00)
  ,log_mu.p.b5=log(1.455521e+00)
)
optim.vals.pam <- optim(par = init.pars.pam
                       , objFXN
                       , fixed.params = tsetse_params()
                       , dat = temps.count
                       , control = list(trace = trace, maxit = 2000, reltol = 10^-7)
                       , method = "Nelder-Mead" # 
                       , hessian = T)
optim.vals.pam 
MLEfits.pam <- optim.vals.pam$par 
exp(MLEfits.pam)
-2*-optim.vals.pam$value+ 2*6 # 1764
#**********************************************************************


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
#*********************Fitted model*******************************
sim <- simPop(parms=tsetse_params(
                                   mud.p=2.356929e-05   #exp(MLEfits.am)[1]
                                  , mu.a.a1=3.365351e-02 #exp(MLEfits.am)[2]
                                  , mu.a.a2=1.167739e-01 #exp(MLEfits.am)[3]
                                 
))

temps.count$A <- sim$A+sim$J
plot.dat <- temps.count[65:length(temps.count$time),] # remove starting runs

#******************************Figure 4**************************
tiff("Fig_4.tiff", height =4, width = 5, units = 'in', compression="lzw", res=400)

ggplot(plot.dat, aes(x=time,y=count)) +
  geom_point(size=0.4,col="red") +
  geom_line(mapping=aes(x=time,y=A),size=0.2,col="black") +
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
#************Figure 3*****************************************
tiff("Fig_3.tif", height = 5, width = 5, units = 'in', compression="lzw", res=500)

args.list <- c(list(amPlotFit
                  ,pmPlotFit
                  ,pdPlotFit
                  ,lrPlotFit)
                  ,list(ncol=2
                        ,nrow=2
                        ,bottom=textGrob(bquote("Temperature"~"(?C)")
                        , gp=gpar(fontsize=7))
                  )
)
do.call(grid.arrange,args.list)
dev.off()


tiff("S2.tif", height = 5, width = 5, units = 'in', compression="lzw", res=500)
amPredfrommod <- sapply(seq(15,34,0.01),function(x){amFunc(
  mu.a.a1=exp(MLEfits.am)[2]
  , mu.a.a2=exp(MLEfits.am)[3]
  ,x)
})
amPred$amPredmod <- amPredfrommod

ggplot(amDat, aes(x=mean_temp,y=mortality)) +  # set up the plot
  xlim(low=15, high=35) +
  geom_point(size=1) +
  labs( y= "Adult mortality rate (per day)"
        , x=bquote("Temperature"~"(°C)")
        ,title="Black line - fitted to mark-recapture data, grey line - parameter estimates from population model fit") + 
  theme_set(theme_bw()) +
  theme( panel.border = element_blank()
         ,axis.line = element_line(color = 'black')
         ,text=element_text(size=6)
         ,plot.margin=unit(c(0.4,0.2,0.2,0.1), "cm")
         ,axis.text=element_text(size=6)) +
  geom_line(data=amPred
            ,mapping=aes(x=temp,y=amPred)) +
  geom_line(data=amPred
            ,mapping=aes(x=temp,y=amPredfrommod),col="grey")

dev.off()

library("tiff")
S1 <- readTIFF("S1.tiff")
writeTIFF(S1,"S1lzw.tiff",compression="LZW")


# temps <- rep(25,length(temps))
# sim.const.temp <- simPop(parms=tsetse_params(
#   mud.p=exp(MLEfits.am)[1]
#   , mu.a.a1=exp(MLEfits.am)[2]
#   , mu.a.a2=exp(MLEfits.am)[3]
#   
# ))
# plot(sim.const.temp$time,sim.const.temp$A)
