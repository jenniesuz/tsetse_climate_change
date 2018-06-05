# This script contains the function defining temperature-dependent
# adult mortality and the model fitted to mark-release-recapture estimates of mortality which can be found in 
# Hargrove (2004) Tsetse population dynamics. In: The Trypanosomiases. Maudlin, I, Holmes, P and Miles, M
# It also produces Figure 3 using the fitted parameters from r_4 the tsetse population dynamics model

require(bbmle)    # required packages
require(ggplot2)

#**************Temperature-dependent adult mortality data*****************
dat <- read.csv("data_adult_mort.csv",header=T)         # read in the data
dat <- dat[!dat$mort %in% NA,]                          # get rid of NA values

# set up the plot
p <- ggplot(dat, aes(x=temp,y=mort)) +
  xlim(low=15, high=35) +
  geom_point(size=2) +
 # ylim(low=0, high=3) +
  labs( y= "Mortality rate"
       , x=bquote("Temperature"~"(°C)")) + 
  theme_set(theme_bw()) +
  theme( panel.border = element_blank()
         ,axis.line = element_line(color = 'black')
         ,text=element_text(size=12)
         ,plot.margin=unit(c(0.2,0.2,0.2,0.1), "cm")
         ,axis.text=element_text(size=12)
  )
#************************************************************************

#*************Adult temperature-dependent mortality function**************
adultmort <- function(mu.a.k1,mu.a.k2,temp){      
if (temp > 25){                                   # if the temperature is >25
  mu.it <- mu.a.k1*exp(mu.a.k2*(temp-25))         # temperature-dependent mortality
}else{                         
  mu.it <- mu.a.k1                               # otherwise no additional temperature-dependent mortality
}
  return(mu.it)
}
#***********************************************************************

#*************Likelihood function**************************************
nll <- function(log_mu.a.k1
                ,log_mu.a.k2
                ,dat=dat){
  mu.a.k1 <- exp(log_mu.a.k1)                                # Parameters to be fitted
  mu.a.k2 <- exp(log_mu.a.k2)
  
  
  mod.mort <- sapply(dat$temp,function(x){ adultmort(        # Model values
     mu.a.k1=mu.a.k1
    ,mu.a.k2=mu.a.k2
    ,x)
  })
  
  obs.mort <- dat$mort                                        # Observed data - fit the logged values
  ll <- sum(dnorm(log(mod.mort),mean=log(obs.mort),log=T))    # Likelihood assuming normal distribution
  return(-ll)
}
#***********************************************************

#************initial parameters***************************
init.pars <- c(log_mu.a.k1=log(0.03)
               ,log_mu.a.k2=log(0.06))
#************************************************************

#***********Fitting******************
objFXN <- function(fit.params
                   ,dat=dat) {
  parms <- fit.params
  nll(fit.params[1],fit.params[2], dat = dat) ## then call likelihood
}
trace <- 3

optim.vals <- optim(par = init.pars           # first stochastic
                    , objFXN
                    , dat = dat
                    , control = list(trace = trace, maxit = 200)
                    , method = "SANN")
exp(optim.vals$par) #                         # second stochastic
optim.vals <- optim(par = optim.vals$par
                    , objFXN
                    , dat = dat
                    , control = list(trace = trace, maxit = 200)
                    , method = "SANN")
exp(optim.vals$par) #


optim.vals <- optim(par = optim.vals$par      # 3rd nelder-mead
                    , objFXN
                    , dat = dat
                    , control = list(trace = trace, maxit = 1000, reltol = 10^-7)
                    , method = "Nelder-Mead" #
                    , hessian = T)
optim.vals # convergence 0 means algorithm converged

MLEfits <- optim.vals$par
exp(MLEfits)

mu.a.k1.fit <- exp(MLEfits)[1]
mu.a.k2.fit <- exp(MLEfits)[2]
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
#********************************************************************


#********************Fitted values as a function of temperature***************
mod.mort <- sapply(seq(15,35,0.01),function(x){ adultmort(             # calcluate model values
  mu.a.k1=mu.a.k1.fit
  ,mu.a.k2=mu.a.k2.fit
  ,x)
})
modmort <- cbind.data.frame(temp = seq(15,35,0.01), mod.mort = mod.mort)
#*************************************************************************

#********************Produce plot*************************
p + 
  geom_line(data=modmort
             ,mapping=aes(x=temp,y=mod.mort))

#******************plot*****************
#tiff("Fig_3_adultmort.tiff", height = 3, width =5, units = 'in', compression="lzw", res=400)
#dev.off()
#***************************************


