# This script contains the functions for fitting the model to the
# data and for plotting the results
# Likelihood function takes the model output and the data and calculates the likelihood
# subsParms just enables multiple parameters to be fitted and substituted back into the parameter list

#***********Likelihood function**********************************
nll.pois <- function(parms=tsetse_params(),dat=count.times,temps.yearmon=b.temps$yearmon){ 
  
  if (parms$mu.a.k1 > 0.04 | parms$mu.a.k1 < 0.01 ) {         # penalise model fits where adult mortality at temperatures =< 25 is > 0.03
    ll <- -1000000000
  }else{                                                      # if mu.a.k1 parameter is < 0.03 proceed to run the model
    
    simulate <- simPop(parms=parms)                           # run the model and save output to object 'simulate'

      mean.tsetse <- round(dat$counts,0)     # round observed counts
      
      sim.date <- cbind.data.frame(A=simulate$A,yearmon=temps.yearmon)
      est <- ddply(sim.date,.(yearmon),summarise,mean.sim=mean(A))  # take mean numbers for the month
      est <- round(est$mean.sim[est$yearmon %in% dat$yearmon],0)
  
      subset <- sim.date$A[(sim.date$yearmon > as.yearmon("August 1970"))& (sim.date$yearmon < as.yearmon("February 1990"))]                       # take simulated counts between ~ 1970 and 1990
        
     # mean.count <- mean(subset)                            # average count for this period
         
      # if(mean.count < 50)                                  # if the average is below 50 then penalise model
      #     { ll <- -10000000
      #   }else{                                             # otherwise go ahead and calculate log likelihood

            ll <- sum(dpois(mean.tsetse,lambda=est,log=T))  # log likelihood assuming data are Poisson distributed
       #  }
  }
  return(-ll)
}
#**********************************************************



#********enable multiple parameters to be fit and substituted back in parameter list****
subsParms <- function(fit.params, fixed.params=tsetse_params())
  
  within(fixed.params, {
    loggedParms <- names(fit.params)[grepl('log_', names(fit.params))]
    
    unloggedParms <- names(fit.params)[!grepl('log_', names(fit.params))]  
    
    for(nm in unloggedParms) assign(nm, as.numeric(fit.params[nm]))
    
    for(nm in loggedParms) assign(gsub('log_','',nm), exp(as.numeric(fit.params[nm])))
    
    rm(nm, loggedParms, unloggedParms)
  })           

                                                                                       ## Make likelihood a function of fixed and fitted parameters.
objFXN <- function(fit.params                                                          ## paramters to fit
                   , fixed.params =tsetse_params()                                     ## fixed paramters
                   , dat=count.times) {
  parms <- subsParms(fit.params, fixed.params)
  nll.pois(parms, dat = dat)                                                           ## then call likelihood function
}
#***************************************************************************************

