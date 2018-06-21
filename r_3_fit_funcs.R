# This script contains the functions for fitting the model to the
# data and for plotting the results
# Likelihood function takes the model output and the data and calculates the likelihood
# subsParms just enables multiple parameters to be fitted and substituted back into the parameter list



#***********Likelihood function**********************************
nll.pois <- function(parms=tsetse_params(),dat=temps.count){ 
  
 # if (parms$mu.a.a1 > 0.04 | parms$mu.a.a1 < 0.01 ) {                                # penalise model fits where adult mortality at temperatures =< 25 is > 0.03
#    ll <- -1000000000
#  }else{                                                      # if mu.a.k1 parameter is < 0.03 proceed to run the model
    
    simulate <- simPop(parms=parms)                           # run the model and save output to object 'simulate'
    
      counts <- dat$count[1:length(dat$temp)]               # observed counts from bioassay catches
  
      times <- seq(from = 0, to = length(counts)*30-1
                  ,by = 30)                                   # select only the model times that we have data for
  
      times <- times[!counts %in% NA]                         # times have data for
  
      counts <- counts[!counts %in% NA]                       # get rid of NA rows in observed data

      mean.tsetse <- round(counts,0)                         # round observed counts
  
      est <- round(simulate[simulate$time %in% times,"A"]+simulate[simulate$time %in% times,"J"],0)  # model output numbers of adult females at timepoints have data for
  
      subset <- simulate[124:363,"A"] + simulate[124:363,"J"]   # take simulated counts between ~ 1970 and 1990
        
      mean.count <- mean(subset)                              # average count for this period
         
   #    if(mean.count < 50)                                  # if the average is below 50 then penalise model
  #        { ll <- -10000000
    #     }else{                                               # otherwise go ahead and calculate log likelihood

            ll <- sum(dpois(mean.tsetse,lambda=est,log=T))  # log likelihood assuming data are Poisson distributed
     #    }
#  }
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
                   , dat=temps.count) {
  parms <- subsParms(fit.params, fixed.params)
  nll.pois(parms, dat = dat)                                                           ## then call likelihood function
}
#***************************************************************************************

