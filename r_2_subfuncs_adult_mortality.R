# Fit temperature-dependent adult mortality using mark-release-recapture estimates of mortality from
# Hargrove (2001) Factors affecting density-independent survival of an island population of tsetse flies in Zimbabwe. Entomologia Experimentalis et Applicata. 100(2), 151-164
# Hargrove (2004) Tsetse population dynamics. In: The Trypanosomiases. Maudlin, I, Holmes, P and Miles, M. Cabi Publishing
# Data are for G. pallidipes

require(minpack.lm) # required packages
require(ggplot2) 

#**************Data***************************************************
amDat <- read.csv("data_adult_mort.csv",header=T)         # read in data

amPlot <- ggplot(amDat, aes(x=mean_temp,y=mortality)) +   # set up plot
  xlim(low=15, high=35) +
  geom_point(size=1) +
  labs( y= "Adult mortality rate (per day)"
       , x=" " 
       ,title="a") + 
  theme_set(theme_bw()) +
  theme( panel.border = element_blank()
         ,axis.line = element_line(color = 'black')
         ,text=element_text(size=6)
         ,plot.margin=unit(c(0.4,0.2,0.2,0.1), "cm")
         ,axis.text=element_text(size=6)
  )

amPlot                                                  # view plot
#*************Temperature-dependent function************************
amFunc <- function(mu.a.a1,mu.a.a2,temp){      
  if (temp > 25){                                   # if the temperature is >25
    mu.it <- mu.a.a1*exp(mu.a.a2*(temp-25))         # temperature-dependent mortality
  }else{                         
    mu.it <- mu.a.a1                                # otherwise no additional mortality
  }
 return(mu.it)
}
#*********Fit function using nonlinear least squares regression************
amFuncFit <- function(params,dat=amDat){            # function to calcluate residuals
    a1 <- params[1]
    a2 <- params[2]
    temp <- dat$mean_temp
    mod.mort <- sapply(temp,function(x){ amFunc(    
      mu.a.a1=a1
      ,mu.a.a2=a2
      ,x)
    })
    res <- log(mod.mort) - log(dat$mortality)
    return(res)
}
amFit <- nls.lm(par=c(a1=0.01,a2=0.06),fn=amFuncFit)
#*****************Fitted parameter values******************************
amFitsum <- summary(amFit)                          # store parameter estimates for use in population model
a1 <- coef(amFit)[1]
a2 <- coef(amFit)[2]
#***************Predicted values of mortality*************************
amPred <- sapply(seq(15,34,0.01),function(x){ amFunc(            
   mu.a.a1=a1
  ,mu.a.a2=a2
  ,x)
})
amPred <- cbind.data.frame(temp = seq(15,34,0.01), amPred = amPred)
#********************Add fit to plot***********************************
amPlotFit <- amPlot + 
  geom_line(data=amPred
             ,mapping=aes(x=temp,y=amPred))
amPlotFit
#**********************************************************************