# This file fits pupal temperature-dependent mortality function using laboratory data available in:
# Phelps 1973 The effect of temperature on fat consumption during the puparial stages of
# Glossina morsitans morsitans Westw. (Dipt., Clossinidae) under laboratory conditions,
# and its implication in the field. Bull. Ent. Res. 62, 423-438
# source("r_2_subfuncs_pupal_duration.R")          # need pupal duration as a function of temperature

#Data for Glossina morsitans
#**************************Data*************************************
pmDat <- read.csv("data_pupal_mort.csv")
pmEmerged <- pmDat$m.emerge + pmDat$f.emerge                                 # total number that emerged
pmDied <- round(pmDat$mortality/100*pmEmerged / (1-pmDat$mortality/100),0)   # estimate number that died
pmTotal <- pmEmerged + pmDied                                                # total
#****Convert proportion dead to instantaneous daily mortality********
pupDur <- 1/pdFunc(temp=pmDat$temp)                  # see r_2_subfuncs_pupal_duration
pmiDat <- -log(1-(pmDied/pmTotal)) / pupDur          # instantaneous mortality over pupal period/ pupal duration at each temperature
pmDat$pmi <- pmiDat
#************************Plot data*********************************
pmPlot <- ggplot(pmDat, aes(x=temp,y=pmi)) +
  xlim(low=15, high=33) +
  geom_point(size=1) +
  ylim(low=0, high=0.04) +
  labs( y="Pupal mortality rate (per day)"
        , x=" " 
        ,title="b") + 
  theme_set(theme_bw()) +
  theme( panel.border = element_blank()
         ,axis.line = element_line(color = 'black')
         ,text=element_text(size=6)
         ,plot.margin=unit(c(0.4,0.2,0.2,0.1), "cm")
         ,axis.text=element_text(size=6)
  )

pmPlot
#**************Mortality function*******************************
pmFunc <- function(mu.p.b1,mu.p.b2,mu.p.b3,mu.p.b4,mu.p.b5,T1=16,T2=32,temp){                                 # This function has 4 parameters to be fitted - k1 to k4
  mu.it <- mu.p.b1 + mu.p.b2*exp(-(mu.p.b3*(temp-T1)))+mu.p.b4*exp(mu.p.b5*(temp-T2))
  return(mu.it)
}
#*********Fit function using nonlinear least squares regression************
pmFit <- nls(pmi ~ mu.p.b1 + mu.p.b2*exp(-(mu.p.b3*(temp-16)))+mu.p.b4*exp(mu.p.b5*(temp-32))
             ,data=pmDat
             ,start=list(mu.p.b1=0.002,mu.p.b2=0.005,mu.p.b3=1.5,mu.p.b4=0.03,mu.p.b5=1.2)
             ,trace=T)
#*****************Fitted parameter values**********************************
pmFitsum <- summary(pmFit)
pmFitsum
b1 <- coef(pmFit)[1]
b2 <- coef(pmFit)[2]
b3 <- coef(pmFit)[3]
b4 <- coef(pmFit)[4]
b5 <- coef(pmFit)[5]
#**************Predicted mortality***************************************
pmPred <- pmFunc(mu.p.b1=b1
                 ,mu.p.b2=b2
                 ,mu.p.b3=b3
                 ,mu.p.b4=b4
                 ,mu.p.b5=b5
                 ,temp=seq(15,34,0.1))
pmPred <- cbind.data.frame(temp = seq(15,34,0.1), pmPred = pmPred)
#***********************Plot********************************************
pmPlotFit <- pmPlot +
  geom_line(data=pmPred
            ,mapping=aes(x=temp,y=pmPred))
pmPlotFit
#***********************************************************************