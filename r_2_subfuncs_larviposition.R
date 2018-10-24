# Fit larviposition rate as a function of temperature
# Data from Hargrove (1994) Reproductive rates of tsetse flies
# in the field in Zimbabwe. Physioloigical Entomology. Vol. 19 (4), 307-318
# and Hargrove (2004) Tsetse Population Dynamics. In: The Trypanosomiases.

lrDat <- read.csv("data_larviposition.csv")

#****************Plot data**************************************
lrPlot <- ggplot(lrDat, aes(x=mean_temp,y=rate,col=stage)) +
  scale_color_manual(values=c("black","grey")) +
  xlim(low=15, high=33) +
  geom_point(size=1) +
  ylim(low=0, high=0.14) +
  labs( y="Larviposition rate (per day)"
        , x=" "
        ,title="d") +
  theme_set(theme_bw()) +
  theme( panel.border = element_blank()
         ,axis.line = element_line(color = 'black')
         ,text=element_text(size=6)
         ,plot.margin=unit(c(0.4,0.2,0.2,0.1), "cm")
         ,axis.text=element_text(size=6)
         ,legend.position =c(0.83,0.2)
         ,legend.title = element_blank()
         ,legend.background = element_blank()
         ,legend.key.size = unit(0.2,"line")
         ,legend.text=element_text(size=4.5)
         ,legend.key.height=unit(0.4,"line")
  )
lrPlot

#******Temperature-dependent function*******************
lrFunc <- function(lr.d1,lr.d2,T1=24,temp){                                 # This function has 4 parameters to be fitted - k1 to k4
  lr <- lr.d1 + lr.d2*(temp-T1)
  return(lr)
}
#*******************************************************************
#******Temperature-dependent function******************************
lrFunc <- function(lr.d1,lr.d2,T1=24,temp){             # two parameters - lr.d1 and lr.d2
  lr <- lr.d1 + lr.d2*(temp-T1)
  return(lr)
}
#**************Parameter estimates********************************
d1=0.1046     # these parameter estimates are taken from the references at the top of this script
d2=0.0052     # store for use in population model
d3=0.061
d4=0.002
#**************Predicted mortality***********************************
lrPred1 <- lrFunc(lr.d1=d1
                  ,lr.d2=d2
                  ,temp=seq(15,34,0.1))
lrPred1 <- cbind.data.frame(temp = seq(15,34,0.1), pred = lrPred1, stage=rep("Subsequent larvipositions",length(lrPred1)))
lrPred2 <- lrFunc(lr.d1=d3
                  ,lr.d2=d4
                  ,temp=seq(15,34,0.1))
lrPred2 <- cbind.data.frame(temp = seq(15,34,0.1), pred = lrPred2,stage=rep("First larviposition",length(lrPred2)))
lrPred <- rbind.data.frame(lrPred1,lrPred2)
#******************Plot fit*****************************************
lrPlotFit <- lrPlot +
  geom_line(data=lrPred
            ,mapping=aes(x=temp,y=pred))
lrPlotFit
#********************************************************************
