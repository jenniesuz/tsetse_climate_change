# Pupal duration as a function of temperature using data and fit provided in 
# Phelps (1969) Puparial duration in Glossina morsitans orientalis under
# conditions of constant temperature. Entomologia Experimentalis et Applicata. Vol 12, 33-43

# Data for G. morsitans
#***********Rate of pupal development as a function of temperature******************
pdFunc <- function(p.dur.c1=0.05884
       ,p.dur.c2=4.8829
       ,p.dur.c3=-0.2159
       ,temp){
  pupDur <- (p.dur.c1/(1 + exp(p.dur.c2 + p.dur.c3*temp)))
  return(pupDur)
}

c1=0.05884                    # parameter estimates from reference at top of script
c2=4.8829
c3=-0.2159
#******************Read in data for plotting******************
pdDat <- read.csv("data_pupal_duration.csv")
#*****************Plot****************************************
pdPlot <- ggplot(pdDat, aes(x=test_temp,y=females_r)) +
  xlim(low=15, high=33) +
  geom_point(size=1) +
  ylim(low=0, high=0.06) +
  labs( y="Pupal emergence rate (per day)"
        , x=" "
        ,title="c") + 
  theme_set(theme_bw()) +
  theme( panel.border = element_blank()
         ,axis.line = element_line(color = 'black')
         ,text=element_text(size=6)
         ,plot.margin=unit(c(0.4,0.2,0.2,0.1), "cm")
         ,axis.text=element_text(size=6)
  )

pdPlot
#**************Predicted development rate************************
pdPred <- pdFunc(temp=seq(15,34,0.1))
pdPred <- cbind.data.frame(temp = seq(15,34,0.1), pdPred = pdPred)
#***********************Plot predicted values*********************
pdPlotFit <- pdPlot +
  geom_line(data=pdPred
            ,mapping=aes(x=temp,y=pdPred))
pdPlotFit
#****************************************************************