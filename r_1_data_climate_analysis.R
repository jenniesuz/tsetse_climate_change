# This script produces the Figure 1 climate data

library("TTR")                                                # required packages
library("tseries")
library("L1pack")
library("forecast")

source("r_1_data_bioassay.R")                                 # read in climate data

#*****create 30 year monthly mean reference period:***********
ref.period <- b.temp[4:363,]                                   # Jan 1960 to Dec 1989

ref.period$month <- rep(1:12,30)                               # create a new column of just months 1 to 12

ref.means <- ddply(ref.period,.(month)                         # calculate the average for each month over the 30 years
                   ,summarise,refMean=mean(temp,na.rm=T))
#************************************************************

#*******attach the monthly reference means to the dataset and calc anomalies*****
b.temp$ref <- c(ref.means[10:12,2]                            # create a new column containing the reference means for each month
                ,rep(ref.means[,2],57)
                ,ref.means[1:6,2])

b.temp$anom <- b.temp$temp - b.temp$ref                      # for each month take the reference mean from the actual mean

anomCts <- ts(data=b.temp$anom                                # create a time series object from the monthly anomalies
              ,start=b.temp$time[1]                           # period 12 as monthly data
              ,end=b.temp$time[length(b.temp$time)]
              ,frequency=12)


anom.smooth <- SMA(anomCts,60)                                # create a five-year running mean for anomalies

anomCts.lm <- tslm(anomCts~season+trend)
#********************************************************************************

#***********plot the climate data***********************************
tiff("Fig_1_temps.tiff", height = 3.5, width = 5, units = 'in', compression="lzw", res=400)
par(mar=c(4,4.2,1,1),mgp=c(3,1,0),mfcol=c(1,1),cex=0.7)

plot(b.temp$time
     ,anom.smooth
     ,bty="n",type="l"
     ,xlab="Date"
     ,xlim=c(as.yearqtr('1960 Q4'),
            as.yearqtr('2017 Q4'))
     ,ylab=expression(paste("Temperature ( ",degree,"C) anomalies")))


dev.off()
#*******************************************************************


#********************************************************************************************************************


#******************By Month***********************************

monthplot(anomCts)
acf(anomCts)
test <- ets(anomCts, model="AAA")

test <- tslm(anomCts~trend + season)

p <- predict(test,new.data=c(1959:2017),interval="confidence")

plot(b.temp$time,b.temp$anom)
plot(1959:2017,p[,1])

par(mfcol=c(1,1))
monthly.trends <- sapply(month.abb,function(x){
  temporary <- b.temp[grep(x,b.temp$time),]
  temp.ts <- ts(data=temporary$temp
                ,start=temporary$time[1]
                ,end=temporary$time[length(temporary$time)]
                ,frequency=1)
  
  test <- tslm(temp.ts~trend)
  print(summary(test))
  predt <- forecast(test,newdata=c(1959:2017), level=0.95)
  predt <- as.data.frame(predt)
  fitdiff <- predt[57,1] - predt[1,1]
  lwrdiff <- predt[57,2] - predt[1,2]
  upprdiff <- predt[57,3] - predt[1,3]
  return(c(fitdiff,lwrdiff,upprdiff))
  #return(predt)
})

tiff("Fig_3.tiff", height = 3.5, width = 5, units = 'in', compression="lzw", res=400)
par(mar=c(4,4.2,1,1),mgp=c(3,1,0),mfcol=c(1,1),cex=0.7)
b <- barplot(monthly.trends[1,1:12]
        ,names.arg=month.abb
        ,xlab="Month"
        ,ylab=expression(paste("Temperature ( ",degree,"C) increase 1959 - 2017"))
        ,ylim=c(0,2.8))
segments(b,monthly.trends[2,1:12],b,monthly.trends[3,1:12])
text(x=2, y=max(monthly.trends[3,2]), "*", pos=3, cex=1.2)


dev.off()








# #***********plot the climate data***********************************
# tiff("Fig_1_temps.tiff", height = 4, width = 5, units = 'in', compression="lzw", res=400)
# par(mar=c(4,4.2,1,1),mgp=c(3,1,0),mfcol=c(2,1),cex=0.5)
# plot(b.temp$time,b.temp$MeanC
#      ,xlab="Date"
#      ,ylab=expression(paste("Mean monthly temperature ( ",degree,"C)"))
#      ,bty="n"
#      # ,ylim=c(10,40)
#      ,type="l")
# 
# plot(b.temp$time
#      ,anom.smooth
#      ,bty="n",type="l"
#      ,xlab="Date"
#      ,ylab=expression(paste("Temperature ( ",degree,"C) anomalies")))
# 
# 
# dev.off()
# #*******************************************************************











