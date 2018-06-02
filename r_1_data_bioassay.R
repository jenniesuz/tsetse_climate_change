# This script:
# 1 reads in the tsetse count and temperature data
# 2 formats the time columns and summarises the temperature data
# 3 binds the count data and temperature data together 

library("zoo")     # required packages
library("plyr")



#***********read in data**********************************
b.counts <- read.csv("data_bioassay_counts.csv",header=T) # bioassay count data

b.temps <- read.csv("data_bioassay_temps.csv",header=T)   # temperature data
#********************************************************

#***********change time column to year-month format*******
b.temps$Date <- as.Date(b.temps$Date,"%d/%m/%Y")          # tell R this is a date


#b.counts$Mon.yr <- as.yearmon(b.counts$Mon.yr,"%d/%m/%Y") # format count data date column to year month 

#b.temps$yearmon <- as.yearmon(b.temps$Date,"%Y/%m/%d")    # format temp data date column to year month

#twenty16 <- b.temps[grep("2016",b.temps$yearmon),]
#plot(twenty16$Date,twenty16$MeanC,bty="n",xlab="Date",ylab="Mean Temp.")
#*********************************************************

#**********summarise to monthly means for modelling**************************
b.temp <- ddply(b.temps,.(yearmon),summarise,mean.temp=mean(MeanC,na.rm=T))  # for each level in year-month column calculate the mean temperature

names(b.temp) <- c("time","temp")                                            # rename the columns in the table
#***************************************************************************



#**********bind data together into a single table*************************************
temps.count <- cbind.data.frame(time=c(b.temp$time)
                                ,temp=c(b.temp$temp)
                                ,count=c(rep(NA,369),b.counts$Mean[1:(length(b.counts$Mean)-1)])
                                ,count.time=as.yearmon(c(rep(NA,369),b.counts$Mon.yr[1:(length(b.counts$Mean)-1)]))
                                )

#extend beginning by 5 years to allow population to stabilise - use a repeat of the first year
temps.count <- rbind.data.frame (temps.count[13:24,]
                                ,temps.count[13:24,]
                                ,temps.count[13:24,]
                                ,temps.count[13:24,]
                                ,temps.count[13:24,]
 ,temps.count
)

#**************************************************************************************


#*********************add in additional data***CURRENTLY OMITTED*************************
# # Pilson data
# temps.count$count[2] <- 101
# temps.count$count[5] <- 22
# temps.count$count[8] <- 96
# temps.count$count[11] <- 60
# # additional from GV
# temps.count$count[85] <- 49
# temps.count$count[86] <- 26
# temps.count$count[97] <- 112
# 
# temps.count$count[140] <- 34
# temps.count$count[142] <- 78
# temps.count$count[143] <- 125
# 
# temps.count$count[205] <- 28
# 
# temps.count$count[284] <- 78
# temps.count$count[285] <- 93
# temps.count$count[286] <- 88
# temps.count$count[287] <- 82
# temps.count$count[288] <- 54
# temps.count$count[289] <- 36
# temps.count$count[290] <- 31
# 
# temps.count$count[292] <- 16
# temps.count$count[293] <- 20
# temps.count$count[294] <- 33
# temps.count$count[295] <- 14
# temps.count$count[296] <- 10
# temps.count$count[297] <- 11
# temps.count$count[298] <- 32

