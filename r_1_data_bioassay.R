# 1 reads in the tsetse count and temperature data
# 2 formats the time columns and summarises the temperature data
# 3 binds the count data and temperature data together 

library("zoo")     # required packages
library("plyr")

#***********read in data**********************************
b.counts <- read.csv("data_bioassay_counts.csv",header=T) # bioassay count data

b.temps <- read.csv("data_bioassay_temps.csv",header=T)   # temperature data
#***********change time column to year-month format*******
b.temps$Date <- as.Date(b.temps$Date,"%d/%m/%Y")          # format to date

b.counts$Mon.yr <- as.yearmon(b.counts$Mon.yr,"%d/%m/%Y") # format count data date column to year month 

b.temps$yearmon <- as.yearmon(b.temps$Date,"%Y/%m/%d")    # format temp data date column to year month
#**********summarise to monthly means for modelling**************************
b.temp <- ddply(b.temps,.(yearmon),summarise,mean.temp=mean(MeanC,na.rm=T))  # for each level in year-month column calculate the mean temperature

names(b.temp) <- c("time","temp")                                            # rename the columns in the table
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
#************************************************************************************
