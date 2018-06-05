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

b.temps$Date <- as.Date(b.temps$Date,"%d/%m/%Y")          # tell R this is a date

#***********add column of year-month format*******
b.counts$Mon.yr <- as.yearmon(b.counts$Mon.yr,"%d/%m/%Y") 
b.temps$yearmon <- as.yearmon(b.temps$Date,"%Y/%m/%d")    
b.temps$times <- c(1:length(b.temps$MeanC))
#***************************************

b.temps <- b.temps[b.temps$yearmon %in% b.counts$Mon.yr,]

count.times <- ddply(b.temps,.(yearmon),summarise,median.times=round(median(times),0))   # find time point in middle of each month
count.times <- count.times[count.times$yearmon %in% b.counts$Mon.yr,]
count.times$counts <- b.counts$Mean
count.times <- count.times[!count.times$counts %in% NA,]
count.times
#***************************************************************