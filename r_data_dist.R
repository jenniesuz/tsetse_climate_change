library(plyr)
raw_dat <- read.csv("data_daily_mean.csv",header=T)
raw_dat <- ddply(raw_dat,.(Year,Month),transform,num=length(Avg))

raw_dat <- raw_dat[raw_dat$num>10,]
names(raw_dat)

years <- unique(raw_dat$Year)
months <- unique(raw_dat$Month)

means.vars <- lapply(unique(raw_dat$Year),function(x){
  dat <- raw_dat[raw_dat$Year %in% x,]
  
  mean.mon <- sapply(unique(dat$Month),function(y){
    sub.dat <- dat[dat$Month %in% y,]
    mean <- mean(sub.dat$Avg)
    return(mean)
  })
  
  var.mon <- sapply(unique(dat$Month),function(y){
    sub.dat <- dat[dat$Month %in% y,]
    var <- var(sub.dat$Avg)
    return(var)
  })
  return(var.mon/mean.mon)
})

varmean<-unlist(means.vars)
hist(varmean)
length(varmean[varmean<1.5])
length(varmean)
91/113
