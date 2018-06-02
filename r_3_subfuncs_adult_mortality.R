# This script contains the function defining temperature-dependent
# adult mortality and the model fitted to mark-release-recapture estimates of mortality
# It also produces Figure 3 using the fitted parameters from r_4
require(bbmle)


dat <- read.csv("data_adult_mort.csv",header=T)


#**************Temperature-dependent adult mortality*****************
temp <- seq(15,35,0.01)                                  # vector of temperatures
                           

mu.at <- exp(-1.47+0.106*temp)/100                       # as per mark-recapture data analysis

#*******************THESE ARE FITTED VALUES FROM r_4_tsetsemode_adultmortpupalmort.R*************
mu.a.k1 <- 3.808201e-02                                  # fitted values from ODE model
mu.a.k2 <- 8.417606e-02
#***********************************************************************************************

adultmort <- function(temp){                             # calculate from populaiton model fitted values for each temperature
if (temp > 25){                                          # if the temperature is >25
  mu.it <- mu.a.k1*exp(mu.a.k2*(temp-25))                # temperature-dependent mortality
}else{                         
  mu.it <- mu.a.k1                                       # otherwise no temperature-dependent mortality
}
  return(mu.it)
}

mu.it <- sapply(temp,adultmort)

adultmortjs <- function(temp){                            # estimates from mark-recapture data
  if (temp > 25){
mu.at <- exp(-1.47+0.106*temp)/100
  }else{
  mu.at <- 0.023
}
}

mu.at <- sapply(temp,adultmortjs)


#******************plot*****************
tiff("Fig_3_adultmort.tiff", height = 3, width =5, units = 'in', compression="lzw", res=400)
par(mar=c(3,4,1,1),mgp=c(2,1,0),cex=0.5)


#***********Plot data******************

plot(dat$temp
     ,dat$mort
     ,bty="n"
     ,pch=19
     ,col="darkgrey"
     ,ylim=c(0.02,0.12)
     ,xlim=c(15,35)
     ,xlab=expression(paste("Temperature (",degree,"C)"))
     ,ylab=expression(paste("Mortality rate (days"^"-1",")"))
)

par(new=T)
plot(temp,mu.at,bty="n",
     type="l",col="grey",xaxt="n",yaxt="n",ylab=" ",xlab=" ",xlim=c(15,35),ylim=c(0.02,0.12))

legend("topleft",legend=c(
                          "Fitted to mark-recapture data","Fitted to population data"),col=c("grey","black")
       ,lty=c(1,1),bty="n")
dev.off()
#***************************************


