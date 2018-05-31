
# Plot Phelps and Burrows (1969) Puparial duration in Glossina morsitans 
# orientalis under conditions of constant temperature

# Add Phelps and Hargrove model fits

phelps.dat <- read.csv("data_pupal_duration_phelps1969.csv")

# phelps fit
k <- 0.05884
a <- 4.8829
b <- -0.2159

# hargs fit
hk <- 0.057
ha <- 5.5
hb <- -0.25

par(mfcol=c(1,2))
plot(phelps.dat$test_temp
     , phelps.dat$females_d
     ,bty="n"
     ,xlab="Temperature"
     ,ylab="Pupal duration")

curve(1/(k/(1+exp(a+b*x)))
      ,add=T
      ,col="red")
curve(1/(hk/(1+exp(ha+hb*x)))
      ,add=T
      ,col="blue")


plot(phelps.dat$test_temp
     , phelps.dat$females_r
     ,bty="n"
     ,xlab="Temperature"
     ,ylab="Pupal development rate")

curve(k/(1+exp(a+b*x))
      ,add=T
      ,col="red")
curve(hk/(1+exp(ha+hb*x))
      ,add=T
      ,col="blue")
legend("topleft",legend=c("Hargrove fit", "Phelps fit")
       ,bty="n"
       ,col=c("blue","red")
       ,lty=1
)

