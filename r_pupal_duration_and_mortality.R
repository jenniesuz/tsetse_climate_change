# Set system run time
# sim.start.time <- Sys.time()

# Set working directory
#setwd(getwd())

# === Load all necessary packages
library(tidyverse)

# Read in files
T_bar <- readxl::read_xlsx("data_tbar_92_93_JH.xlsx", col_names = TRUE)
colnames(T_bar) <- c("Date","TempScreen")

# Calculate Temperature at Site of Pupal deposit
# linear y = mx + c
m <- 0.88
c <- 1.2



# view T_bar
# head(T_bar,5)
# attributes(T_bar)

# Formula to calculate the rate at which pupal develops depending on the temperature of the day
# Rate formula 
# r(T) = k/(1+exp(a + bT))
# fixed variables

# Mortality
# Full model Phelps 
k0 =	0.00202
k1 =	0.0053
k2 =	1.5524
k3 =	0.0297
k4 =	1.2175
T1 =	16
T2 =	32 

# Additional columns
T_bar$PupalMortalityDailyScreen <- NA
T_bar$PupalDurationScreen <- NA
T_bar$MeanTemperatureScreen <- NA
T_bar$InstantPupalMortalityScreen <- NA
T_bar$PupalSurvivalScreen <- NA
T_bar$TempSite <- NA
T_bar$PupalDurationSite <- NA
T_bar$MeanTemperatureSite <- NA
T_bar$InstantPupalMortalitySite <- NA
T_bar$PupalSurvivalSite <- NA
T_bar$Simulation <- NA

Temp.Duration.Function <- function(temp.data, date.data, row.position, # data
                                   k = 0.057, a = 5.5, b = -0.25) # rate equation coefficients
{
  pupae.duration = 0 
  mean.temp = 0
  counter = 0
  temp.total = 0
  start.date <- date.data[row.position]

  PupalDuration <- NA
  MeanTemperature <- NA

  for (j in row.position:length(temp.data))
  {
    temp <- temp.data[j]
    
    if(pupae.duration < 1) 
    {
      temp.total <- temp.total + temp
      counter <- counter + 1
      rate_T <- k / (1 + exp( a + b * temp)) # Determine r(T) for each temperature of the day

      pupae.duration <- pupae.duration + rate_T  
    }else{
      temp.total <- temp.total - temp.data[j-1]
      # calculate the total number of days for the pupal duration
      PupalDuration <- date.data[j] - start.date - 1
      # calculate the mean temperature over the pupal duration
      MeanTemperature <- temp.total/(counter - 1)
      break
    }
  }
  return(c(PupalDuration,MeanTemperature))
}


# Single run
T_bar$TempSite <- (m * T_bar$TempScreen) + c
for (i in 1:nrow(T_bar))
{

  function.output.screen <-  Temp.Duration.Function(T_bar$TempScreen, T_bar$Date, i)
  T_bar$PupalMortalityDailyScreen[i] <- k0 + (k1 * exp(-k2*(T_bar$TempScreen[i] - T1))) + (k3 * exp(k4 * (T_bar$TempScreen[i] - T2)))
  T_bar$PupalDurationScreen[i] <- function.output.screen[1]
  T_bar$MeanTemperatureScreen[i] <- function.output.screen[2]
  T_bar$InstantPupalMortalityScreen[i] <- k0 + (k1 * exp(-k2*(T_bar$MeanTemperatureScreen[i] - T1))) + (k3 * exp(k4 * (T_bar$MeanTemperatureScreen[i] - T2)))
  T_bar$PupalSurvivalScreen[i] <- exp(-T_bar$PupalDurationScreen[i] * T_bar$InstantPupalMortalityScreen[i])*100
  
  function.output.site <- Temp.Duration.Function(T_bar$TempSite, T_bar$Date, i)
  T_bar$PupalDurationSite[i] <- function.output.site[1]
  T_bar$MeanTemperatureSite[i] <- function.output.site[2]
  T_bar$InstantPupalMortalitySite[i] <- k0 + (k1 * exp(-k2*(  T_bar$MeanTemperatureSite[i] - T1))) + (k3 * exp(k4 * (T_bar$MeanTemperatureSite[i] - T2)))
  T_bar$PupalSurvivalSite[i] <- exp(-T_bar$PupalDurationSite[i] * T_bar$InstantPupalMortalitySite[i])*100
}

summary(T_bar)
# Graphs - add legends, titles, axis limits, etc
plot(T_bar$Date,T_bar$PupalDurationScreen, type = "l" )
lines(T_bar$Date, T_bar$TempScreen, col = "blue")
lines(T_bar$Date, T_bar$MeanTemperatureScreen, col = "red", lwd = 2)
# Graphs - add legends, titles, axis limits, etc
plot(T_bar$Date, T_bar$MeanTemperatureScreen, col = "blue", type = "l", ylim = c(10,50))
lines(T_bar$Date, T_bar$MeanTemperatureSite, col = "red", lwd = 2)
lines(T_bar$Date, T_bar$PupalDurationScreen, col = "blue", lwd = 2, lty = 2)
lines(T_bar$Date, T_bar$PupalDurationSite, col = "red", lwd = 2, lty = 2)

write.csv(T_bar, "Pupal_Duration.csv")






#==== ignore below === work in progress =======#
# For a multiple run with an error term added on to the TempSite (temperature at site of pupal deposit)
# To calculate the mean , standard deviation , standard error and CIs

for (j in 1:1000)
{
  print(Sys.time())
  T_bar$TempSite <- (m * T_bar$TempScreen) + c + rnorm(1,mean=0,sd=1)
  T_bar$Simulation <- paste0("Sim_",j)
  for (i in 1:nrow(T_bar))
  {
    function.output.screen <-  Temp.Duration.Function(T_bar$TempScreen, T_bar$Date, i)
    T_bar$PupalMortalityDailyScreen[i] <- k0 + (k1 * exp(-k2*(T_bar$TempScreen[i] - T1))) + (k3 * exp(k4 * (T_bar$TempScreen[i] - T2)))
    T_bar$PupalDurationScreen[i] <- function.output.screen[1]
    T_bar$MeanTemperatureScreen[i] <- function.output.screen[2]
    T_bar$InstantPupalMortalityScreen[i] <- k0 + (k1 * exp(-k2*(T_bar$MeanTemperatureScreen[i] - T1))) + (k3 * exp(k4 * (T_bar$MeanTemperatureScreen[i] - T2)))
    T_bar$PupalSurvivalScreen[i] <- exp(-T_bar$PupalDurationScreen[i] * T_bar$InstantPupalMortalityScreen[i])*100
    
    function.output.site <- Temp.Duration.Function(T_bar$TempSite, T_bar$Date, i)
    T_bar$PupalDurationSite[i] <- function.output.site[1]
    T_bar$MeanTemperatureSite[i] <- function.output.site[2]
    T_bar$InstantPupalMortalitySite[i] <- k0 + (k1 * exp(-k2*(  T_bar$MeanTemperatureSite[i] - T1))) + (k3 * exp(k4 * (T_bar$MeanTemperatureSite[i] - T2)))
    T_bar$PupalSurvivalSite[i] <- exp(-T_bar$PupalDurationSite[i] * T_bar$InstantPupalMortalitySite[i])*100
  }
  Sims <- bind_rows(Sims, T_bar)
}

# Using Sims, we can now calculate the mean, se and CIs
# Incomplete - Work in progress
            # Mean.df <- data.frame(matrix(NA, nrow = nrow(T_bar),ncol=13))
            # colnames(Mean.df) <- colnames(T_bar)
            # Mean.df$Date <- T_bar$Date
            # for (z in 2:(ncol(T_bar)-1))
            # {
            #   aa  <- Sims %>%
            #     group_by(Date) %>% print(Sims[[2]])
            #     summarise(avg = mean(Sims[[2]])) %>%
            #     arrange(Date)
            #   Mean.df[z] <- aa$avg
            # }
