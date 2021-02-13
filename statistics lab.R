# Program for statistics lab
# Michael Sigler
# May 19, 2019, Revised November 22, 2019, March 24, 2020, February 12, 2021

# Northeast US continental shelf plankton data, NOAA Northeast Fisheries Science Center
# Sampling years
#   1977-1987, MARMAP
#   1988-1998, some sampling occurred, but less balanced seasonally
#   1999-2017, Ecosystem Monitoring (EcoMon)
#   2018-2019, sampling occurred but data currently unavailable publicly
# Focus on 6 fish species (ichthyoplankton) that are primary prey of terns

# Call packages (libraries) used in data analysis; run every time
library(dplyr)
library(ggplot2)
library(tidyr)
library(knitr)
library(lubridate)
library(mgcv)
library(fields)
library(maps)
library(mapdata)
library(mapproj)

############################################################################################################
# Read plankton file
# This file is compiled from the original database available at:
# https://www.ncei.noaa.gov/metadata/geoportal/rest/metadata/item/gov.noaa.nodc%3A0187513/html
# Key variables
# lat = latitude
# lon = longitude
# date = date (GMT)
# sfc_temp = surface temperature
# sfc_salt = surface salinity
# volume_100m3 = Zooplankton Displacement Volume (ml)  per 1m2 of surface area
# zoop_100m3 = Zooplankton count per 100 m3 of volume
# fish_100m3 = Fish count per 100 m3 of volume
# polvir_100m3 = pollock count per 100 m3 of volume
# cluhar_100m3 = Atlantic herring count per 100 m3 of volume
# ammspp_100m3 = sandlance count per 100 m3 of volume
# pepspp_100m3 = butterfish count per 100 m3 of volume
# scosco_100m3 = Atlantic mackerel count per 100 m3 of volume
# urospp_100m3 = red and white hake count per 100 m3 of volume
# e.g., log_cluhar = log 10 transformed (value +1) of cluhar_100m3
# distance = distance (km) from sample location to Appledore Island 

rm(list = ls())   # clear memory

# Appledore Island location
Lat.SML <- 42.987727
Lon.SML <- -70.613940
appledore.location <- c(Lon.SML,Lat.SML)

##################################
#1. Read data and cross tabulate.#
##################################
# set project directories
projdir <- "D:\\Documents\\Shoals\\2021 class\\Analysis\\plankton\\"
projdir.dat <- "D:\\Documents\\Shoals\\2021 class\\Analysis\\data\\"
projdir.results <- "D:\\Documents\\Shoals\\2021 class\\Analysis\\results\\"

# read data
filename <- paste(projdir.dat,"Plankton.csv",sep="") # filename
Catch <- as_tibble(read.csv(filename,header=T))      # read file as a "tibble"
Catch$date <- as.Date(Catch$date, "%Y-%m-%d")		     # Assign Date format 'yyyy-mm-dd'

# Crosstabulate the data by year and month
Catch %>%
  group_by(year, month)%>%                      # group by year and month
  summarise(n=n())%>%                           # summarise by count n=n()
  spread(month, n)%>%                           # spread month as header of table
  kable()                                       # format table

####################################################################################
#2. Fit General Additive Model (GAM) to log-transformed values, log10(value+1), include year and month terms#
# One spatial map for each species covering all years (i.e., average location)     #
#  run time about 1 minute                                                         #
####################################################################################

start.year <- 1977  # all years
Catch.ich <- subset(Catch, year>=start.year)         

# First GAM, fit as a function of location (lon, lat), year and month
# Atlantic herring (Clupea harengus)
# run GAM
fit.cluhar <- gam(log_cluhar ~ s(lon, lat) + s(year) + s(month,bs='cc'), data=Catch.ich) 

# 2a. diagnostics
summary(fit.cluhar)$r.sq                                   # R-squared: proportion of variance explained
hist(residuals(fit.cluhar, type = "deviance"),             # histogram of residuals
     xlab = "Residuals", main = "Histogram of residuals")

# 2b. plot results
plot(fit.cluhar,pages=0,seWithMean=TRUE,select=2,
     ylim=c(-0.3,0.3),main="herring")                      # year plot
plot(fit.cluhar,pages=0,seWithMean=TRUE,select=3,
     ylim=c(-0.3,0.3),main="herring")                      # month plot
vis.gam(fit.cluhar, plot.type="contour", too.far=0.15,     # map
        main="Atlantic herring", xlab="Longitude",ylab="Latitude",color="topo")
   points(Lon.SML, Lat.SML, pch=16, col = "red")
   map('worldHires',fill=T,add=T, col="grey")

# 3. Run GAMs for the remaining species
fit.urospp <- gam(log_urospp ~ s(lon, lat) + s(year) + s(month,bs='cc'), data=Catch.ich)
fit.polvir <- gam(log_polvir ~ s(lon, lat) + s(year) + s(month,bs='cc'), data=Catch.ich)
fit.scosco <- gam(log_scosco ~ s(lon, lat) + s(year) + s(month,bs='cc'), data=Catch.ich)
fit.pepspp <- gam(log_pepspp ~ s(lon, lat) + s(year) + s(month,bs='cc'), data=Catch.ich)
fit.ammspp <- gam(log_ammspp ~ s(lon, lat) + s(year) + s(month,bs='cc'), data=Catch.ich)

# Output all plotted results to a file
# Look at this file after running this code
pdf(paste(projdir.results,"GAM 3 output.pdf",sep=""))
par(mfrow = c(2,3))

#by year
plot(fit.cluhar,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.3,0.3),main="herring")
plot(fit.urospp,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.2,0.2),main="red/white hake")
plot(fit.polvir,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.1,0.1),main="pollock")
plot(fit.scosco,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.1,0.1),main="mackerel")
plot(fit.pepspp,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.1,0.1),main="butterfish")
plot(fit.ammspp,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.4,0.4),main="sandlance")

#by month
plot(fit.cluhar,pages=0,seWithMean=TRUE,select=3,ylim=c(-0.3,0.3),main="herring")
plot(fit.urospp,pages=0,seWithMean=TRUE,select=3,ylim=c(-0.4,0.4),main="red/white hake")
plot(fit.polvir,pages=0,seWithMean=TRUE,select=3,ylim=c(-0.3,0.3),main="pollock")
plot(fit.scosco,pages=0,seWithMean=TRUE,select=3,ylim=c(-0.3,0.3),main="mackerel")
plot(fit.pepspp,pages=0,seWithMean=TRUE,select=3,ylim=c(-0.3,0.3),main="butterfish")
plot(fit.ammspp,pages=0,seWithMean=TRUE,select=3,ylim=c(-0.6,0.6),main="sandlance")

# spatial plots
par(mfrow = c(1,1))
vis.gam(fit.cluhar, plot.type="contour", too.far=0.15, main="Atlantic herring", xlab="Longitude",ylab="Latitude",color="topo")
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
vis.gam(fit.urospp, plot.type="contour", too.far=0.15, main="red white hake", xlab="Longitude",ylab="Latitude",color="topo")
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
vis.gam(fit.polvir, plot.type="contour", too.far=0.15, main="pollock", xlab="Longitude",ylab="Latitude",color="topo")
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
vis.gam(fit.scosco, plot.type="contour", too.far=0.15, main="Atlantic mackerel", xlab="Longitude",ylab="Latitude",color="topo")
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
vis.gam(fit.pepspp, plot.type="contour", too.far=0.15, main="butterfish", xlab="Longitude",ylab="Latitude",color="topo")
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
vis.gam(fit.ammspp, plot.type="contour", too.far=0.15, main="sandlance", xlab="Longitude",ylab="Latitude",color="topo")
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")

dev.off()

############################################################
#4. Fit GAM to log-transformed values, accounting for month#
# Spatial model by year                                    #
#   More detailed model than previous model.
#   This model fits a map for each year 
#     whereas the previous model fitted one model for all years

# read data
filename <- paste(projdir.dat,"Plankton.csv",sep="") # filename
Catch <- as_tibble(read.csv(filename,header=T))      # read file as a "tibble"
Catch$date <- as.Date(Catch$date, "%Y-%m-%d")		     # Assign Date format 'yyyy-mm-dd'

# Choose analysis years
start.year <- 1999                               # Ecomon years: 1999-present
end.year <- 2017                                 # Last year in data set
num.year <- end.year - start.year + 1
Catch.ich <- subset(Catch, year>=start.year) 

# Run year-specific map model for butterfish 
pdf(paste(projdir.results,"GAM 4 butterfish output.pdf",sep=""))
par(mfrow = c(2,3))
  fit.pepspp <- gam(Catch.ich$log_pepspp ~ s(lon, lat, by=as.factor(year),id=0) +
                 s(month,bs='cc'),
                 data=Catch.ich, dist ="normal")   # fit GAM to butterfish data  
  summary(fit.pepspp)$r.sq                                   # R-squared: proportion of variance explained
  plot(fit.pepspp,pages=0,seWithMean=TRUE,select=num.year+1,ylim=c(-0.4,0.4),
       main=paste("butterfish by month"))     # by month plot
  
# by year group plots  
  for (k in start.year:end.year) {
    vis.gam(fit.spp, view=c("lon","lat"), plot.type="contour", too.far=0.15, cond=list(year=k),
            main=paste("butterfish year ", k), xlab="Longitude",ylab="Latitude",color="topo") 
    points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
  }
par(mfrow = c(1,1))
dev.off()

############################################################
#5. Fit GAM to log-transformed values, accounting for month#
# Spatial model by year                                    #
#   More detailed model than previous model.
#   This model fits a map for each year 
#     whereas the previous model fitted one model for all years
### DO NOT RUN ### TAKES SEVERAL MINUTES ON A "FAST" COMPUTER
# ONLY REVIEW THE OUTPUT IN THE FILE "GAM 5 all species output.pdf"

# Run year-specific map model for remaining species 
fit.cluhar <- gam(Catch.ich$log_cluhar ~ s(lon, lat, by=as.factor(year),id=0) +
                    s(month,bs='cc'),
                  data=Catch.ich, dist ="normal")          # fit GAM   
fit.urospp <- gam(Catch.ich$log_urospp ~ s(lon, lat, by=as.factor(year),id=0) +
                    s(month,bs='cc'),
                  data=Catch.ich, dist ="normal")          # fit GAM   
fit.polvir <- gam(Catch.ich$log_polvir ~ s(lon, lat, by=as.factor(year),id=0) +
                    s(month,bs='cc'),
                  data=Catch.ich, dist ="normal")          # fit GAM   
fit.scosco <- gam(Catch.ich$log_scosco ~ s(lon, lat, by=as.factor(year),id=0) +
                    s(month,bs='cc'),
                  data=Catch.ich, dist ="normal")          # fit GAM   
fit.ammspp <- gam(Catch.ich$log_ammspp ~ s(lon, lat, by=as.factor(year),id=0) +
                 s(month,bs='cc'),
               data=Catch.ich, dist ="normal")          # fit GAM   
# by year group plots  
pdf(paste(projdir.results,"GAM 5 all species output.pdf",sep=""))
par(mfrow = c(2,3))
for (k in start.year:end.year) {
  vis.gam(fit.cluhar, view=c("lon","lat"), plot.type="contour", too.far=0.15, cond=list(year=k),
          main=paste("herring year ", k), xlab="Longitude",ylab="Latitude",color="topo") 
  points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
  vis.gam(fit.urospp, view=c("lon","lat"), plot.type="contour", too.far=0.15, cond=list(year=k),
          main=paste("red white hake year ", k), xlab="Longitude",ylab="Latitude",color="topo") 
  points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
  vis.gam(fit.polvir, view=c("lon","lat"), plot.type="contour", too.far=0.15, cond=list(year=k),
          main=paste("pollock year ", k), xlab="Longitude",ylab="Latitude",color="topo") 
  points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
  vis.gam(fit.scosco, view=c("lon","lat"), plot.type="contour", too.far=0.15, cond=list(year=k),
          main=paste("mackerel year ", k), xlab="Longitude",ylab="Latitude",color="topo") 
  points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
  vis.gam(fit.pepspp, view=c("lon","lat"), plot.type="contour", too.far=0.15, cond=list(year=k),
          main=paste("butterfish year ", k), xlab="Longitude",ylab="Latitude",color="topo") 
  points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
  vis.gam(fit.ammspp, view=c("lon","lat"), plot.type="contour", too.far=0.15, cond=list(year=k),
          main=paste("sand lance year ", k), xlab="Longitude",ylab="Latitude",color="topo") 
  points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
}
par(mfrow = c(1,1))
dev.off()

# some extra code
# Table of species names
name.tbl <- as_tibble(data.frame(spp.num = seq(26,31),
                       spp.name = c("herring","red.white.hake","pollock","mackerel","butterfish","sandlance"),
                       spp.var = c("Catch.ich$log_cluhar","Catch.ich$log_urospp",
                                   "Catch.ich$log_polvir","Catch.ich$log_scosco",
                                   "Catch.ich$log_pepspp","lCatch.ich$og_ammspp")))

####################################
####################################
####################################
# CORRELATION COEFFICIENT
####################################
####################################
####################################
install.packages("dplyr")
install.packages("PerformanceAnalytics")
install.packages("Hmisc")
library("Hmisc")
library("PerformanceAnalytics")
library("dplyr")

##################################
#1. Read data and cross tabulate.#
##################################
Catch <- as.data.frame(read.csv("Plankton.csv",header=T))   # The output file "Plankton.csv" has an added first column (record number)
Catch$date <- as.Date(Catch$date, "%Y-%m-%d")		            # Assign Date format 'yyyy-mm-dd'
xtabs(~ year + month, data = Catch)                         # Crosstabulate the data by year and month

# 2. compute yearly means by species
Catch.year <- Catch[,c(23,26:31)]                           # Select columns with log(ichthyoplankton values), include year
Catch.byyear <- Catch.year %>%                              # Compute average CPUE by year and species
  group_by(year) %>%
  summarise_all(list(mean=mean), na.rm = TRUE)

###########################################################
# Correlation coefficients for two cases (lecture example)#
###########################################################
plot(Catch.byyear$year,Catch.byyear$log_urospp_mean,ylab="silver hake CPUE",xlab="year")
par(mfrow = c(1,2))
plot(Catch.byyear$log_urospp_mean,Catch.byyear$log_merbil_mean,ylab="silver hake",xlab="Red/white hake")
plot(Catch.byyear$log_urospp_mean,Catch.byyear$log_cluhar_mean,ylab="herring",xlab="Red/white hake")
rcorr(Catch.byyear$log_urospp_mean,Catch.byyear$log_merbil_mean, type = "pearson")
rcorr(Catch.byyear$log_urospp_mean,Catch.byyear$log_cluhar_mean, type = "pearson")
par(mfrow = c(1,1))

####################################################################
# 3. Plot data and compute correlation coefficients for all species#
####################################################################
pdf("correlation coefficient 3 output.pdf")
my_data <- Catch.byyear[,c(2:7)]                           # Select columns with log(ichthyoplankton values), include year
chart.Correlation(my_data, histogram=TRUE, pch=19)
dev.off()

# http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software
#In the above plot:
#  The distribution of each variable is shown on the diagonal.
# On the bottom of the diagonal : the bivariate scatter plots with a fitted line are displayed
# On the top of the diagonal : the value of the correlation plus the significance level as stars
# Each significance level is associated to a symbol : p-values(0, 0.001, 0.01, 0.05, 0.1, 1) <=> symbols(“***”, “**”, “*”, “.”, " “)

####################################
####################################
####################################
# STATISTICAL FORECASTING
####################################
####################################
####################################
install.packages("dplyr")
install.packages("PerformanceAnalytics")
install.packages("ggplot2")
install.packages("gridExtra")
library("dplyr")
library("PerformanceAnalytics")
library("ggplot2")
library("gridExtra")


###############
#1. Read data #
###############
my_data <- as.data.frame(read.csv("TempPlanktonCommontern.csv",header=T))   # The file "TempPlanktonCommontern.csv" has an added first column (record number)
my_data <- subset(my_data[,c(2:10)])                                        # Delete first column, which is extraneous

###########################
# 1a. Plot variable pairs #
###########################
pairs(my_data[,c(2:9)],upper.panel=NULL)
pdf("Statistical forecasting 1a output.pdf")
chart.Correlation(my_data, histogram=TRUE, pch=19)
dev.off()

##################################################################################
# 2. Linear regression between butterfish (pepspp) CPUE and surface temperature. #
##################################################################################
# first plot data
plot(my_data$sfc_temp_mean,my_data$log_pepspp_mean, xlab="surface temperature (deg C)", ylab="butterfish log10 CPUE")
cor(my_data$sfc_temp_mean,my_data$log_pepspp_mean)            # and correlation coefficient

#######################
# 2a linear regression#
#######################
# linear model
lm_pepspp <- lm(log_pepspp_mean ~ sfc_temp_mean, my_data)
summary(lm_pepspp)
AIC(lm_pepspp)   # lowest AIC value, best model

# test nonlinear model (add squared term)
my_data$sfc_temp_mean_squared <- my_data$sfc_temp_mean^2
lm_pepspp_nonlinear <- lm(log_pepspp_mean ~ sfc_temp_mean + sfc_temp_mean_squared, my_data)
summary(lm_pepspp_nonlinear)
AIC(lm_pepspp_nonlinear)   # higher AIC value, added parameter does not improve fit enough

##############################
# 2b plot fitted relationship#
##############################
pdf("Statistical forecasting 2b output.pdf")
par(mfrow = c(1,2))
# Plot observations and predicted relationship
plot(my_data$sfc_temp_mean,my_data$log_pepspp_mean, col = "blue", 
     xlab="surface temperature (deg C)", ylab="butterfish log10 CPUE",
     main="Butterfish pred vs surface temperature")
abline(lm_pepspp, col = "red") #Add a regression line
# Plot residuals and predicted relationship
plot(my_data$sfc_temp_mean,lm_pepspp$residuals, col = "blue", 
     xlab="surface temperature (deg C)", ylab="butterfish log10 CPUE residual",
     main="Butterfish residuals vs surface temperature")
abline(a=0,b=0,col="red")
dev.off()

# You can display the regression coefficients like so: 
coefficients(lm_pepspp)[1]                           # intercept
coefficients(lm_pepspp)[2]                           # slope
# and let's name them for use in the prediction
b0 <- coefficients(lm_pepspp)[1]                     # intercept
b1 <- coefficients(lm_pepspp)[2]                     # slope

# Examine variability of fit, for use in prediction
sd_pepspp <- sd(lm_pepspp$residuals)                 # sd of residuals, for use in prediction
                                                     # sd_pepspp
                                                     # [1] 0.02695914

###########################################################
# 3a. Forecast surface temperature and log10 butterfish CPUE#
############################################################################
# Generate random values for predicted temperature and log10 butterfish CPUE, both normally distributed
# Assume 0.03 deg C average increase per year 
temp_rate <- 0.03                                # rate of temperature (deg C) change per year
temp_sd <- 0.3                                   # standard deviation of the forecast temperature

# Check distributions of random numbers
hist(rnorm(1000, mean=temp_rate, sd=temp_sd))    # adjusted SD to approximate distribution of annual changes
                                                 # (Fig. 1, Pershing) (by eye, no published table of exact values)
                                                 # i.e., +- 0.5 deg C annual change is common
hist(rnorm(1000, mean=0, sd=sd_pepspp))          # examine SD for log10 butterfish CPUE residuals
                                                 # looks reasonable, simulated values +-0.05, similar to residuals

# Set up forecast
num_obs <- length(my_data$year)                      # number of observations in forecast
last_temp <- my_data$sfc_temp_mean[num_obs]          # last observation of surface temperature
last_year <- my_data$year[num_obs]                   # last observed year (available)

# Create table "pred" for one 30-year set of predictions
num_pred <- 30                                                          # number of prediction years
pred <- data.frame(matrix(nrow=num_pred,ncol=3,byrow=TRUE))             # create data frame (aka table) for predicted values
colnames(pred) <- c("year","sfc_temp_mean","log_pepspp_mean")           # label columns
pred$year <- seq(from = last_year+1, to = last_year+num_pred, by = 1)   # fill year column with prediction years

# Create a second table "sim" for predicted vaues for all simulations
num_sim <- 1000                                                                      # number of simulations
sim <- data.frame(matrix(nrow=num_sim*num_pred,ncol=4,byrow=TRUE))                  # create data frame (aka table) for predicted values
colnames(sim) <- c("sim","year","sfc_temp_mean","log_pepspp_mean")                  # label columns
sim$year <- rep(seq(from = last_year+1, to = last_year+num_pred, by = 1),num_sim)   # fill year column with prediction years

# Run simulations
for (i in 1:num_sim) {
  randvar_temp <- rnorm(num_pred, mean=temp_rate, sd=temp_sd)           # generate random variation for annual temperature changes
  randvar_pepspp <- rnorm(num_pred, mean=0, sd=sd_pepspp)               # generate random variation for annual butterfish CPUE

  sim$sim[(i-1)*30+1] <- i                                              # track simulation number 

# first year predictions
  pred$sfc_temp_mean[1] <- last_temp + randvar_temp[1]                  # year 1 time series based on last observation plus 
                                                                        #   first predicted annual change and some random variation
  pred$log_pepspp_mean[1] <- b0 + b1 * pred$sfc_temp_mean[1] +          # predicted log10 pepspp based on predicted temperature 
                              randvar_pepspp[1]                         #   plus random variation 

  sim$sfc_temp_mean[(i-1)*30+1] <- pred$sfc_temp_mean[1]                # surface temperature prediction stored in sim table
  sim$log_pepspp_mean[(i-1)*30+1] <- pred$log_pepspp_mean[1]            # butterfish prediction stored in sim table

# years 2 and onward predictions  
  for (j in 2:num_pred) {                                               # repeat computations for years 2-30
    sim$sim[(i-1)*30+j] <- i                                            # track simulation number
    pred$sfc_temp_mean[j] <- pred$sfc_temp_mean[j-1] + randvar_temp[j]  # predicted surface temperature
    pred$log_pepspp_mean[j] <- b0 + b1 * pred$sfc_temp_mean[j] +        # predicted log10 pepspp  
                                randvar_pepspp[j]                       
    
    sim$sfc_temp_mean[(i-1)*30+j] <- pred$sfc_temp_mean[j]              # store predicted values in the sim table
    sim$log_pepspp_mean[(i-1)*30+j] <- pred$log_pepspp_mean[j]
    }
}

# Compute mean and sd of predicted values by year
sim_temp_mean <- (with(sim, tapply(sfc_temp_mean,list(year=year),mean)))          # mean temperature
sim_temp_sd <- (with(sim, tapply(sfc_temp_mean,list(year=year),sd)))              # sd temperature
sim_pepspp_mean <- (with(sim, tapply(log_pepspp_mean,list(year=year),mean)))      # mean log10 butterfish CPUE
sim_pepspp_sd <- (with(sim, tapply(log_pepspp_mean,list(year=year),sd)))          # sd log10 butterfish CPUE

# Check that simulated values correctly mimic the assumed level of variability
mean(sim_temp_sd)   # > mean(sim_temp_sd)
                    # [1] 1.094438  (However the first year value is ~0.3, similar to sd of random numbers)

mean(sim_pepspp_sd)   # ~ mean(sim_pepspp_sd)
                      # [1] 0.02798031  (similar to sd of random numbers, which is ~0.03)

# Plot simulation results
# First, combine simulation results into a single table
temp_plus_sd <- sim_temp_mean + sim_temp_sd                       # add columns for +- one SD, surface temperature
temp_minus_sd <- sim_temp_mean - sim_temp_sd 
pepspp_plus_sd <- sim_pepspp_mean + sim_pepspp_sd                 # add columns for +- one SD, log10 butterfish CPUE
pepspp_minus_sd <- sim_pepspp_mean - sim_pepspp_sd 
sim_table <- data.frame(cbind(pred$year,sim_temp_mean, temp_plus_sd, temp_minus_sd,
                              sim_pepspp_mean, pepspp_plus_sd, pepspp_minus_sd))
colnames(sim_table) <- c("year","temp.mean","temp.plus.sd","temp.minus.sd","pepspp.mean","pepspp.plus.sd","pepspp.minus.sd")

# Use ggplot to display results (more flexible)
#   Plot mean and sd of predicted surface temperature
theme_set(theme_bw())  # pre-set the bw theme.
p1 <- ggplot(sim_table, aes(year, temp.mean)) + 
  geom_line(aes(y = temp.mean)) +
  geom_line(aes(y = temp.plus.sd), linetype="dashed") +
  geom_line(aes(y = temp.minus.sd), linetype="dashed") +
  labs(y = "Surface temperature (deg C)", x = "Year", title = "Simulations") +
  theme(
    axis.text.x = element_text(color="black", size=12, face="plain"),
    axis.text.y = element_text(color="black", size=12, face="plain"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold")
  )
p1

# Plot mean and sd of predicted log10 butterfish CPUE
p2 <- ggplot(sim_table, aes(year, pepspp.mean)) + 
  geom_line(aes(y = pepspp.mean)) +
  geom_line(aes(y = pepspp.plus.sd), linetype="dashed") +
  geom_line(aes(y = pepspp.minus.sd), linetype="dashed") +
  labs(y = "log10 butterfish CPUE", x = "Year", title = "Simulations") +
  theme(
    axis.text.x = element_text(color="black", size=12, face="plain"),
    axis.text.y = element_text(color="black", size=12, face="plain"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold")
  )
p2

# Plot one simulation, including observations
# First select columns and then combine observations with first simulation
# surface temperature
my_data_temp <- subset(my_data,select = c(year,sfc_temp_mean))       # select year and surface temperature from data table
colnames(my_data_temp) <- c("year","sfc_temp_mean")                  # label columns
one_sim_temp <- subset(sim, sim==1)                                  # select first simulation
one_sim_temp <- subset(one_sim_temp,select = c(year,sfc_temp_mean))  # select year and surface temperature
one_sim_temp <- rbind(my_data_temp,one_sim_temp)                     # combine observations and first simulation year

# log10 butterfish CPUE
my_data_pepspp <- subset(my_data,select = c(year,log_pepspp_mean))         # select year and log10 butterfish CPUE from data table
colnames(my_data_pepspp) <- c("year","log_pepspp_mean")                    # label columns
one_sim_pepspp <- subset(sim, sim==1)                                      # select first simulation
one_sim_pepspp <- subset(one_sim_pepspp,select = c(year,log_pepspp_mean))  # select year and log10 butterfish CPUE
one_sim_pepspp <- rbind(my_data_pepspp,one_sim_pepspp)                     # combine observatins and first simulation year

# Plot surface temperature
# Observations are much more variable than the predictions
# Likely an effect of the temporal and spatial variability of the observations
theme_set(theme_bw())  # pre-set the bw theme.
p3 <- ggplot(one_sim_temp, aes(year, sfc_temp_mean)) + 
  geom_line(aes(y = sfc_temp_mean)) +
  labs(y = "Surface temperature (deg C)", x = "Year", title = "Observed and one simulation") +
  theme(
    axis.text.x = element_text(color="black", size=12, face="plain"),
    axis.text.y = element_text(color="black", size=12, face="plain"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold")
  )
p3

# Plot log10 butterfish CPUE
theme_set(theme_bw())  # pre-set the bw theme.
p4 <- ggplot(one_sim_pepspp, aes(year, log_pepspp_mean)) + 
  geom_line(aes(y = log_pepspp_mean)) +
  labs(y = "log10 butterfish CPUE", x = "Year", title = "Observed and one simulation") +
  theme(
    axis.text.x = element_text(color="black", size=12, face="plain"),
    axis.text.y = element_text(color="black", size=12, face="plain"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold")
  )
p4

# Plot observations and all simulation results together
# First select columns and then combine observations with simulation average
my_data_temp <- subset(my_data,select = c(year,sfc_temp_mean))     # select year and surface temperature from observations
colnames(my_data_temp) <- c("year","temp.mean")                    # label columns
my_data_pepspp <- subset(my_data,select = c(log_pepspp_mean))      # select year and log10 butterfish CPUE from observations
colnames(my_data_pepspp) <- c("pepspp.mean")                       # label columns

temp1 <- data.frame(matrix(nrow=18,ncol=2,byrow=TRUE))             # create placeholder values for sd temperature (all NA values)
                                                                   # Need these placeholder columns for matching format of sim_table (below)
colnames(temp1) <- c("temp.plus.sd","temp.minus.sd")               # label columns

temp2 <- data.frame(matrix(nrow=18,ncol=2,byrow=TRUE))             # create placeholder values for sd butterfish CPUE (all NA values)
colnames(temp2) <- c("pepspp.plus.sd","pepspp.minus.sd")           # label columns

my_data_temp_pepspp <- cbind(my_data_temp,temp1,my_data_pepspp,temp2) # combine observation tables
obs_sim <- rbind(my_data_temp_pepspp,sim_table)                       # combine observations with simulations

# Plot temperature predictions
theme_set(theme_bw())  # pre-set the bw theme.
p5 <- ggplot(obs_sim, aes(year, temp.mean)) + 
  geom_line(aes(y = temp.mean)) +
  geom_line(aes(y = temp.plus.sd), linetype="dashed") +
  geom_line(aes(y = temp.minus.sd), linetype="dashed") +
  labs(y = "Surface temperature (deg C)", x = "Year", title = "Observed and all simulations") +
  theme(
    axis.text.x = element_text(color="black", size=12, face="plain"),
    axis.text.y = element_text(color="black", size=12, face="plain"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold")
  )
p5
# FYI, the warning messages are due to the placeholder (mising) values for sd of the observations

# Plot lgo10 butterfish CPUE predictions
theme_set(theme_bw())  # pre-set the bw theme.
p6 <- ggplot(obs_sim, aes(year, pepspp.mean)) + 
  geom_line(aes(y = pepspp.mean)) +
  geom_line(aes(y = pepspp.plus.sd), linetype="dashed") +
  geom_line(aes(y = pepspp.minus.sd), linetype="dashed") +
  labs(y = "log10 butterfish CPUE", x = "Year", title = "Observed and all simulations") +
  theme(
    axis.text.x = element_text(color="black", size=12, face="plain"),
    axis.text.y = element_text(color="black", size=12, face="plain"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold")
  )
p6


# Plot everything together!
pdf("Statistical forecasting 3a output.pdf")
grid.arrange(p3, p4, p5, p6, ncol = 2)
dev.off()

##########
# OPTIONAL
####################################
####################################
####################################
#PRINCIPAL COMPONENTS ANALYSIS (PCA)
####################################
####################################
####################################
install.packages("dplyr")
library(dplyr)

##################################
#1. Read data and cross tabulate.#
##################################
Catch <- as.data.frame(read.csv("Plankton.csv",header=T))   # The output file "Plankton.csv" has an added first column (record number)
Catch$date <- as.Date(Catch$date, "%Y-%m-%d")		            # Assign Date format 'yyyy-mm-dd'
xtabs(~ year + month, data = Catch)                         # Crosstabulate the data by year and month

#####################################
# 2. compute yearly means by species#
#####################################
Catch.year <- Catch[,c(23,26:31)]                           # Select columns with log(ichthyoplankton values), include year
Catch.byyear <- Catch.year %>%                              # Compute average CPUE by year and species
  group_by(year) %>%
  summarise_all(list(mean=mean), na.rm = TRUE)

####################################################################
# 3. PCA                                                           #
####################################################################
pairs(Catch.byyear[,c(2:7)],lower.panel=NULL)               # Plot data pairs
PCA <- prcomp(Catch.byyear[,c(2:7)], scale=TRUE)            # Run PCA
summary(PCA)                                                # PCA summary
biplot(PCA, scale=0)                                        # Plot first two principal components

pdf("PCA.pdf")
par(mfrow = c(1,1))
biplot(PCA, scale=0)                                        # Plot first two principal components
chart.Correlation(my_data, histogram=TRUE, pch=19)          # Correlation coefficients
dev.off()

############
# END OF LAB
############

####################################
####################################
####################################
# PREPARE FILE FOR STATISTICAL FORECASTING
####################################
####################################
####################################
install.packages("dplyr")
install.packages("PerformanceAnalytics")
library("dplyr")
library("PerformanceAnalytics")

##################################
#Read data and cross tabulate.#
##################################
Catch <- as.data.frame(read.csv("Plankton.csv",header=T))   # The file "Plankton.csv" has an added first column (record number)
Catch$date <- as.Date(Catch$date, "%Y-%m-%d")		            # Assign Date format 'yyyy-mm-dd'
xtabs(~ year + month, data = Catch)                         # Crosstabulate the data by year and month

Tern <- as.data.frame(read.csv("Tern productivity.csv",header=T))   # Tern productivity

##########################
#Compute yearly means#
##########################
Catch.year <- Catch[,c(23,26:31)]                           # Select columns with log(ichthyoplankton values), include year

# Compute average CPUE by year and species
Catch.byyear <- Catch.year %>%                              
  group_by(year) %>%
  summarise_all(list(mean=mean), na.rm = TRUE)

SprTemp.year <- subset(Catch,month>=4 &month<=6)            # Select spring only (April-June)
xtabs(~ year + month, data = SprTemp.year)                  # Crosstabulate the data by year and month, to check subset
SprTemp.year <- subset(SprTemp.year[,c(9,23)])              # Select surface temperature, include year

# Compute average spring temperature by year
SprTemp.byyear <- SprTemp.year %>%                          
  group_by(year) %>%
  summarise_all(list(mean=mean), na.rm = TRUE)
colnames(SprTemp.byyear) <- c("year","sfc_temp_mean")       # Change column name from "mean" to "sfc_temp_mean" to be consistent with ichthyoplankton naming

# Merge files by year
my_data <- merge(Catch.byyear,Tern, by.x="year", by.y="Year",all.x=TRUE)
my_data <- subset(my_data[,c(1:7,10)])                      # Keep only common tern number fledged per pair
my_data <- merge(my_data,SprTemp.byyear, by="year",all.x=TRUE)
my_data <- subset(my_data,year>=1999)

my_data <- merge(SprTemp.byyear, Catch.byyear, by="year",all.x=TRUE)
my_data <- merge(my_data,Tern, by.x="year", by.y="Year",all.x=TRUE)
my_data <- subset(my_data[,c(1:8,11)])                      # Keep only common tern number fledged per pair
my_data <- subset(my_data,year>=1999)

pairs(my_data[,c(2:9)],upper.panel=NULL)
chart.Correlation(my_data, histogram=TRUE, pch=19)

# Write output file (which will be used as input to lab for statistical forecasting)
write.csv(my_data, file = "TempPlanktonCommontern.csv")

####################################
####################################
####################################
# PREPARE FILE FOR CHRIS'S LECTURE AND LAB
####################################
####################################
####################################
install.packages("dplyr")
install.packages("PerformanceAnalytics")
library("dplyr")
library("PerformanceAnalytics")

##################################
#Read data and cross tabulate.#
##################################
Catch <- as.data.frame(read.csv("Plankton.csv",header=T))   # The file "Plankton.csv" has an added first column (record number)
Catch$date <- as.Date(Catch$date, "%Y-%m-%d")		            # Assign Date format 'yyyy-mm-dd'
xtabs(~ year + month, data = Catch)                         # Crosstabulate the data by year and month

Tern <- as.data.frame(read.csv("Tern productivity.csv",header=T))   # Tern productivity

##########################
#Compute yearly means#
##########################
Catch.year <- Catch[,c(23,26:31)]                           # Select columns with log(ichthyoplankton values), include year

# Compute average CPUE by year and species
Catch.byyear <- Catch.year %>%                              
  group_by(year) %>%
  summarise_all(list(mean=mean), na.rm = TRUE)

SprTemp.year <- subset(Catch,month>=4 &month<=6)            # Select spring only (April-June)
xtabs(~ year + month, data = SprTemp.year)                  # Crosstabulate the data by year and month, to check subset
SprTemp.year <- subset(SprTemp.year[,c(9,23)])              # Select surface temperature, include year

# Compute average spring temperature by year
SprTemp.byyear <- SprTemp.year %>%                          
  group_by(year) %>%
  summarise_all(list(mean=mean), na.rm = TRUE)
colnames(SprTemp.byyear) <- c("year","sfc_temp_mean")       # Change column name from "mean" to "sfc_temp_mean" to be consistent with ichthyoplankton naming

# Merge files by year
my_data <- merge(Catch.byyear,Tern, by.x="year", by.y="Year",all.x=TRUE)
my_data <- subset(my_data[,c(1:7,10,14)])                      # Keep only number fledged per pair for both common and roseate tern
my_data <- merge(my_data,SprTemp.byyear, by="year",all.x=TRUE)
my_data <- subset(my_data,year>=1999)

my_data <- merge(SprTemp.byyear, Catch.byyear, by="year",all.x=TRUE)
my_data <- merge(my_data,Tern, by.x="year", by.y="Year",all.x=TRUE)
my_data <- subset(my_data[,c(1:8,11,15)])                      # Keep only number fledged per pair for common tern and roseate tern
my_data <- subset(my_data,year>=1999)

pairs(my_data[,c(2:9)],upper.panel=NULL)
chart.Correlation(my_data, histogram=TRUE, pch=19)

# Write output file (which will be used as input to Chris's lab exercises)
write.csv(my_data, file = "TempPlanktonCommonRoseatetern.csv")


##################################
##################################
##################################
# Experiment with GENERAL ADDITIVE MODEL (GAM)
##################################
##################################
##################################
library(fields)
library(dplyr)
library(mgcv)
library(mapdata)
library(mapproj)

# Appledore Island location
Lat.SML <- 42.987727
Lon.SML <- -70.613940

##################################
#1. Read data and cross tabulate.#
##################################
Catch <- as.data.frame(read.csv("Plankton.csv",header=T))   # The output file "Plankton.csv" has an added first column (record number)
Catch$date <- as.Date(Catch$date, "%Y-%m-%d")		            # Assign Date format 'yyyy-mm-dd'
xtabs(~ year + month, data = Catch)                         # Crosstabulate the data by year and month

start.year <- 1989  # all years, can change start year if desired (e.g., Ecomon years: 1989-present)
Catch.ich <- subset(Catch, year>=start.year)         

#1.5 Estimate ichthyoplankton abundance estimates 
Catch.ich <- subset(Catch.ich, distance < 150)   # select samples within 150 km of Appledore Island, i.e., the Gulf of Maine

# TESTING: no spatial component
# For each species, fit GAM as a function of year and month; summarize results and plot diagnostics (gam.check)
# year is a smoothed effect
par(mfrow = c(2,2))   # format plot layout (2 rows and 2 columns of plots per page)
fit.cluhar <- gam(log_cluhar ~ s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.cluhar); gam.check(fit.cluhar)
fit.urospp <- gam(log_urospp ~ s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.urospp); gam.check(fit.urospp)
fit.merbil <- gam(log_merbil ~ s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.merbil); gam.check(fit.merbil)
fit.scosco <- gam(log_scosco ~ s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.scosco); gam.check(fit.scosco)
fit.pepspp <- gam(log_pepspp ~ s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.pepspp); gam.check(fit.pepspp)
fit.ammspp <- gam(log_ammspp ~ s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.ammspp); gam.check(fit.ammspp)
par(mfrow = c(1,1))  # reset plot layout

# TESTING: year as random effect
par(mfrow = c(2,2))   # format plot layout (2 rows and 2 columns of plots per page)
fit.cluhar <- gam(log_cluhar ~ s(year,bs="re") + s(month,bs='cc'), data=Catch.ich); summary(fit.cluhar); gam.check(fit.cluhar) # nonsensical result
par(mfrow = c(1,1))  # reset plot layout

# Output plotted results to a file
pdf("GAM 1.5 output.pdf")
par(mfrow = c(2,3))

newd <- data.frame(year=(1989:2016),month=rep.int(6,2016-1989+1))
year <- (1989:2016)
pred.cluhar <- predict.gam(fit.cluhar,newd)
pred.urospp <- predict.gam(fit.urospp,newd)
pred.merbil <- predict.gam(fit.merbil,newd)
pred.scosco <- predict.gam(fit.scosco,newd)
pred.pepspp <- predict.gam(fit.pepspp,newd)
pred.ammspp <- predict.gam(fit.ammspp,newd)

plot(year,pred.cluhar)
plot(year,pred.urospp)
plot(year,pred.merbil)
plot(year,pred.scosco)
plot(year,pred.pepspp)
plot(year,pred.ammspp)

#by year
plot(fit.cluhar,pages=0,seWithMean=TRUE,select=1,ylim=c(-0.3,0.3),main="herring")
plot(fit.urospp,pages=0,seWithMean=TRUE,select=1,ylim=c(-0.2,0.2),main="red/white hake")
plot(fit.merbil,pages=0,seWithMean=TRUE,select=1,ylim=c(-0.1,0.1),main="silver hake")
plot(fit.scosco,pages=0,seWithMean=TRUE,select=1,ylim=c(-0.1,0.1),main="mackerel")
plot(fit.pepspp,pages=0,seWithMean=TRUE,select=1,ylim=c(-0.1,0.1),main="butterfish")
plot(fit.ammspp,pages=0,seWithMean=TRUE,select=1,ylim=c(-0.4,0.4),main="sandlance")

#by month
plot(fit.cluhar,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.3,0.3),main="herring")
plot(fit.urospp,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.4,0.4),main="red/white hake")
plot(fit.merbil,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.3,0.3),main="silver hake")
plot(fit.scosco,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.3,0.3),main="mackerel")
plot(fit.pepspp,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.3,0.3),main="butterfish")
plot(fit.ammspp,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.6,0.6),main="sandlance")

dev.off()

###############################
# Modify code from Ingrid
library("PerformanceAnalytics")

pdf("Yearly_predictions.pdf")
par(mfrow = c(6,2))

year <- seq(1989,2016)
pred.allspp <- data.frame(cbind(year,pred.ammspp,pred.cluhar,pred.merbil,pred.pepspp,pred.scosco,pred.urospp))
colnames(pred.allspp) <- c("Year","sandlance","herring","silver hake","butterfish","mackerel","red/white hake")
my_data <- pred.allspp[,c(2:7)]                           # Select columns with ichthyoplankton values
str(pred.allspp)

chart.Correlation(my_data, histogram=TRUE, pch=19)

# Add tern data
Tern <- as.data.frame(read.csv("Tern productivity.csv",header=T))   # Tern productivity
Tern.years <- subset(Tern,Year>=1999&Year<=2016)
pred.allspp <- subset(pred.allspp,Year>=1999&Year<=2016)
Tern.pred <- cbind(Tern.years,pred.allspp)
str(Tern.pred)

my_data <- Tern.pred[,c("COTE_numberfledged","sandlance","herring","silver hake","butterfish","mackerel","red/white hake")]
str(my_data)

chart.Correlation(my_data, histogram=TRUE, pch=19)

# Code sent by Ingrid for me to test, on 7/5/2020
pdf("Yearly Predictions .pdf")
par(mfrow = c(6,2))

# This line misidentifies year as character, which causes correlation coefficient to fail
# year <- c("1989","1990","1991","1992","1993","1994","1995","1996","1997","1998","1999","2000","2001","2002","2003",
#          "2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016")
# Instead use
year <- seq(1989,2016)

pred.allspp <- data.frame(cbind(year, pred.ammspp,pred.cluhar,pred.merbil,pred.pepspp,pred.scosco,pred.urospp))
colnames(pred.allspp) <- c("year", "sandlance","herring","silver hake","butterfish","mackeral","red/white hake")

pred.allspp

my_data <- pred.allspp[c(2:7)]  #selecting columns with ichthyoplankton data
str(pred.allspp)

chart.Correlation(my_data, histogram=TRUE, pch=19)

dev.off()

#############################################################################################################
# TESTING: include spatial component, draw polygons on spatial output
# From Chris Rooper, 6/22/2020: If you’re doing it in R, a way to do it would be to
# 1)     Make a spatial point for the center of your circle
# 2)     Buffer the point with a radius of 150 km (turning it into a spatial polygon)
# 3)     Use the raster::extract function to pull a mean from the overlapping polygon and CPUE rasters
# Items 1 and 2 can be done using the sp package I think. I probably have some code that would get you close.

library(fields)
library(dplyr)
library(mgcv)
library(mapdata)
library(mapproj)
library(rgeos)
library(sp)
library(raster)
library(rgdal)

# Appledore Island location
Lat.SML <- 42.987727
Lon.SML <- -70.613940

# Code from Chris Rooper
Points<-SpatialPoints(cbind(-122.5,58.5),proj4string=CRS("+proj=longlat +datum=WGS84"))
points.project <- spTransform(Points, "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-80 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 
                              +units=m +no_defs")

#Bufferedpoint<-gBuffer(points.project,byid=TRUE,width=150,capStyle=”ROUND”) # 150 m buffer
Bufferedpoint<-gBuffer(points.project,byid=TRUE,width=150) # 150 m buffer

Meanvalue<-raster::extract(yourraster,Bufferedpoint,fun=”mean”) # Then you take your raster layer and extract the mean value inside the buffered point

# My adapted code
Points<-SpatialPoints(cbind(Lon.SML,Lat.SML),
                      proj4string=CRS("+proj=longlat +datum=WGS84")) # create object of class "SpatialPoints-class" (irregularly spaced points)
                                                                     # CRS = coordinate reference system
points.project <- spTransform(Points, "+proj=aea +lat_1=45 +lat_2=55 +lat_0=40 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 
                              +units=m +no_defs")                    # Transforms "Points" to an albers equal area projection.
                                                                     #   This will allow you to use a 150 m buffer.
Bufferedpoint<-gBuffer(points.project,byid=TRUE,width=150)           # 150 m buffer

Meanvalue<-raster::extract(raster.gam,Bufferedpoint,fun=mean)      # Then you take your raster layer and extract the mean value inside the buffered point


##########################################
###TO DO NEXT: figure out how to convert the spatial gam predictions to a raster, presumably using the same projection
# Create raster layer with a prediction layer basedon a fitted model object
# See https://www.rdocumentation.org/packages/raster/versions/1.7-18/topics/predict
# and https://rspatial.org/raster/spatial/8-rastermanip.html

# Output from a simple, spatial GAM (see below)
vis.gam(fit.cluhar, plot.type="contour", too.far=0.15, main="Atlantic herring", xlab="Longitude",ylab="Latitude",color="topo")
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")

raster.gam<-vis.gam(fit.cluhar, plot.type="contour", too.far=0.15, main="Atlantic herring", xlab="Longitude",ylab="Latitude",color="topo")
raster.gam
#NULL
raster.gam<-predict(fit.cluhar)
str(raster.gam)

# GAM by year, month, and location (all years)
##################################
#1. Read data and cross tabulate.#
##################################
Catch <- as.data.frame(read.csv("Plankton.csv",header=T))   # The output file "Plankton.csv" has an added first column (record number)
Catch$date <- as.Date(Catch$date, "%Y-%m-%d")		            # Assign Date format 'yyyy-mm-dd'
xtabs(~ year + month, data = Catch)                         # Crosstabulate the data by year and month

start.year <- 1989  # all years, can change start year if desired (e.g., Ecomon years: 1989-present)
Catch.ich <- subset(Catch, year>=start.year)         

####################################################################################
#2. Fit GAM to log-transformed values, log10(value+1), include year and month terms#
# One spatial map for each species covering all years (i.e., average location)     #
#  run time about 1 minute                                                         #
####################################################################################

start.year <- 1989  # all years, can change start year if desired (e.g., Ecomon years: 1989-present)
Catch.ich <- subset(Catch, year>=start.year)         

# For each species, fit GAM as a function of location, year, and month; summarize results and plot diagnostics (gam.check)
par(mfrow = c(2,2))   # format plot layout (2 rows and 2 columns of plots per page)
fit.cluhar <- gam(log_cluhar ~ s(lon, lat) + s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.cluhar); gam.check(fit.cluhar)
fit.urospp <- gam(log_urospp ~ s(lon, lat) + s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.urospp); gam.check(fit.urospp)
fit.merbil <- gam(log_merbil ~ s(lon, lat) + s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.merbil); gam.check(fit.merbil)
fit.scosco <- gam(log_scosco ~ s(lon, lat) + s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.scosco); gam.check(fit.scosco)
fit.pepspp <- gam(log_pepspp ~ s(lon, lat) + s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.pepspp); gam.check(fit.pepspp)
fit.ammspp <- gam(log_ammspp ~ s(lon, lat) + s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.ammspp); gam.check(fit.ammspp)
par(mfrow = c(1,1))  # reset plot layout

# Output plotted results to a file
pdf("GAM 2 output.pdf")
par(mfrow = c(2,3))

#by year
plot(fit.cluhar,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.3,0.3),main="herring")
plot(fit.urospp,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.2,0.2),main="red/white hake")
plot(fit.merbil,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.1,0.1),main="silver hake")
plot(fit.scosco,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.1,0.1),main="mackerel")
plot(fit.pepspp,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.1,0.1),main="butterfish")
plot(fit.ammspp,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.4,0.4),main="sandlance")

#by month
plot(fit.cluhar,pages=0,seWithMean=TRUE,select=3,ylim=c(-0.3,0.3),main="herring")
plot(fit.urospp,pages=0,seWithMean=TRUE,select=3,ylim=c(-0.4,0.4),main="red/white hake")
plot(fit.merbil,pages=0,seWithMean=TRUE,select=3,ylim=c(-0.3,0.3),main="silver hake")
plot(fit.scosco,pages=0,seWithMean=TRUE,select=3,ylim=c(-0.3,0.3),main="mackerel")
plot(fit.pepspp,pages=0,seWithMean=TRUE,select=3,ylim=c(-0.3,0.3),main="butterfish")
plot(fit.ammspp,pages=0,seWithMean=TRUE,select=3,ylim=c(-0.6,0.6),main="sandlance")

# spatial plots
vis.gam(fit.cluhar, plot.type="contour", too.far=0.15, main="Atlantic herring", xlab="Longitude",ylab="Latitude",color="topo")
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
vis.gam(fit.urospp, plot.type="contour", too.far=0.15, main="red white hake", xlab="Longitude",ylab="Latitude",color="topo")
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
vis.gam(fit.merbil, plot.type="contour", too.far=0.15, main="silver hake", xlab="Longitude",ylab="Latitude",color="topo")
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
vis.gam(fit.scosco, plot.type="contour", too.far=0.15, main="Atlantic mackerel", xlab="Longitude",ylab="Latitude",color="topo")
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
vis.gam(fit.pepspp, plot.type="contour", too.far=0.15, main="butterfish", xlab="Longitude",ylab="Latitude",color="topo")
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
vis.gam(fit.ammspp, plot.type="contour", too.far=0.15, main="sandlance", xlab="Longitude",ylab="Latitude",color="topo")
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")

par(mfrow = c(1,1))
dev.off()


=======
# Program for statistics lab
# Michael Sigler
# May 19, 2019, Revised November 22, 2019, March 24, 2020, February 12, 2021

# Northeast US continental shelf plankton data, NOAA Northeast Fisheries Science Center
# Sampling years
#   1977-1987, MARMAP
#   1988-1998, some sampling occurred, but less balanced seasonally
#   1999-2017, Ecosystem Monitoring (EcoMon)
#   2018-2019, sampling occurred but data currently unavailable publicly
#   In some years, sampling was monthly, but in recent years, sampling focused on April-June and August-November
# Focus on 6 fish species (ichthyoplankton) that are primary prey of terns

# Call packages (libraries) used in data analysis; run every time
library(dplyr)
library(ggplot2)
library(tidyr)
library(knitr)
library(lubridate)
library(mgcv)
library(fields)
library(mapdata)
library(mapproj)

############################################################################################################
# Read plankton file
# This file is compiled from the original database available at:
# https://www.ncei.noaa.gov/metadata/geoportal/rest/metadata/item/gov.noaa.nodc%3A0187513/html
# Key variables
# lat = latitude
# lon = longitude
# date = date (GMT)
# sfc_temp = surface temperature
# sfc_salt = surface salinity
# volume_100m3 = Zooplankton Displacement Volume (ml)  per 1m2 of surface area
# zoop_100m3 = Zooplankton count per 100 m3 of volume
# fish_100m3 = Fish count per 100 m3 of volume
# merbil_100m3 = silver hake count per 100 m3 of volume
# cluhar_100m3 = Atlantic herring count per 100 m3 of volume
# ammspp_100m3 = sandlance count per 100 m3 of volume
# pepspp_100m3 = butterfish count per 100 m3 of volume
# scosco_100m3 = Atlantic mackerel count per 100 m3 of volume
# urospp_100m3 = red and white hake count per 100 m3 of volume
# e.g., log_cluhar = log 10 transformed (value +1) of cluhar_100m3
# distance = distance (km) from sample location to Appledore Island 

rm(list = ls())

# Appledore Island location
Lat.SML <- 42.987727
Lon.SML <- -70.613940
appledore.location <- c(Lon.SML,Lat.SML)

##################################
#1. Read data and cross tabulate.#
##################################
# set project directories
projdir <- "D:\\Documents\\Shoals\\2021 class\\Analysis\\plankton\\"
projdir.dat <- "D:\\Documents\\Shoals\\2021 class\\Analysis\\data\\"
projdir.results <- "D:\\Documents\\Shoals\\2021 class\\Analysis\\results\\"

# read data
filename <- paste(projdir.dat,"Plankton.csv",sep="") # filename
Catch <- as_tibble(read.csv(filename,header=T))      # read file as a "tibble"
Catch$date <- as.Date(Catch$date, "%Y-%m-%d")		     # Assign Date format 'yyyy-mm-dd'

# Crosstabulate the data by year and month
Catch %>%
  group_by(year, month)%>%                      # group by year and month
  summarise(n=n())%>%                           # summarise by count n=n()
  spread(month, n)%>%                           # spread month as header of table
  kable()                                       # format table

####################################################################################
#2. Fit General Additive Model (GAM) to log-transformed values, log10(value+1), include year and month terms#
# One spatial map for each species covering all years (i.e., average location)     #
#  run time about 1 minute                                                         #
####################################################################################

start.year <- 1977  # all years
Catch.ich <- subset(Catch, year>=start.year)         

# First GAM, fit as a function of location, year and month
# Atlantic herring (Clupea harengus)
# run GAM
fit.cluhar <- gam(log_cluhar ~ s(lon, lat) + s(year) + s(month,bs='cc'), data=Catch.ich) 

# diagnostics
summary(fit.cluhar)$r.sq                                   # R-squared: proportion of variance explained
hist(residuals(fit.cluhar, type = "deviance"),             # histogram of residuals
     xlab = "Residuals", main = "Histogram of residuals")

# plot results
plot(fit.cluhar,pages=0,seWithMean=TRUE,select=2,
     ylim=c(-0.3,0.3),main="herring")                      # year plot
plot(fit.cluhar,pages=0,seWithMean=TRUE,select=3,
     ylim=c(-0.3,0.3),main="herring")                      # month plot
vis.gam(fit.cluhar, plot.type="contour", too.far=0.15,     # map
        main="Atlantic herring", xlab="Longitude",ylab="Latitude",color="topo")
points(Lon.SML, Lat.SML, pch=16, col = "red")
map('worldHires',fill=T,add=T, col="grey")


fit.urospp <- gam(log_urospp ~ s(lon, lat) + s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.urospp); gam.check(fit.urospp)
fit.merbil <- gam(log_merbil ~ s(lon, lat) + s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.merbil); gam.check(fit.merbil)
fit.scosco <- gam(log_scosco ~ s(lon, lat) + s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.scosco); gam.check(fit.scosco)
fit.pepspp <- gam(log_pepspp ~ s(lon, lat) + s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.pepspp); gam.check(fit.pepspp)
fit.ammspp <- gam(log_ammspp ~ s(lon, lat) + s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.ammspp); gam.check(fit.ammspp)
par(mfrow = c(1,1))  # reset plot layout

# Output plotted results to a file
pdf("GAM 2 output.pdf")
par(mfrow = c(2,3))

#by year
plot(fit.cluhar,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.3,0.3),main="herring")
plot(fit.urospp,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.2,0.2),main="red/white hake")
plot(fit.merbil,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.1,0.1),main="silver hake")
plot(fit.scosco,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.1,0.1),main="mackerel")
plot(fit.pepspp,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.1,0.1),main="butterfish")
plot(fit.ammspp,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.4,0.4),main="sandlance")

#by month
plot(fit.cluhar,pages=0,seWithMean=TRUE,select=3,ylim=c(-0.3,0.3),main="herring")
plot(fit.urospp,pages=0,seWithMean=TRUE,select=3,ylim=c(-0.4,0.4),main="red/white hake")
plot(fit.merbil,pages=0,seWithMean=TRUE,select=3,ylim=c(-0.3,0.3),main="silver hake")
plot(fit.scosco,pages=0,seWithMean=TRUE,select=3,ylim=c(-0.3,0.3),main="mackerel")
plot(fit.pepspp,pages=0,seWithMean=TRUE,select=3,ylim=c(-0.3,0.3),main="butterfish")
plot(fit.ammspp,pages=0,seWithMean=TRUE,select=3,ylim=c(-0.6,0.6),main="sandlance")

# spatial plots
vis.gam(fit.cluhar, plot.type="contour", too.far=0.15, main="Atlantic herring", xlab="Longitude",ylab="Latitude",color="topo")
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
vis.gam(fit.urospp, plot.type="contour", too.far=0.15, main="red white hake", xlab="Longitude",ylab="Latitude",color="topo")
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
vis.gam(fit.merbil, plot.type="contour", too.far=0.15, main="silver hake", xlab="Longitude",ylab="Latitude",color="topo")
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
vis.gam(fit.scosco, plot.type="contour", too.far=0.15, main="Atlantic mackerel", xlab="Longitude",ylab="Latitude",color="topo")
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
vis.gam(fit.pepspp, plot.type="contour", too.far=0.15, main="butterfish", xlab="Longitude",ylab="Latitude",color="topo")
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
vis.gam(fit.ammspp, plot.type="contour", too.far=0.15, main="sandlance", xlab="Longitude",ylab="Latitude",color="topo")
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")

par(mfrow = c(1,1))
dev.off()

############################################################
#3. Fit GAM to log-transformed values, accounting for month#
# Spatial model by year groups (here 1-year groups)        #
############################################################
############################################################
############################################################
# I NEED TO SHORTEN THIS EXERCISE SO THAT THE CODE RUNS FASTER ON A STUDENT COMPUTER
############################################################
############################################################
############################################################

# Read data
Catch <- as.data.frame(read.csv("Plankton.csv",header=T))   # The output file "Plankton.csv" has an added first column (record number)
Catch$date <- as.Date(Catch$date, "%Y-%m-%d")		            # Assign Date format 'yyyy-mm-dd'
xtabs(~ year + month, data = Catch)                         # Crosstabulate the data by year and month

# Table of species names
name.tbl <- data.frame(spp.num = seq(26,31),spp.name = c("herring","red white hake","silver hake","mackerel","butterfish","sandlance"))

# Choose analysis years
start.year <- 1999                                          # Ecomon years: 1999-present
end.year <- 2016                                            # Last year in data set
Catch.ich <- subset(Catch, year>=start.year) 
xtabs(~ year + month, data = Catch.ich)                     # Crosstabulate the data by year and month

# Year grouping
#year.int <- 3                                                        # Size of year group (1, 3, 5, 10)
year.int <- 1                                                         # Size of year group (1, 3, 5, 10)

# Compute some indices to control the next program steps (GAM and plotting)
Catch.ich$year.group <- ceiling((Catch.ich$year-1971+1)/year.int)     # create year grouping (e.g., 1973 is a member of year.group 1 (1971-1975))
table(Catch.ich$year.group,Catch.ich$month)                           # tabulate year group by month
min.year.group <- min(Catch.ich$year.group)
max.year.group <- max(Catch.ich$year.group)
num.year.group <- max.year.group - min.year.group + 1                 # determine the number of year groups

# Run model (run time is several minutes)
pdf("GAM 3 output.pdf")
par(mfrow = c(2,3))

# This for loop goes through each species to run the gam, e.g., 26=cluhar, 27=urospp, etc.
for (i in 26:31) {  
  j <- i - 25   # counter for species names
  fit.spp <- gam(Catch.ich[,i] ~ s(lon, lat, by=as.factor(year.group),id=0) + s(month,bs='cc'), data=Catch.ich, dist ="normal")   # GAM
  summary(fit.spp); # gam.check(fit.spp)
  plot(fit.spp,pages=0,seWithMean=TRUE,select=num.year.group+1,ylim=c(-0.4,0.4),main=paste(name.tbl[j,2],"by month"))     # by month plot
  
  # by year group plots  
  for (k in min.year.group:max.year.group) {
    vis.gam(fit.spp, view=c("lon","lat"), plot.type="contour", too.far=0.15, cond=list(year.group=k),
            main=paste(name.tbl[j,2]," year group ",k*year.int+1971-1), xlab="Longitude",ylab="Latitude",color="topo") 
    points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
  }
}
par(mfrow = c(1,1))
dev.off()

####################################
####################################
####################################
# CORRELATION COEFFICIENT
####################################
####################################
####################################
install.packages("dplyr")
install.packages("PerformanceAnalytics")
install.packages("Hmisc")
library("Hmisc")
library("PerformanceAnalytics")
library("dplyr")

##################################
#1. Read data and cross tabulate.#
##################################
Catch <- as.data.frame(read.csv("Plankton.csv",header=T))   # The output file "Plankton.csv" has an added first column (record number)
Catch$date <- as.Date(Catch$date, "%Y-%m-%d")		            # Assign Date format 'yyyy-mm-dd'
xtabs(~ year + month, data = Catch)                         # Crosstabulate the data by year and month

# 2. compute yearly means by species
Catch.year <- Catch[,c(23,26:31)]                           # Select columns with log(ichthyoplankton values), include year
Catch.byyear <- Catch.year %>%                              # Compute average CPUE by year and species
  group_by(year) %>%
  summarise_all(list(mean=mean), na.rm = TRUE)

###########################################################
# Correlation coefficients for two cases (lecture example)#
###########################################################
plot(Catch.byyear$year,Catch.byyear$log_urospp_mean,ylab="silver hake CPUE",xlab="year")
par(mfrow = c(1,2))
plot(Catch.byyear$log_urospp_mean,Catch.byyear$log_merbil_mean,ylab="silver hake",xlab="Red/white hake")
plot(Catch.byyear$log_urospp_mean,Catch.byyear$log_cluhar_mean,ylab="herring",xlab="Red/white hake")
rcorr(Catch.byyear$log_urospp_mean,Catch.byyear$log_merbil_mean, type = "pearson")
rcorr(Catch.byyear$log_urospp_mean,Catch.byyear$log_cluhar_mean, type = "pearson")
par(mfrow = c(1,1))

####################################################################
# 3. Plot data and compute correlation coefficients for all species#
####################################################################
pdf("correlation coefficient 3 output.pdf")
my_data <- Catch.byyear[,c(2:7)]                           # Select columns with log(ichthyoplankton values), include year
chart.Correlation(my_data, histogram=TRUE, pch=19)
dev.off()

# http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software
#In the above plot:
#  The distribution of each variable is shown on the diagonal.
# On the bottom of the diagonal : the bivariate scatter plots with a fitted line are displayed
# On the top of the diagonal : the value of the correlation plus the significance level as stars
# Each significance level is associated to a symbol : p-values(0, 0.001, 0.01, 0.05, 0.1, 1) <=> symbols(“***”, “**”, “*”, “.”, " “)

####################################
####################################
####################################
# STATISTICAL FORECASTING
####################################
####################################
####################################
install.packages("dplyr")
install.packages("PerformanceAnalytics")
install.packages("ggplot2")
install.packages("gridExtra")
library("dplyr")
library("PerformanceAnalytics")
library("ggplot2")
library("gridExtra")


###############
#1. Read data #
###############
my_data <- as.data.frame(read.csv("TempPlanktonCommontern.csv",header=T))   # The file "TempPlanktonCommontern.csv" has an added first column (record number)
my_data <- subset(my_data[,c(2:10)])                                        # Delete first column, which is extraneous

###########################
# 1a. Plot variable pairs #
###########################
pairs(my_data[,c(2:9)],upper.panel=NULL)
pdf("Statistical forecasting 1a output.pdf")
chart.Correlation(my_data, histogram=TRUE, pch=19)
dev.off()

##################################################################################
# 2. Linear regression between butterfish (pepspp) CPUE and surface temperature. #
##################################################################################
# first plot data
plot(my_data$sfc_temp_mean,my_data$log_pepspp_mean, xlab="surface temperature (deg C)", ylab="butterfish log10 CPUE")
cor(my_data$sfc_temp_mean,my_data$log_pepspp_mean)            # and correlation coefficient

#######################
# 2a linear regression#
#######################
# linear model
lm_pepspp <- lm(log_pepspp_mean ~ sfc_temp_mean, my_data)
summary(lm_pepspp)
AIC(lm_pepspp)   # lowest AIC value, best model

# test nonlinear model (add squared term)
my_data$sfc_temp_mean_squared <- my_data$sfc_temp_mean^2
lm_pepspp_nonlinear <- lm(log_pepspp_mean ~ sfc_temp_mean + sfc_temp_mean_squared, my_data)
summary(lm_pepspp_nonlinear)
AIC(lm_pepspp_nonlinear)   # higher AIC value, added parameter does not improve fit enough

##############################
# 2b plot fitted relationship#
##############################
pdf("Statistical forecasting 2b output.pdf")
par(mfrow = c(1,2))
# Plot observations and predicted relationship
plot(my_data$sfc_temp_mean,my_data$log_pepspp_mean, col = "blue", 
     xlab="surface temperature (deg C)", ylab="butterfish log10 CPUE",
     main="Butterfish pred vs surface temperature")
abline(lm_pepspp, col = "red") #Add a regression line
# Plot residuals and predicted relationship
plot(my_data$sfc_temp_mean,lm_pepspp$residuals, col = "blue", 
     xlab="surface temperature (deg C)", ylab="butterfish log10 CPUE residual",
     main="Butterfish residuals vs surface temperature")
abline(a=0,b=0,col="red")
dev.off()

# You can display the regression coefficients like so: 
coefficients(lm_pepspp)[1]                           # intercept
coefficients(lm_pepspp)[2]                           # slope
# and let's name them for use in the prediction
b0 <- coefficients(lm_pepspp)[1]                     # intercept
b1 <- coefficients(lm_pepspp)[2]                     # slope

# Examine variability of fit, for use in prediction
sd_pepspp <- sd(lm_pepspp$residuals)                 # sd of residuals, for use in prediction
                                                     # sd_pepspp
                                                     # [1] 0.02695914

###########################################################
# 3a. Forecast surface temperature and log10 butterfish CPUE#
############################################################################
# Generate random values for predicted temperature and log10 butterfish CPUE, both normally distributed
# Assume 0.03 deg C average increase per year 
temp_rate <- 0.03                                # rate of temperature (deg C) change per year
temp_sd <- 0.3                                   # standard deviation of the forecast temperature

# Check distributions of random numbers
hist(rnorm(1000, mean=temp_rate, sd=temp_sd))    # adjusted SD to approximate distribution of annual changes
                                                 # (Fig. 1, Pershing) (by eye, no published table of exact values)
                                                 # i.e., +- 0.5 deg C annual change is common
hist(rnorm(1000, mean=0, sd=sd_pepspp))          # examine SD for log10 butterfish CPUE residuals
                                                 # looks reasonable, simulated values +-0.05, similar to residuals

# Set up forecast
num_obs <- length(my_data$year)                      # number of observations in forecast
last_temp <- my_data$sfc_temp_mean[num_obs]          # last observation of surface temperature
last_year <- my_data$year[num_obs]                   # last observed year (available)

# Create table "pred" for one 30-year set of predictions
num_pred <- 30                                                          # number of prediction years
pred <- data.frame(matrix(nrow=num_pred,ncol=3,byrow=TRUE))             # create data frame (aka table) for predicted values
colnames(pred) <- c("year","sfc_temp_mean","log_pepspp_mean")           # label columns
pred$year <- seq(from = last_year+1, to = last_year+num_pred, by = 1)   # fill year column with prediction years

# Create a second table "sim" for predicted vaues for all simulations
num_sim <- 1000                                                                      # number of simulations
sim <- data.frame(matrix(nrow=num_sim*num_pred,ncol=4,byrow=TRUE))                  # create data frame (aka table) for predicted values
colnames(sim) <- c("sim","year","sfc_temp_mean","log_pepspp_mean")                  # label columns
sim$year <- rep(seq(from = last_year+1, to = last_year+num_pred, by = 1),num_sim)   # fill year column with prediction years

# Run simulations
for (i in 1:num_sim) {
  randvar_temp <- rnorm(num_pred, mean=temp_rate, sd=temp_sd)           # generate random variation for annual temperature changes
  randvar_pepspp <- rnorm(num_pred, mean=0, sd=sd_pepspp)               # generate random variation for annual butterfish CPUE

  sim$sim[(i-1)*30+1] <- i                                              # track simulation number 

# first year predictions
  pred$sfc_temp_mean[1] <- last_temp + randvar_temp[1]                  # year 1 time series based on last observation plus 
                                                                        #   first predicted annual change and some random variation
  pred$log_pepspp_mean[1] <- b0 + b1 * pred$sfc_temp_mean[1] +          # predicted log10 pepspp based on predicted temperature 
                              randvar_pepspp[1]                         #   plus random variation 

  sim$sfc_temp_mean[(i-1)*30+1] <- pred$sfc_temp_mean[1]                # surface temperature prediction stored in sim table
  sim$log_pepspp_mean[(i-1)*30+1] <- pred$log_pepspp_mean[1]            # butterfish prediction stored in sim table

# years 2 and onward predictions  
  for (j in 2:num_pred) {                                               # repeat computations for years 2-30
    sim$sim[(i-1)*30+j] <- i                                            # track simulation number
    pred$sfc_temp_mean[j] <- pred$sfc_temp_mean[j-1] + randvar_temp[j]  # predicted surface temperature
    pred$log_pepspp_mean[j] <- b0 + b1 * pred$sfc_temp_mean[j] +        # predicted log10 pepspp  
                                randvar_pepspp[j]                       
    
    sim$sfc_temp_mean[(i-1)*30+j] <- pred$sfc_temp_mean[j]              # store predicted values in the sim table
    sim$log_pepspp_mean[(i-1)*30+j] <- pred$log_pepspp_mean[j]
    }
}

# Compute mean and sd of predicted values by year
sim_temp_mean <- (with(sim, tapply(sfc_temp_mean,list(year=year),mean)))          # mean temperature
sim_temp_sd <- (with(sim, tapply(sfc_temp_mean,list(year=year),sd)))              # sd temperature
sim_pepspp_mean <- (with(sim, tapply(log_pepspp_mean,list(year=year),mean)))      # mean log10 butterfish CPUE
sim_pepspp_sd <- (with(sim, tapply(log_pepspp_mean,list(year=year),sd)))          # sd log10 butterfish CPUE

# Check that simulated values correctly mimic the assumed level of variability
mean(sim_temp_sd)   # > mean(sim_temp_sd)
                    # [1] 1.094438  (However the first year value is ~0.3, similar to sd of random numbers)

mean(sim_pepspp_sd)   # ~ mean(sim_pepspp_sd)
                      # [1] 0.02798031  (similar to sd of random numbers, which is ~0.03)

# Plot simulation results
# First, combine simulation results into a single table
temp_plus_sd <- sim_temp_mean + sim_temp_sd                       # add columns for +- one SD, surface temperature
temp_minus_sd <- sim_temp_mean - sim_temp_sd 
pepspp_plus_sd <- sim_pepspp_mean + sim_pepspp_sd                 # add columns for +- one SD, log10 butterfish CPUE
pepspp_minus_sd <- sim_pepspp_mean - sim_pepspp_sd 
sim_table <- data.frame(cbind(pred$year,sim_temp_mean, temp_plus_sd, temp_minus_sd,
                              sim_pepspp_mean, pepspp_plus_sd, pepspp_minus_sd))
colnames(sim_table) <- c("year","temp.mean","temp.plus.sd","temp.minus.sd","pepspp.mean","pepspp.plus.sd","pepspp.minus.sd")

# Use ggplot to display results (more flexible)
#   Plot mean and sd of predicted surface temperature
theme_set(theme_bw())  # pre-set the bw theme.
p1 <- ggplot(sim_table, aes(year, temp.mean)) + 
  geom_line(aes(y = temp.mean)) +
  geom_line(aes(y = temp.plus.sd), linetype="dashed") +
  geom_line(aes(y = temp.minus.sd), linetype="dashed") +
  labs(y = "Surface temperature (deg C)", x = "Year", title = "Simulations") +
  theme(
    axis.text.x = element_text(color="black", size=12, face="plain"),
    axis.text.y = element_text(color="black", size=12, face="plain"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold")
  )
p1

# Plot mean and sd of predicted log10 butterfish CPUE
p2 <- ggplot(sim_table, aes(year, pepspp.mean)) + 
  geom_line(aes(y = pepspp.mean)) +
  geom_line(aes(y = pepspp.plus.sd), linetype="dashed") +
  geom_line(aes(y = pepspp.minus.sd), linetype="dashed") +
  labs(y = "log10 butterfish CPUE", x = "Year", title = "Simulations") +
  theme(
    axis.text.x = element_text(color="black", size=12, face="plain"),
    axis.text.y = element_text(color="black", size=12, face="plain"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold")
  )
p2

# Plot one simulation, including observations
# First select columns and then combine observations with first simulation
# surface temperature
my_data_temp <- subset(my_data,select = c(year,sfc_temp_mean))       # select year and surface temperature from data table
colnames(my_data_temp) <- c("year","sfc_temp_mean")                  # label columns
one_sim_temp <- subset(sim, sim==1)                                  # select first simulation
one_sim_temp <- subset(one_sim_temp,select = c(year,sfc_temp_mean))  # select year and surface temperature
one_sim_temp <- rbind(my_data_temp,one_sim_temp)                     # combine observations and first simulation year

# log10 butterfish CPUE
my_data_pepspp <- subset(my_data,select = c(year,log_pepspp_mean))         # select year and log10 butterfish CPUE from data table
colnames(my_data_pepspp) <- c("year","log_pepspp_mean")                    # label columns
one_sim_pepspp <- subset(sim, sim==1)                                      # select first simulation
one_sim_pepspp <- subset(one_sim_pepspp,select = c(year,log_pepspp_mean))  # select year and log10 butterfish CPUE
one_sim_pepspp <- rbind(my_data_pepspp,one_sim_pepspp)                     # combine observatins and first simulation year

# Plot surface temperature
# Observations are much more variable than the predictions
# Likely an effect of the temporal and spatial variability of the observations
theme_set(theme_bw())  # pre-set the bw theme.
p3 <- ggplot(one_sim_temp, aes(year, sfc_temp_mean)) + 
  geom_line(aes(y = sfc_temp_mean)) +
  labs(y = "Surface temperature (deg C)", x = "Year", title = "Observed and one simulation") +
  theme(
    axis.text.x = element_text(color="black", size=12, face="plain"),
    axis.text.y = element_text(color="black", size=12, face="plain"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold")
  )
p3

# Plot log10 butterfish CPUE
theme_set(theme_bw())  # pre-set the bw theme.
p4 <- ggplot(one_sim_pepspp, aes(year, log_pepspp_mean)) + 
  geom_line(aes(y = log_pepspp_mean)) +
  labs(y = "log10 butterfish CPUE", x = "Year", title = "Observed and one simulation") +
  theme(
    axis.text.x = element_text(color="black", size=12, face="plain"),
    axis.text.y = element_text(color="black", size=12, face="plain"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold")
  )
p4

# Plot observations and all simulation results together
# First select columns and then combine observations with simulation average
my_data_temp <- subset(my_data,select = c(year,sfc_temp_mean))     # select year and surface temperature from observations
colnames(my_data_temp) <- c("year","temp.mean")                    # label columns
my_data_pepspp <- subset(my_data,select = c(log_pepspp_mean))      # select year and log10 butterfish CPUE from observations
colnames(my_data_pepspp) <- c("pepspp.mean")                       # label columns

temp1 <- data.frame(matrix(nrow=18,ncol=2,byrow=TRUE))             # create placeholder values for sd temperature (all NA values)
                                                                   # Need these placeholder columns for matching format of sim_table (below)
colnames(temp1) <- c("temp.plus.sd","temp.minus.sd")               # label columns

temp2 <- data.frame(matrix(nrow=18,ncol=2,byrow=TRUE))             # create placeholder values for sd butterfish CPUE (all NA values)
colnames(temp2) <- c("pepspp.plus.sd","pepspp.minus.sd")           # label columns

my_data_temp_pepspp <- cbind(my_data_temp,temp1,my_data_pepspp,temp2) # combine observation tables
obs_sim <- rbind(my_data_temp_pepspp,sim_table)                       # combine observations with simulations

# Plot temperature predictions
theme_set(theme_bw())  # pre-set the bw theme.
p5 <- ggplot(obs_sim, aes(year, temp.mean)) + 
  geom_line(aes(y = temp.mean)) +
  geom_line(aes(y = temp.plus.sd), linetype="dashed") +
  geom_line(aes(y = temp.minus.sd), linetype="dashed") +
  labs(y = "Surface temperature (deg C)", x = "Year", title = "Observed and all simulations") +
  theme(
    axis.text.x = element_text(color="black", size=12, face="plain"),
    axis.text.y = element_text(color="black", size=12, face="plain"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold")
  )
p5
# FYI, the warning messages are due to the placeholder (mising) values for sd of the observations

# Plot lgo10 butterfish CPUE predictions
theme_set(theme_bw())  # pre-set the bw theme.
p6 <- ggplot(obs_sim, aes(year, pepspp.mean)) + 
  geom_line(aes(y = pepspp.mean)) +
  geom_line(aes(y = pepspp.plus.sd), linetype="dashed") +
  geom_line(aes(y = pepspp.minus.sd), linetype="dashed") +
  labs(y = "log10 butterfish CPUE", x = "Year", title = "Observed and all simulations") +
  theme(
    axis.text.x = element_text(color="black", size=12, face="plain"),
    axis.text.y = element_text(color="black", size=12, face="plain"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold")
  )
p6


# Plot everything together!
pdf("Statistical forecasting 3a output.pdf")
grid.arrange(p3, p4, p5, p6, ncol = 2)
dev.off()

##########
# OPTIONAL
####################################
####################################
####################################
#PRINCIPAL COMPONENTS ANALYSIS (PCA)
####################################
####################################
####################################
install.packages("dplyr")
library(dplyr)

##################################
#1. Read data and cross tabulate.#
##################################
Catch <- as.data.frame(read.csv("Plankton.csv",header=T))   # The output file "Plankton.csv" has an added first column (record number)
Catch$date <- as.Date(Catch$date, "%Y-%m-%d")		            # Assign Date format 'yyyy-mm-dd'
xtabs(~ year + month, data = Catch)                         # Crosstabulate the data by year and month

#####################################
# 2. compute yearly means by species#
#####################################
Catch.year <- Catch[,c(23,26:31)]                           # Select columns with log(ichthyoplankton values), include year
Catch.byyear <- Catch.year %>%                              # Compute average CPUE by year and species
  group_by(year) %>%
  summarise_all(list(mean=mean), na.rm = TRUE)

####################################################################
# 3. PCA                                                           #
####################################################################
pairs(Catch.byyear[,c(2:7)],lower.panel=NULL)               # Plot data pairs
PCA <- prcomp(Catch.byyear[,c(2:7)], scale=TRUE)            # Run PCA
summary(PCA)                                                # PCA summary
biplot(PCA, scale=0)                                        # Plot first two principal components

pdf("PCA.pdf")
par(mfrow = c(1,1))
biplot(PCA, scale=0)                                        # Plot first two principal components
chart.Correlation(my_data, histogram=TRUE, pch=19)          # Correlation coefficients
dev.off()

############
# END OF LAB
############

####################################
####################################
####################################
# PREPARE FILE FOR STATISTICAL FORECASTING
####################################
####################################
####################################
install.packages("dplyr")
install.packages("PerformanceAnalytics")
library("dplyr")
library("PerformanceAnalytics")

##################################
#Read data and cross tabulate.#
##################################
Catch <- as.data.frame(read.csv("Plankton.csv",header=T))   # The file "Plankton.csv" has an added first column (record number)
Catch$date <- as.Date(Catch$date, "%Y-%m-%d")		            # Assign Date format 'yyyy-mm-dd'
xtabs(~ year + month, data = Catch)                         # Crosstabulate the data by year and month

Tern <- as.data.frame(read.csv("Tern productivity.csv",header=T))   # Tern productivity

##########################
#Compute yearly means#
##########################
Catch.year <- Catch[,c(23,26:31)]                           # Select columns with log(ichthyoplankton values), include year

# Compute average CPUE by year and species
Catch.byyear <- Catch.year %>%                              
  group_by(year) %>%
  summarise_all(list(mean=mean), na.rm = TRUE)

SprTemp.year <- subset(Catch,month>=4 &month<=6)            # Select spring only (April-June)
xtabs(~ year + month, data = SprTemp.year)                  # Crosstabulate the data by year and month, to check subset
SprTemp.year <- subset(SprTemp.year[,c(9,23)])              # Select surface temperature, include year

# Compute average spring temperature by year
SprTemp.byyear <- SprTemp.year %>%                          
  group_by(year) %>%
  summarise_all(list(mean=mean), na.rm = TRUE)
colnames(SprTemp.byyear) <- c("year","sfc_temp_mean")       # Change column name from "mean" to "sfc_temp_mean" to be consistent with ichthyoplankton naming

# Merge files by year
my_data <- merge(Catch.byyear,Tern, by.x="year", by.y="Year",all.x=TRUE)
my_data <- subset(my_data[,c(1:7,10)])                      # Keep only common tern number fledged per pair
my_data <- merge(my_data,SprTemp.byyear, by="year",all.x=TRUE)
my_data <- subset(my_data,year>=1999)

my_data <- merge(SprTemp.byyear, Catch.byyear, by="year",all.x=TRUE)
my_data <- merge(my_data,Tern, by.x="year", by.y="Year",all.x=TRUE)
my_data <- subset(my_data[,c(1:8,11)])                      # Keep only common tern number fledged per pair
my_data <- subset(my_data,year>=1999)

pairs(my_data[,c(2:9)],upper.panel=NULL)
chart.Correlation(my_data, histogram=TRUE, pch=19)

# Write output file (which will be used as input to lab for statistical forecasting)
write.csv(my_data, file = "TempPlanktonCommontern.csv")

####################################
####################################
####################################
# PREPARE FILE FOR CHRIS'S LECTURE AND LAB
####################################
####################################
####################################
install.packages("dplyr")
install.packages("PerformanceAnalytics")
library("dplyr")
library("PerformanceAnalytics")

##################################
#Read data and cross tabulate.#
##################################
Catch <- as.data.frame(read.csv("Plankton.csv",header=T))   # The file "Plankton.csv" has an added first column (record number)
Catch$date <- as.Date(Catch$date, "%Y-%m-%d")		            # Assign Date format 'yyyy-mm-dd'
xtabs(~ year + month, data = Catch)                         # Crosstabulate the data by year and month

Tern <- as.data.frame(read.csv("Tern productivity.csv",header=T))   # Tern productivity

##########################
#Compute yearly means#
##########################
Catch.year <- Catch[,c(23,26:31)]                           # Select columns with log(ichthyoplankton values), include year

# Compute average CPUE by year and species
Catch.byyear <- Catch.year %>%                              
  group_by(year) %>%
  summarise_all(list(mean=mean), na.rm = TRUE)

SprTemp.year <- subset(Catch,month>=4 &month<=6)            # Select spring only (April-June)
xtabs(~ year + month, data = SprTemp.year)                  # Crosstabulate the data by year and month, to check subset
SprTemp.year <- subset(SprTemp.year[,c(9,23)])              # Select surface temperature, include year

# Compute average spring temperature by year
SprTemp.byyear <- SprTemp.year %>%                          
  group_by(year) %>%
  summarise_all(list(mean=mean), na.rm = TRUE)
colnames(SprTemp.byyear) <- c("year","sfc_temp_mean")       # Change column name from "mean" to "sfc_temp_mean" to be consistent with ichthyoplankton naming

# Merge files by year
my_data <- merge(Catch.byyear,Tern, by.x="year", by.y="Year",all.x=TRUE)
my_data <- subset(my_data[,c(1:7,10,14)])                      # Keep only number fledged per pair for both common and roseate tern
my_data <- merge(my_data,SprTemp.byyear, by="year",all.x=TRUE)
my_data <- subset(my_data,year>=1999)

my_data <- merge(SprTemp.byyear, Catch.byyear, by="year",all.x=TRUE)
my_data <- merge(my_data,Tern, by.x="year", by.y="Year",all.x=TRUE)
my_data <- subset(my_data[,c(1:8,11,15)])                      # Keep only number fledged per pair for common tern and roseate tern
my_data <- subset(my_data,year>=1999)

pairs(my_data[,c(2:9)],upper.panel=NULL)
chart.Correlation(my_data, histogram=TRUE, pch=19)

# Write output file (which will be used as input to Chris's lab exercises)
write.csv(my_data, file = "TempPlanktonCommonRoseatetern.csv")


##################################
##################################
##################################
# Experiment with GENERAL ADDITIVE MODEL (GAM)
##################################
##################################
##################################
library(fields)
library(dplyr)
library(mgcv)
library(mapdata)
library(mapproj)

# Appledore Island location
Lat.SML <- 42.987727
Lon.SML <- -70.613940

##################################
#1. Read data and cross tabulate.#
##################################
Catch <- as.data.frame(read.csv("Plankton.csv",header=T))   # The output file "Plankton.csv" has an added first column (record number)
Catch$date <- as.Date(Catch$date, "%Y-%m-%d")		            # Assign Date format 'yyyy-mm-dd'
xtabs(~ year + month, data = Catch)                         # Crosstabulate the data by year and month

start.year <- 1989  # all years, can change start year if desired (e.g., Ecomon years: 1989-present)
Catch.ich <- subset(Catch, year>=start.year)         

#1.5 Estimate ichthyoplankton abundance estimates 
Catch.ich <- subset(Catch.ich, distance < 150)   # select samples within 150 km of Appledore Island, i.e., the Gulf of Maine

# TESTING: no spatial component
# For each species, fit GAM as a function of year and month; summarize results and plot diagnostics (gam.check)
# year is a smoothed effect
par(mfrow = c(2,2))   # format plot layout (2 rows and 2 columns of plots per page)
fit.cluhar <- gam(log_cluhar ~ s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.cluhar); gam.check(fit.cluhar)
fit.urospp <- gam(log_urospp ~ s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.urospp); gam.check(fit.urospp)
fit.merbil <- gam(log_merbil ~ s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.merbil); gam.check(fit.merbil)
fit.scosco <- gam(log_scosco ~ s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.scosco); gam.check(fit.scosco)
fit.pepspp <- gam(log_pepspp ~ s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.pepspp); gam.check(fit.pepspp)
fit.ammspp <- gam(log_ammspp ~ s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.ammspp); gam.check(fit.ammspp)
par(mfrow = c(1,1))  # reset plot layout

# TESTING: year as random effect
par(mfrow = c(2,2))   # format plot layout (2 rows and 2 columns of plots per page)
fit.cluhar <- gam(log_cluhar ~ s(year,bs="re") + s(month,bs='cc'), data=Catch.ich); summary(fit.cluhar); gam.check(fit.cluhar) # nonsensical result
par(mfrow = c(1,1))  # reset plot layout

# Output plotted results to a file
pdf("GAM 1.5 output.pdf")
par(mfrow = c(2,3))

newd <- data.frame(year=(1989:2016),month=rep.int(6,2016-1989+1))
year <- (1989:2016)
pred.cluhar <- predict.gam(fit.cluhar,newd)
pred.urospp <- predict.gam(fit.urospp,newd)
pred.merbil <- predict.gam(fit.merbil,newd)
pred.scosco <- predict.gam(fit.scosco,newd)
pred.pepspp <- predict.gam(fit.pepspp,newd)
pred.ammspp <- predict.gam(fit.ammspp,newd)

plot(year,pred.cluhar)
plot(year,pred.urospp)
plot(year,pred.merbil)
plot(year,pred.scosco)
plot(year,pred.pepspp)
plot(year,pred.ammspp)

#by year
plot(fit.cluhar,pages=0,seWithMean=TRUE,select=1,ylim=c(-0.3,0.3),main="herring")
plot(fit.urospp,pages=0,seWithMean=TRUE,select=1,ylim=c(-0.2,0.2),main="red/white hake")
plot(fit.merbil,pages=0,seWithMean=TRUE,select=1,ylim=c(-0.1,0.1),main="silver hake")
plot(fit.scosco,pages=0,seWithMean=TRUE,select=1,ylim=c(-0.1,0.1),main="mackerel")
plot(fit.pepspp,pages=0,seWithMean=TRUE,select=1,ylim=c(-0.1,0.1),main="butterfish")
plot(fit.ammspp,pages=0,seWithMean=TRUE,select=1,ylim=c(-0.4,0.4),main="sandlance")

#by month
plot(fit.cluhar,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.3,0.3),main="herring")
plot(fit.urospp,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.4,0.4),main="red/white hake")
plot(fit.merbil,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.3,0.3),main="silver hake")
plot(fit.scosco,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.3,0.3),main="mackerel")
plot(fit.pepspp,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.3,0.3),main="butterfish")
plot(fit.ammspp,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.6,0.6),main="sandlance")

dev.off()

###############################
# Modify code from Ingrid
library("PerformanceAnalytics")

pdf("Yearly_predictions.pdf")
par(mfrow = c(6,2))

year <- seq(1989,2016)
pred.allspp <- data.frame(cbind(year,pred.ammspp,pred.cluhar,pred.merbil,pred.pepspp,pred.scosco,pred.urospp))
colnames(pred.allspp) <- c("Year","sandlance","herring","silver hake","butterfish","mackerel","red/white hake")
my_data <- pred.allspp[,c(2:7)]                           # Select columns with ichthyoplankton values
str(pred.allspp)

chart.Correlation(my_data, histogram=TRUE, pch=19)

# Add tern data
Tern <- as.data.frame(read.csv("Tern productivity.csv",header=T))   # Tern productivity
Tern.years <- subset(Tern,Year>=1999&Year<=2016)
pred.allspp <- subset(pred.allspp,Year>=1999&Year<=2016)
Tern.pred <- cbind(Tern.years,pred.allspp)
str(Tern.pred)

my_data <- Tern.pred[,c("COTE_numberfledged","sandlance","herring","silver hake","butterfish","mackerel","red/white hake")]
str(my_data)

chart.Correlation(my_data, histogram=TRUE, pch=19)

# Code sent by Ingrid for me to test, on 7/5/2020
pdf("Yearly Predictions .pdf")
par(mfrow = c(6,2))

# This line misidentifies year as character, which causes correlation coefficient to fail
# year <- c("1989","1990","1991","1992","1993","1994","1995","1996","1997","1998","1999","2000","2001","2002","2003",
#          "2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016")
# Instead use
year <- seq(1989,2016)

pred.allspp <- data.frame(cbind(year, pred.ammspp,pred.cluhar,pred.merbil,pred.pepspp,pred.scosco,pred.urospp))
colnames(pred.allspp) <- c("year", "sandlance","herring","silver hake","butterfish","mackeral","red/white hake")

pred.allspp

my_data <- pred.allspp[c(2:7)]  #selecting columns with ichthyoplankton data
str(pred.allspp)

chart.Correlation(my_data, histogram=TRUE, pch=19)

dev.off()

#############################################################################################################
# TESTING: include spatial component, draw polygons on spatial output
# From Chris Rooper, 6/22/2020: If you’re doing it in R, a way to do it would be to
# 1)     Make a spatial point for the center of your circle
# 2)     Buffer the point with a radius of 150 km (turning it into a spatial polygon)
# 3)     Use the raster::extract function to pull a mean from the overlapping polygon and CPUE rasters
# Items 1 and 2 can be done using the sp package I think. I probably have some code that would get you close.

library(fields)
library(dplyr)
library(mgcv)
library(mapdata)
library(mapproj)
library(rgeos)
library(sp)
library(raster)
library(rgdal)

# Appledore Island location
Lat.SML <- 42.987727
Lon.SML <- -70.613940

# Code from Chris Rooper
Points<-SpatialPoints(cbind(-122.5,58.5),proj4string=CRS("+proj=longlat +datum=WGS84"))
points.project <- spTransform(Points, "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-80 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 
                              +units=m +no_defs")

#Bufferedpoint<-gBuffer(points.project,byid=TRUE,width=150,capStyle=”ROUND”) # 150 m buffer
Bufferedpoint<-gBuffer(points.project,byid=TRUE,width=150) # 150 m buffer

Meanvalue<-raster::extract(yourraster,Bufferedpoint,fun=”mean”) # Then you take your raster layer and extract the mean value inside the buffered point

# My adapted code
Points<-SpatialPoints(cbind(Lon.SML,Lat.SML),
                      proj4string=CRS("+proj=longlat +datum=WGS84")) # create object of class "SpatialPoints-class" (irregularly spaced points)
                                                                     # CRS = coordinate reference system
points.project <- spTransform(Points, "+proj=aea +lat_1=45 +lat_2=55 +lat_0=40 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 
                              +units=m +no_defs")                    # Transforms "Points" to an albers equal area projection.
                                                                     #   This will allow you to use a 150 m buffer.
Bufferedpoint<-gBuffer(points.project,byid=TRUE,width=150)           # 150 m buffer

Meanvalue<-raster::extract(raster.gam,Bufferedpoint,fun=mean)      # Then you take your raster layer and extract the mean value inside the buffered point


##########################################
###TO DO NEXT: figure out how to convert the spatial gam predictions to a raster, presumably using the same projection
# Create raster layer with a prediction layer basedon a fitted model object
# See https://www.rdocumentation.org/packages/raster/versions/1.7-18/topics/predict
# and https://rspatial.org/raster/spatial/8-rastermanip.html

# Output from a simple, spatial GAM (see below)
vis.gam(fit.cluhar, plot.type="contour", too.far=0.15, main="Atlantic herring", xlab="Longitude",ylab="Latitude",color="topo")
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")

raster.gam<-vis.gam(fit.cluhar, plot.type="contour", too.far=0.15, main="Atlantic herring", xlab="Longitude",ylab="Latitude",color="topo")
raster.gam
#NULL
raster.gam<-predict(fit.cluhar)
str(raster.gam)

# GAM by year, month, and location (all years)
##################################
#1. Read data and cross tabulate.#
##################################
Catch <- as.data.frame(read.csv("Plankton.csv",header=T))   # The output file "Plankton.csv" has an added first column (record number)
Catch$date <- as.Date(Catch$date, "%Y-%m-%d")		            # Assign Date format 'yyyy-mm-dd'
xtabs(~ year + month, data = Catch)                         # Crosstabulate the data by year and month

start.year <- 1989  # all years, can change start year if desired (e.g., Ecomon years: 1989-present)
Catch.ich <- subset(Catch, year>=start.year)         

####################################################################################
#2. Fit GAM to log-transformed values, log10(value+1), include year and month terms#
# One spatial map for each species covering all years (i.e., average location)     #
#  run time about 1 minute                                                         #
####################################################################################

start.year <- 1989  # all years, can change start year if desired (e.g., Ecomon years: 1989-present)
Catch.ich <- subset(Catch, year>=start.year)         

# For each species, fit GAM as a function of location, year, and month; summarize results and plot diagnostics (gam.check)
par(mfrow = c(2,2))   # format plot layout (2 rows and 2 columns of plots per page)
fit.cluhar <- gam(log_cluhar ~ s(lon, lat) + s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.cluhar); gam.check(fit.cluhar)
fit.urospp <- gam(log_urospp ~ s(lon, lat) + s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.urospp); gam.check(fit.urospp)
fit.merbil <- gam(log_merbil ~ s(lon, lat) + s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.merbil); gam.check(fit.merbil)
fit.scosco <- gam(log_scosco ~ s(lon, lat) + s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.scosco); gam.check(fit.scosco)
fit.pepspp <- gam(log_pepspp ~ s(lon, lat) + s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.pepspp); gam.check(fit.pepspp)
fit.ammspp <- gam(log_ammspp ~ s(lon, lat) + s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.ammspp); gam.check(fit.ammspp)
par(mfrow = c(1,1))  # reset plot layout

# Output plotted results to a file
pdf("GAM 2 output.pdf")
par(mfrow = c(2,3))

#by year
plot(fit.cluhar,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.3,0.3),main="herring")
plot(fit.urospp,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.2,0.2),main="red/white hake")
plot(fit.merbil,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.1,0.1),main="silver hake")
plot(fit.scosco,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.1,0.1),main="mackerel")
plot(fit.pepspp,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.1,0.1),main="butterfish")
plot(fit.ammspp,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.4,0.4),main="sandlance")

#by month
plot(fit.cluhar,pages=0,seWithMean=TRUE,select=3,ylim=c(-0.3,0.3),main="herring")
plot(fit.urospp,pages=0,seWithMean=TRUE,select=3,ylim=c(-0.4,0.4),main="red/white hake")
plot(fit.merbil,pages=0,seWithMean=TRUE,select=3,ylim=c(-0.3,0.3),main="silver hake")
plot(fit.scosco,pages=0,seWithMean=TRUE,select=3,ylim=c(-0.3,0.3),main="mackerel")
plot(fit.pepspp,pages=0,seWithMean=TRUE,select=3,ylim=c(-0.3,0.3),main="butterfish")
plot(fit.ammspp,pages=0,seWithMean=TRUE,select=3,ylim=c(-0.6,0.6),main="sandlance")

# spatial plots
vis.gam(fit.cluhar, plot.type="contour", too.far=0.15, main="Atlantic herring", xlab="Longitude",ylab="Latitude",color="topo")
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
vis.gam(fit.urospp, plot.type="contour", too.far=0.15, main="red white hake", xlab="Longitude",ylab="Latitude",color="topo")
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
vis.gam(fit.merbil, plot.type="contour", too.far=0.15, main="silver hake", xlab="Longitude",ylab="Latitude",color="topo")
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
vis.gam(fit.scosco, plot.type="contour", too.far=0.15, main="Atlantic mackerel", xlab="Longitude",ylab="Latitude",color="topo")
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
vis.gam(fit.pepspp, plot.type="contour", too.far=0.15, main="butterfish", xlab="Longitude",ylab="Latitude",color="topo")
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
vis.gam(fit.ammspp, plot.type="contour", too.far=0.15, main="sandlance", xlab="Longitude",ylab="Latitude",color="topo")
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")

par(mfrow = c(1,1))
dev.off()


>>>>>>> 15f0978613261279179a1de444f7bae29b9fb5e1
