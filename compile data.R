# Program to compile Gulf of Maine plankton data (EcoMon surveys)
# In some years, sampling was monthly, but in recent years, sampling focused on April-June and August-November
# Michael Sigler
# May 19, 2019
# Revised November 22, 2019, July 19, 2020 (to add data through 2017)
# Data source: https://www.ncei.noaa.gov/metadata/geoportal/rest/metadata/item/gov.noaa.nodc%3A0187513/html

options(scipen=999)  # turn-off scientific notation like 1e+48
library(fields)
library(dplyr)
library(mgcv)
library(mapdata)
library(mapproj)

############################################################################################################
# Read plankton file 
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
# merbil_100m3 = silver hake count per 100 m3 of volume, *** drop from analysis and substitute polvir (pollock) to match tern diet database
# cluhar_100m3 = Atlantic herring count per 100 m3 of volume
# ammspp_100m3 = sandlance count per 100 m3 of volume
# pepspp_100m3 = butterfish count per 100 m3 of volume
# scosco_100m3 = Atlantic mackerel count per 100 m3 of volume
# urospp_100m3 = red and white count per 100 m3 of volume

Catch <- as.data.frame(read.csv("EcoMon 2017.csv",header=T))
Catch$date <- as.Date(Catch$date, "%Y-%m-%d")		# State Date format 'yyyy-mm-dd'
Catch$year <- lubridate::year(Catch$date)				# Assign year
Catch$month <- lubridate::month(Catch$date)			# Assign month
xtabs(~ year + month, data = Catch)

# Distance calculation, longitude then latitude
# Declare matrices for intermediate calculations
y1 <- data.frame(matrix(nrow=1,ncol=2,byrow=TRUE));y1
y2 <- data.frame(matrix(nrow=1,ncol=2,byrow=TRUE));y2

# Appledore Island location
Lat.SML <- 42.987727
Lon.SML <- -70.613940

# Compute great circle distance from each station to Appledore Island
nsamp <- length(Catch$year)											# 31,351 samples
dist <- data.frame(matrix(nrow=nsamp,ncol=1)); colnames(dist) = "distance"
for (j in 1:nsamp)
{
  y1[1:2] <- c(Lon.SML,Lat.SML)   									# longitude then latitude
  y2[1:2] <- c(Catch$lon[j],Catch$lat[j])   							# longitude then latitude
  dist$distance[j] <- rdist.earth(y1[1:2], y2[1:2], miles = FALSE, R = NULL)		
}
Catch <- cbind(Catch,dist)

# Tally occurrence of zooplankton and ichthyoplankton by area and volume
#   Sometimes only one taxa group is tallied for a cruise
#   1 = tallied, 0 = not tallied
sum.catch <- data.frame(matrix(data=NA,nrow=nsamp,ncol=8))
colnames(sum.catch) = c("zooparea","zoopvol","ichthyoarea","ichthyovol",
                        "zooparea.present","zoopvol.present",
                        "ichthyoarea.present","ichthyovol.present")

for (j in 1:nsamp)
{
# Sum catches
  sum.catch[j,1] <- Catch[j,14]   # already summed in input file
  sum.catch[j,2] <- Catch[j,106]  # already summed in input file
  sum.catch[j,3] <- sum(Catch[j,199:243],na.rm=F)   									
  sum.catch[j,4] <- sum(Catch[j,245:289],na.rm=F)   

# Tally whether catch was processed (or not)  
  if(is.na(sum.catch[j,1])) {sum.catch[j,5] <- 0} else {sum.catch[j,5] <- 1}
  if(is.na(sum.catch[j,2])) {sum.catch[j,6] <- 0} else {sum.catch[j,6] <- 1}
  if(is.na(sum.catch[j,3])) {sum.catch[j,7] <- 0} else {sum.catch[j,7] <- 1}
  if(is.na(sum.catch[j,4])) {sum.catch[j,8] <- 0} else {sum.catch[j,8] <- 1}
}

Catch <- cbind(Catch,sum.catch)

# Check that unsampled (NA) rows are assigned '0' for *.present columns
Catch[23912,]   # zooparea.present = 0
Catch[27511,]   # zooparea.present = 0
Catch[27521,]   # zooparea.present = 0
Catch[27539,]   # ichthyoarea.present = 0
Catch[27547,]   # ichthyoarea.present = 0
Catch[28760,]   # ichthyoarea.present = 0

# Compare zooplankton and ichthyoplankton counts
xtabs(~ year, data = subset(Catch,zoopvol.present==1))
xtabs(~ year, data = subset(Catch,ichthyovol.present==1))
# for example, counts for 2016 differ: 
#   In input file, one NA for zooplankton and 6 more samples for zooplankton
#   i.e., looks okay 

# Write tern prey, total ichthyo, total zoop to file
Catch.rpt <- subset.data.frame(Catch,select = c(cruise_name,station,
                                                lat,lon,date,time,depth,
                                                sfc_temp,sfc_salt,
                                                btm_temp,btm_salt,
                                                zoopvol,zoopvol.present,
                                                ichthyovol,ichthyovol.present,
                                                cluhar_100m3,urospp_100m3,
                                                polvir_100m3,scosco_100m3,
                                                pepspp_100m3,ammspp_100m3,
                                                year,month,distance))

# Log-transform ichthyoplankton catch data
Catch.rpt$log_cluhar <- log10(Catch.rpt$cluhar_100m3+1)
Catch.rpt$log_urospp <- log10(Catch.rpt$urospp_100m3+1)
Catch.rpt$log_polvir <- log10(Catch.rpt$polvir_100m3+1)
Catch.rpt$log_scosco <- log10(Catch.rpt$scosco_100m3+1)
Catch.rpt$log_pepspp <- log10(Catch.rpt$pepspp_100m3+1)
Catch.rpt$log_ammspp <- log10(Catch.rpt$ammspp_100m3+1)

# Select only records where ichthyoplankton volume was measured (27,917 records)
Catch.rpt <- subset(Catch.rpt,ichthyovol.present==1)
nrow(Catch.rpt)
str(Catch.rpt)

write.csv(Catch.rpt, file = "Plankton.csv")

# Data compilation complete

############################################################################
# THIS HAS BEEN COMPLETED FOR DATA THROUGH 2015 BUT NOT THROUGH 2017 
# Start data analysis
#   Follow methods list in "Some ideas for introducing statistical methods"
# 1. Use log-transformed ichthyoplankton catch data to compute annual averages (uncorrected by sampling month)

# Select rows with tallied ichthyoplankton and columns with log(ichthyoplankton values), include year
Catch.log <- Catch.rpt[,c(22,25:30)]

# dplyr uses the operator %>% to denote taking what is on the left and putting it into the function on the right
Catch.log.byyear <- Catch.log %>% 
  group_by(year) %>%
  summarise_all(list(mean=mean,sd=sd), na.rm = TRUE)

par(mfrow=c(2,3))
plot(log_cluhar_mean ~ year,data=Catch.log.byyear,ylab="herring")
plot(log_urospp_mean ~ year,data=Catch.log.byyear,ylab="red/white hake")
plot(log_merbil_mean ~ year,data=Catch.log.byyear,ylab="silver hake")
plot(log_scosco_mean ~ year,data=Catch.log.byyear,ylab="mackerel")
plot(log_pepspp_mean ~ year,data=Catch.log.byyear,ylab="butterfish")
plot(log_ammspp_mean ~ year,data=Catch.log.byyear,ylab="sandlance")
par(mfrow=c(1,1))

# 2. Limit region 
# Select rows within 150 km of Appledore Island
Catch.SML.log <- subset(Catch.rpt,distance<=150)
Catch.SML.log <- Catch.SML.log[,c(22,25:30)]

# dplyr uses the operator %>% to denote taking what is on the left 
#  and putting it into the function on the right
Catch.SML.log.byyear <- Catch.SML.log %>% 
  group_by(year) %>%
  summarise_all(list(mean=mean,sd=sd), na.rm = TRUE)

par(mfrow=c(2,3))
plot(log_cluhar_mean ~ year,data=Catch.SML.log.byyear,ylab="herring")
plot(log_urospp_mean ~ year,data=Catch.SML.log.byyear,ylab="red/white hake")
plot(log_merbil_mean ~ year,data=Catch.SML.log.byyear,ylab="silver hake")
plot(log_scosco_mean ~ year,data=Catch.SML.log.byyear,ylab="mackerel")
plot(log_pepspp_mean ~ year,data=Catch.SML.log.byyear,ylab="butterfish")
plot(log_ammspp_mean ~ year,data=Catch.SML.log.byyear,ylab="sandlance")
par(mfrow=c(1,1))

# Compute correlation coefficient for all regions versus SML values
Catch.SML.log.byyear$year # missing 1995-1998
x1<- subset(Catch.log.byyear,year<1995)
x2<- subset(Catch.log.byyear,year>1998)
Catch.log.byyear.subset <- rbind(x1,x2) #exclude 1995-1998 data

corr.byspecies <- data.frame(matrix(nrow=6,ncol=1))
colnames(corr.byspecies) <- "corr"
rownames(corr.byspecies) <- c("cluhar","urospp","merbil","scosco","pepspp","ammspp")

corr.byspecies[1,1] <- cor(Catch.log.byyear.subset$log_cluhar_mean,Catch.SML.log.byyear$log_cluhar_mean)
corr.byspecies[2,1] <- cor(Catch.log.byyear.subset$log_urospp_mean,Catch.SML.log.byyear$log_urospp_mean)
corr.byspecies[3,1] <- cor(Catch.log.byyear.subset$log_merbil_mean,Catch.SML.log.byyear$log_merbil_mean)
corr.byspecies[4,1] <- cor(Catch.log.byyear.subset$log_scosco_mean,Catch.SML.log.byyear$log_scosco_mean)
corr.byspecies[5,1] <- cor(Catch.log.byyear.subset$log_pepspp_mean,Catch.SML.log.byyear$log_pepspp_mean)
corr.byspecies[6,1] <- cor(Catch.log.byyear.subset$log_ammspp_mean,Catch.SML.log.byyear$log_ammspp_mean)

corr.byspecies

# 3. compute monthly global mean and then annual monthly residuals
# Select rows with log(ichthyoplankton values), include month
Catch.log.month <- Catch.rpt[,c(23,25:30)]

# Compute average catch by month
Catch.log.bymonth <- Catch.log.month %>% 
  group_by(month) %>%
  summarise_all(list(mean=mean), na.rm = TRUE)

# plot monthly averages
par(mfrow=c(2,3))
plot(Catch.log.bymonth$month,Catch.log.bymonth$log_cluhar_mean,ylab="herring",xlab="month")
plot(Catch.log.bymonth$month,Catch.log.bymonth$log_urospp_mean,ylab="red/white hake",xlab="month")
plot(Catch.log.bymonth$month,Catch.log.bymonth$log_merbil_mean,ylab="silver hake",xlab="month")
plot(Catch.log.bymonth$month,Catch.log.bymonth$log_scosco_mean,ylab="mackerel",xlab="month")
plot(Catch.log.bymonth$month,Catch.log.bymonth$log_pepspp_mean,ylab="butterfish",xlab="month")
plot(Catch.log.bymonth$month,Catch.log.bymonth$log_ammspp_mean,ylab="sandlance",xlab="month")
par(mfrow=c(1,1))

# use the merge command to merge monthly averages (Catch.log.bymonth) with Catch.rpt based on month matching
# compute residuals
# repeat computation of annual averages

Catch.resid<-merge(Catch.rpt,Catch.log.bymonth,by="month")

Catch.resid$cluhar_resid.log <- Catch.resid$log_cluhar-Catch.resid$log_cluhar_mean
Catch.resid$urospp_resid.log <- Catch.resid$log_urospp-Catch.resid$log_urospp_mean
Catch.resid$merbil_resid.log <- Catch.resid$log_merbil-Catch.resid$log_merbil_mean
Catch.resid$scosco_resid.log <- Catch.resid$log_scosco-Catch.resid$log_scosco_mean
Catch.resid$pepspp_resid.log <- Catch.resid$log_pepspp-Catch.resid$log_pepspp_mean
Catch.resid$ammspp_resid.log <- Catch.resid$log_ammspp-Catch.resid$log_ammspp_mean

# Select columns with  residual ichthyoplankton values, include year
Catch.resid.log.year <- Catch.resid[,c(23,37:42)]

hist(Catch.resid.log.year$cluhar_resid.log)   # Log transform distribution has left limb whereas 4th root does not

# Compute average residual by year
Catch.resid.log.byyear <- Catch.resid.log.year %>% 
  group_by(year) %>%
  summarise_all(list(mean=mean), na.rm = TRUE)

# plot yearly averages
par(mfrow=c(2,3))
plot(Catch.resid.log.byyear$year,Catch.resid.log.byyear$cluhar_resid.log_mean,ylab="herring",xlab="month")
plot(Catch.resid.log.byyear$year,Catch.resid.log.byyear$urospp_resid.log_mean,ylab="red/white hake",xlab="month")
plot(Catch.resid.log.byyear$year,Catch.resid.log.byyear$merbil_resid.log_mean,ylab="silver hake",xlab="month")
plot(Catch.resid.log.byyear$year,Catch.resid.log.byyear$scosco_resid.log_mean,ylab="mackerel",xlab="month")
plot(Catch.resid.log.byyear$year,Catch.resid.log.byyear$pepspp_resid.log_mean,ylab="butterfish",xlab="month")
plot(Catch.resid.log.byyear$year,Catch.resid.log.byyear$ammspp_resid.log_mean,ylab="sandlance",xlab="month")
par(mfrow=c(1,1))

# compare average log values to average residuals (also on a log scale) by correlation
corr.byspecies.resid <- data.frame(matrix(nrow=6,ncol=1))
colnames(corr.byspecies.resid) <- "corr"
rownames(corr.byspecies.resid) <- c("herring","red/white hake","silver hake","mackerel","butterfish","sandlance")

corr.byspecies.resid[1,1] <- cor(Catch.log.byyear$log_cluhar_mean,Catch.resid.log.byyear$cluhar_resid.log_mean)
corr.byspecies.resid[2,1] <- cor(Catch.log.byyear$log_urospp_mean,Catch.resid.log.byyear$urospp_resid.log_mean)
corr.byspecies.resid[3,1] <- cor(Catch.log.byyear$log_merbil_mean,Catch.resid.log.byyear$merbil_resid.log_mean)
corr.byspecies.resid[4,1] <- cor(Catch.log.byyear$log_scosco_mean,Catch.resid.log.byyear$scosco_resid.log_mean)
corr.byspecies.resid[5,1] <- cor(Catch.log.byyear$log_pepspp_mean,Catch.resid.log.byyear$pepspp_resid.log_mean)
corr.byspecies.resid[6,1] <- cor(Catch.log.byyear$log_ammspp_mean,Catch.resid.log.byyear$ammspp_resid.log_mean)

corr.byspecies.resid

# 4. compute monthly global mean and then annual monthly residuals, within 150 km of Shoals
# Start with the Plankton.csv file
Catch.rpt <- as.data.frame(read.csv("Plankton.csv",header=T))   # The output file "Plankton.csv" has an added first column (record number)
Catch.rpt$date <- as.Date(Catch.rpt$date, "%Y-%m-%d")		            # Assign Date format 'yyyy-mm-dd'
xtabs(~ year + month, data = Catch.rpt)                         # Crosstabulate the data by year and month

# Select rows within 150 km of Appledore Island
Catch.SML <- subset(Catch.rpt,distance<=150)
xtabs(~ year + month, data = Catch.SML)                         # Crosstabulate the data by year and month

# Select rows with log(ichthyoplankton values), include month
Catch.log.month <- Catch.SML[,c(23:24,26:31)]

# Compute average catch by month
Catch.log.bymonth <- Catch.log.month %>% 
  group_by(month) %>%
  summarise_all(list(mean=mean), na.rm = TRUE)

# plot monthly averages
par(mfrow=c(2,3))
plot(Catch.log.bymonth$month,Catch.log.bymonth$log_cluhar_mean,ylab="herring",xlab="month")
plot(Catch.log.bymonth$month,Catch.log.bymonth$log_urospp_mean,ylab="red/white hake",xlab="month")
plot(Catch.log.bymonth$month,Catch.log.bymonth$log_merbil_mean,ylab="silver hake",xlab="month")
plot(Catch.log.bymonth$month,Catch.log.bymonth$log_scosco_mean,ylab="mackerel",xlab="month")
plot(Catch.log.bymonth$month,Catch.log.bymonth$log_pepspp_mean,ylab="butterfish",xlab="month")
plot(Catch.log.bymonth$month,Catch.log.bymonth$log_ammspp_mean,ylab="sandlance",xlab="month")
par(mfrow=c(1,1))

# use the merge command to merge monthly averages (Catch.log.bymonth) with Catch.SML based on month matching
# compute residuals
# repeat computation of annual averages

Catch.resid<-merge(Catch.SML,Catch.log.bymonth,by="month")

Catch.resid$cluhar_resid.log <- Catch.resid$log_cluhar-Catch.resid$log_cluhar_mean
Catch.resid$urospp_resid.log <- Catch.resid$log_urospp-Catch.resid$log_urospp_mean
Catch.resid$merbil_resid.log <- Catch.resid$log_merbil-Catch.resid$log_merbil_mean
Catch.resid$scosco_resid.log <- Catch.resid$log_scosco-Catch.resid$log_scosco_mean
Catch.resid$pepspp_resid.log <- Catch.resid$log_pepspp-Catch.resid$log_pepspp_mean
Catch.resid$ammspp_resid.log <- Catch.resid$log_ammspp-Catch.resid$log_ammspp_mean

# Select columns with  residual ichthyoplankton values, include year
Catch.resid.log.year <- Catch.resid[,c(24,39:44)]

hist(Catch.resid.log.year$cluhar_resid.log)   # Log transform distribution has left limb whereas 4th root does not

# Compute average residual by year
Catch.resid.log.byyear <- Catch.resid.log.year %>% 
  group_by(year) %>%
  summarise_all(list(mean=mean), na.rm = TRUE)

# plot yearly averages
par(mfrow=c(2,3))
plot(Catch.resid.log.byyear$year,Catch.resid.log.byyear$cluhar_resid.log_mean,ylab="herring",xlab="month")
plot(Catch.resid.log.byyear$year,Catch.resid.log.byyear$urospp_resid.log_mean,ylab="red/white hake",xlab="month")
plot(Catch.resid.log.byyear$year,Catch.resid.log.byyear$merbil_resid.log_mean,ylab="silver hake",xlab="month")
plot(Catch.resid.log.byyear$year,Catch.resid.log.byyear$scosco_resid.log_mean,ylab="mackerel",xlab="month")
plot(Catch.resid.log.byyear$year,Catch.resid.log.byyear$pepspp_resid.log_mean,ylab="butterfish",xlab="month")
plot(Catch.resid.log.byyear$year,Catch.resid.log.byyear$ammspp_resid.log_mean,ylab="sandlance",xlab="month")
par(mfrow=c(1,1))

### STOP REVISING HERE
# compare average log values to average residuals (also on a log scale) by correlation
corr.byspecies.resid <- data.frame(matrix(nrow=6,ncol=1))
colnames(corr.byspecies.resid) <- "corr"
rownames(corr.byspecies.resid) <- c("herring","red/white hake","silver hake","mackerel","butterfish","sandlance")

corr.byspecies.resid[1,1] <- cor(Catch.log.byyear$log_cluhar_mean,Catch.resid.log.byyear$cluhar_resid.log_mean)
corr.byspecies.resid[2,1] <- cor(Catch.log.byyear$log_urospp_mean,Catch.resid.log.byyear$urospp_resid.log_mean)
corr.byspecies.resid[3,1] <- cor(Catch.log.byyear$log_merbil_mean,Catch.resid.log.byyear$merbil_resid.log_mean)
corr.byspecies.resid[4,1] <- cor(Catch.log.byyear$log_scosco_mean,Catch.resid.log.byyear$scosco_resid.log_mean)
corr.byspecies.resid[5,1] <- cor(Catch.log.byyear$log_pepspp_mean,Catch.resid.log.byyear$pepspp_resid.log_mean)
corr.byspecies.resid[6,1] <- cor(Catch.log.byyear$log_ammspp_mean,Catch.resid.log.byyear$ammspp_resid.log_mean)

corr.byspecies.resid

#######################################
# add spatial structure using a GAM, #5

# Fit GAM to model values, all years, using monthly residuals (#5a)
# All Hessians positive definite, depth significant except for scosco and ammspp, needed to drop depth as factor for pepspp to be positive definite
Catch.ich <- Catch.resid
fit.cluhar <- gam(cluhar_resid.log ~ s(lon, lat) + s(depth,k=3), data=Catch.ich); summary(fit.cluhar); gam.check(fit.cluhar)
fit.urospp <- gam(urospp_resid.log ~ s(lon, lat) + s(depth,k=3), data=Catch.ich); summary(fit.urospp); gam.check(fit.urospp)
fit.merbil <- gam(merbil_resid.log ~ s(lon, lat) + s(depth,k=3), data=Catch.ich); summary(fit.merbil); gam.check(fit.merbil)
fit.scosco <- gam(scosco_resid.log ~ s(lon, lat) , data=Catch.ich); summary(fit.scosco); gam.check(fit.scosco)
fit.pepspp <- gam(pepspp_resid.log ~ s(lon, lat) , data=Catch.ich); summary(fit.pepspp); gam.check(fit.pepspp)
fit.ammspp <- gam(ammspp_resid.log ~ s(lon, lat) , data=Catch.ich); summary(fit.ammspp); gam.check(fit.ammspp)
plot(fit.cluhar,pages=1,seWithMean=TRUE)
plot(fit.urospp,pages=1,seWithMean=TRUE)
plot(fit.merbil,pages=1,seWithMean=TRUE)
plot(fit.scosco,pages=1,seWithMean=TRUE)
plot(fit.pepspp,pages=1,seWithMean=TRUE)
plot(fit.ammspp,pages=1,seWithMean=TRUE)

# Add year term
# All Hessians positive definite
par(mfrow=c(1,1))
fit.cluhar <- gam(cluhar_resid.log ~ s(lon, lat) + s(depth,k=3) + s(year), data=Catch.ich); summary(fit.cluhar); gam.check(fit.cluhar)
fit.urospp <- gam(urospp_resid.log ~ s(lon, lat) + s(depth,k=3) + s(year), data=Catch.ich); summary(fit.urospp); gam.check(fit.urospp)
fit.merbil <- gam(merbil_resid.log ~ s(lon, lat) + s(depth,k=3) + s(year), data=Catch.ich); summary(fit.merbil); gam.check(fit.merbil)
fit.scosco <- gam(scosco_resid.log ~ s(lon, lat) + s(year), data=Catch.ich); summary(fit.scosco); gam.check(fit.scosco)
fit.pepspp <- gam(pepspp_resid.log ~ s(lon, lat) + s(year), data=Catch.ich); summary(fit.pepspp); gam.check(fit.pepspp)
fit.ammspp <- gam(ammspp_resid.log ~ s(lon, lat) + s(year), data=Catch.ich); summary(fit.ammspp); gam.check(fit.ammspp)

plot(fit.cluhar,pages=0,seWithMean=TRUE,select=3,ylim=c(-0.3,0.3),main="herring")
plot(fit.urospp,pages=0,seWithMean=TRUE,select=3,ylim=c(-0.2,0.2),main="red/white hake")
plot(fit.merbil,pages=0,seWithMean=TRUE,select=3,ylim=c(-0.1,0.1),main="silver hake")
plot(fit.scosco,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.1,0.1),main="mackerel")
plot(fit.pepspp,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.1,0.1),main="butterfish")
plot(fit.ammspp,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.4,0.4),main="sandlance")

# Add temperature term, not all terms are significant, makes minor differences in year trend
# but does not make sense to add temperature term given that sampling occurs throughout year
# analysis deleted

# plot gam documentation
# https://www.rdocumentation.org/packages/mgcv/versions/1.8-31/topics/plot.gam
# b<-gam(y~x0+s(x1)+s(x2)+s(x3))
# plot(b,pages=1,residuals=TRUE,all.terms=TRUE,shade=TRUE,shade.col=2)
# plot(b,pages=1,seWithMean=TRUE) ## better coverage intervals

par(mfrow = c(2,3))
vis.gam(fit.cluhar, plot.type="contour", too.far=0.15, main="Atlantic herring", xlab="Longitude",ylab="Latitude",color="topo")
#points(Catch.ich$lon, Catch.ich$lat,cex=0.2); 
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
vis.gam(fit.urospp, plot.type="contour", too.far=0.15, main="red white hake", xlab="Longitude",ylab="Latitude",color="topo")
#points(Catch.ich$lon, Catch.ich$lat,cex=0.2); 
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
vis.gam(fit.merbil, plot.type="contour", too.far=0.15, main="silver hake", xlab="Longitude",ylab="Latitude",color="topo")
#points(Catch.ich$lon, Catch.ich$lat,cex=0.2); 
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
vis.gam(fit.scosco, plot.type="contour", too.far=0.15, main="Atlantic mackerel", xlab="Longitude",ylab="Latitude",color="topo")
#points(Catch.ich$lon, Catch.ich$lat,cex=0.2); 
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
vis.gam(fit.pepspp, plot.type="contour", too.far=0.15, main="butterfish", xlab="Longitude",ylab="Latitude",color="topo")
#points(Catch.ich$lon, Catch.ich$lat,cex=0.2); 
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
vis.gam(fit.ammspp, plot.type="contour", too.far=0.15, main="sandlance", xlab="Longitude",ylab="Latitude",color="topo")
#points(Catch.ich$lon, Catch.ich$lat,cex=0.2); 
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")

########################################################################
# Fit GAM to log-transformed values, include year and month terms  (#5b)
# For these forms, Hessians were positive definite
Catch.ich <- Catch.resid
fit.cluhar <- gam(log_cluhar ~ s(lon, lat) + s(depth,k=3) + s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.cluhar); gam.check(fit.cluhar)
fit.urospp <- gam(log_urospp ~ s(lon, lat) + s(depth,k=3) + s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.urospp); gam.check(fit.urospp)
fit.merbil <- gam(log_merbil ~ s(lon, lat) + s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.merbil); gam.check(fit.merbil)
fit.scosco <- gam(log_scosco ~ s(lon, lat) + s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.scosco); gam.check(fit.scosco)
fit.pepspp <- gam(log_pepspp ~ s(lon, lat) + s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.pepspp); gam.check(fit.pepspp)
fit.ammspp <- gam(log_ammspp ~ s(lon, lat) + s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.ammspp); gam.check(fit.ammspp)

par(mfrow = c(2,3))
#by year
plot(fit.cluhar,pages=0,seWithMean=TRUE,select=3,ylim=c(-0.3,0.3),main="herring")
plot(fit.urospp,pages=0,seWithMean=TRUE,select=3,ylim=c(-0.2,0.2),main="red/white hake")
plot(fit.merbil,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.1,0.1),main="silver hake")
plot(fit.scosco,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.1,0.1),main="mackerel")
plot(fit.pepspp,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.1,0.1),main="butterfish")
plot(fit.ammspp,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.4,0.4),main="sandlance")

#by month
plot(fit.cluhar,pages=0,seWithMean=TRUE,select=4,ylim=c(-0.3,0.3),main="herring")
plot(fit.urospp,pages=0,seWithMean=TRUE,select=4,ylim=c(-0.4,0.4),main="red/white hake")
plot(fit.merbil,pages=0,seWithMean=TRUE,select=3,ylim=c(-0.3,0.3),main="silver hake")
plot(fit.scosco,pages=0,seWithMean=TRUE,select=3,ylim=c(-0.3,0.3),main="mackerel")
plot(fit.pepspp,pages=0,seWithMean=TRUE,select=3,ylim=c(-0.3,0.3),main="butterfish")
plot(fit.ammspp,pages=0,seWithMean=TRUE,select=3,ylim=c(-0.6,0.6),main="sandlance")

# spatial plots
vis.gam(fit.cluhar, plot.type="contour", too.far=0.15, main="Atlantic herring", xlab="Longitude",ylab="Latitude",color="topo")
#points(Catch.ich$lon, Catch.ich$lat,cex=0.2); 
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
vis.gam(fit.urospp, plot.type="contour", too.far=0.15, main="red white hake", xlab="Longitude",ylab="Latitude",color="topo")
#points(Catch.ich$lon, Catch.ich$lat,cex=0.2); 
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
vis.gam(fit.merbil, plot.type="contour", too.far=0.15, main="silver hake", xlab="Longitude",ylab="Latitude",color="topo")
#points(Catch.ich$lon, Catch.ich$lat,cex=0.2); 
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
vis.gam(fit.scosco, plot.type="contour", too.far=0.15, main="Atlantic mackerel", xlab="Longitude",ylab="Latitude",color="topo")
#points(Catch.ich$lon, Catch.ich$lat,cex=0.2); 
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
vis.gam(fit.pepspp, plot.type="contour", too.far=0.15, main="butterfish", xlab="Longitude",ylab="Latitude",color="topo")
#points(Catch.ich$lon, Catch.ich$lat,cex=0.2); 
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
vis.gam(fit.ammspp, plot.type="contour", too.far=0.15, main="sandlance", xlab="Longitude",ylab="Latitude",color="topo")
#points(Catch.ich$lon, Catch.ich$lat,cex=0.2); 
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")

#####################################################################################################################
# Log link function (#5c)
# https://www.theanalysisfactor.com/the-difference-between-link-functions-and-data-transformations/
# https://www.r-bloggers.com/generalized-linear-models-understanding-the-link-function/
# https://stat.ethz.ch/R-manual/R-devel/library/stats/html/family.html

# Add link function, gaussian distribution, log link
# Does not converge after 30 minutes, even when dropping by=year
# fit.cluhar <- gam(cluhar_100m3 ~ s(lon, lat, by=year) + s(depth) + s(year) + s(month,bs='cc'), family=gaussian(link="log"),data=Catch.ich)
# fit.cluhar <- gam(cluhar_100m3 ~ s(lon, lat) + s(depth) + s(year) + s(month,bs='cc'), family=gaussian(link="log"),data=Catch.ich)

# log link, quasi distribution, same error structure as Wood 2005, for sole eggs
# converges after 7 minutes
# results nonsensical, likely not converged (see google presentation)
 fit.cluhar <- gam(cluhar_100m3 ~ s(lon, lat, by=year) + s(depth) + s(year) + s(month,bs='cc'), family=quasi(link="log",variance = "mu"),data=Catch.ich)
 warnings() # step size truncated due to divergence

# drop by=year
fit.cluhar <- gam(cluhar_100m3 ~ s(lon, lat) + s(depth) + s(year) + s(month,bs='cc'), family=quasi(link="log",variance = "mu"),data=Catch.ich)
warnings()
# Iteration limit reached without full convergence - check carefully

#####################################################################################################################
# Add spatial plots BY YEAR
# Use log-transformed data (rather than residuals), which seems more straightforward
# Case 1: full (year, month, location, depth for two species) 
#####################################################################################################################
# Table of names
name.tbl <- data.frame(spp.num = seq(25,30),spp.name = c("herring","red white hake","silver hake","mackerel","butterfish","sandlance"))

# Years to plot (1988-2016, sampling effort from 1988-1998 is less balanced seasonally)
#  (MARMAP, 1977 - 1987) and Ecosystem Monitoring (EcoMon, 1999 - present)
table(Catch.rpt$year,Catch.rpt$month)

####################
# Case 1: Full model
# by= In the numeric vector case the elements multiply the smooth, evaluated at the corresponding covariate values 
# (a `varying coefficient model' results).
# Excluded depth as a factor, but could include for first two species (k=3) with if ... else statement
start.year <- 1977   # all years, can change start year if desired (e.g., Ecomon years: 1989-present)
Catch.ich <- subset(Catch.rpt, year>=start.year)

pdf("gam spatial year month model.pdf")

opar <- par(mfrow = c(2,2), cex=1.0, omi=c(0.4,0.4,0,0),mar=c(3,3,3,1), pch=16)
for (i in 25:30) {
    j <- i - 24   # counter for species names
    fit.spp <- gam(Catch.ich[,i] ~ s(lon, lat, by=year) + s(year) + s(month,bs='cc'), data=Catch.ich, dist ="normal")
  plot(fit.spp,pages=0,seWithMean=TRUE,select=2,ylim=c(-0.4,0.4),main=paste(name.tbl[j,2],"by year"))
  plot(fit.spp,pages=0,seWithMean=TRUE,select=3,ylim=c(-0.4,0.4),main=paste(name.tbl[j,2],"by month"))

  # by year plots  
  for (k in start.year:2016) {
    vis.gam(fit.spp, plot.type="contour", too.far=0.15, cond=list(year=k),main=paste(name.tbl[j,2],k), xlab="Longitude",ylab="Latitude",color="topo") 
    points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
  }
}
par(opar)
dev.off()

####################################################################
# Case 3: Treat year as a factor instead of as a continuous variable
# Reduced model (eliminate depth as an independent variable)
# Developmental model

# From https://www.rdocumentation.org/packages/mgcv/versions/1.8-31/topics/s
# 'by' is a numeric or factor variable of the same dimension as each covariate. In the numeric vector case the elements multiply the smooth, 
# evaluated at the corresponding covariate values (a `varying coefficient model' results). For the numeric by variable case the resulting 
# smooth is not usually subject to a centering constraint (so the by variable should not be added as an additional main effect). 
# In the factor by variable case a replicate of the smooth is produced for each factor level (these smooths will be centered, 
# so the factor usually needs to be added as a main effect as well). See gam.models for further details. A by variable may also be 
# a matrix if covariates are matrices: in this case implements linear functional of a smooth (see gam.models and linear.functional.terms for details).
# id = 1 forces all s(as.factor(year)) terms to have the same smoothing param ### runs a little faster, plots qualitatively similar
# Also see: https://www.rdocumentation.org/packages/mgcv/versions/1.8-31/topics/gam.models

# Select data set and configure year grouping
# data set starts in 1977; EcoMon survey began in 1989
start.year <- 1977
Catch.ich <- subset(Catch.rpt, year>=start.year)                                # Note use of Catch.rpt
                                                                                # only contains raw and log-transformed values
year.int <- 5                                                                   # Size of year group (1, 3, 5, 10)
Catch.ich$year.group <- ceiling((Catch.ich$year-1971+1)/year.int)               # create year grouping (e.g., 1977 = year.group 1 (1971-1975))
table(Catch.ich$year.group,Catch.ich$month)                                     # tabulate year group by month
min.year.group <- min(Catch.ich$year.group)
max.year.group <- max(Catch.ich$year.group)
num.year.group <- max.year.group - min.year.group + 1                           # determine the number of year groups

# Run model
# All Hessians were positive definite
pdf("gam spatial year group month model.pdf")
opar <- par(mfrow = c(2,2), cex=1.0, omi=c(0.4,0.4,0,0),mar=c(3,3,3,1), pch=16)
for (i in 25:30) {       
    j <- i - 24   # counter for species names

    fit.spp <- gam(Catch.ich[,i] ~ s(lon, lat, by=as.factor(year.group),id=0) + s(year) + s(month,bs='cc'), data=Catch.ich, dist ="normal")
  summary(fit.spp); gam.check(fit.spp)

  plot(fit.spp,pages=0,seWithMean=TRUE,select=num.year.group+1,ylim=c(-0.4,0.4),main=paste(name.tbl[j,2],"by year"))
  plot(fit.spp,pages=0,seWithMean=TRUE,select=num.year.group+2,ylim=c(-0.4,0.4),main=paste(name.tbl[j,2],"by month"))

  # by year group plots  
  for (k in min.year.group:max.year.group) {
          vis.gam(fit.spp, view=c("lon","lat"), plot.type="contour", too.far=0.15, cond=list(year.group=k),
            main=paste(name.tbl[j,2]," year group ",k*year.int+1971-5), xlab="Longitude",ylab="Latitude",color="topo") 
    points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
  }
}
par(opar)
dev.off()
####################################################################
# FYI, last year's maps were based on non-zeros (i.e., presence only model: looks quite different)
# this year's maps look more sensible: patterns are more related to bathymetry and, to some extent, temperature 
# (e.g, sand lance to the north)

##############################################################################################
##############################################################################################
##############################################################################################
# Program for statistics lab
# Gulf of Maine plankton data (EcoMon surveys)
# In some years, sampling was monthly, but in recent years, sampling focused on April-June and August-November
# Michael Sigler
# May 19, 2019, Revised November 22, 2019, March 24, 2020

options(scipen=999)  # turn-off scientific notation like 1e+48
library(fields)
library(dplyr)
library(mgcv)
library(mapdata)
library(mapproj)

############################################################################################################
# Read plankton file
# This file is compiled from the original database available at:
# XXX
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
# urospp_100m3 = red and white count per 100 m3 of volume
# e.g., log_cluhar = log 10 transformed (value +1) of cluhar_100m3
# distance = distance (km) from sample location to Appledore Island 
# Appledore Island location
Lat.SML <- 42.987727
Lon.SML <- -70.613940

Catch <- as.data.frame(read.csv("Plankton.csv",header=T))
Catch$date <- as.Date(Catch$date, "%Y-%m-%d")		# Assign Date format 'yyyy-mm-dd'
xtabs(~ year + month, data = Catch)
head(Catch)   # note that the output file "Plankton.csv" has a record number added to the first column

# Years to plot (1988-2016, sampling effort from 1988-1998 is less balanced seasonally)
#  (MARMAP, 1977 - 1987) and Ecosystem Monitoring (EcoMon, 1999 - present)
start.year <- 1977  # all years, can change start year if desired (e.g., Ecomon years: 1989-present)
Catch.ich <- subset(Catch, year>=start.year)         

########################################################################
# Fit GAM to log-transformed values, include year and month terms  (#5b)
# One spatial model for all years
par(mfrow = c(2,2))
fit.cluhar <- gam(log_cluhar ~ s(lon, lat) + s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.cluhar); gam.check(fit.cluhar)
fit.urospp <- gam(log_urospp ~ s(lon, lat) + s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.urospp); gam.check(fit.urospp)
fit.merbil <- gam(log_merbil ~ s(lon, lat) + s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.merbil); gam.check(fit.merbil)
fit.scosco <- gam(log_scosco ~ s(lon, lat) + s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.scosco); gam.check(fit.scosco)
fit.pepspp <- gam(log_pepspp ~ s(lon, lat) + s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.pepspp); gam.check(fit.pepspp)
fit.ammspp <- gam(log_ammspp ~ s(lon, lat) + s(year) + s(month,bs='cc'), data=Catch.ich); summary(fit.ammspp); gam.check(fit.ammspp)
par(mfrow = c(1,1))

pdf("gam spatial model.pdf")
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
#points(Catch.ich$lon, Catch.ich$lat,cex=0.2); 
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
vis.gam(fit.urospp, plot.type="contour", too.far=0.15, main="red white hake", xlab="Longitude",ylab="Latitude",color="topo")
#points(Catch.ich$lon, Catch.ich$lat,cex=0.2); 
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
vis.gam(fit.merbil, plot.type="contour", too.far=0.15, main="silver hake", xlab="Longitude",ylab="Latitude",color="topo")
#points(Catch.ich$lon, Catch.ich$lat,cex=0.2); 
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
vis.gam(fit.scosco, plot.type="contour", too.far=0.15, main="Atlantic mackerel", xlab="Longitude",ylab="Latitude",color="topo")
#points(Catch.ich$lon, Catch.ich$lat,cex=0.2); 
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
vis.gam(fit.pepspp, plot.type="contour", too.far=0.15, main="butterfish", xlab="Longitude",ylab="Latitude",color="topo")
#points(Catch.ich$lon, Catch.ich$lat,cex=0.2); 
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
vis.gam(fit.ammspp, plot.type="contour", too.far=0.15, main="sandlance", xlab="Longitude",ylab="Latitude",color="topo")
#points(Catch.ich$lon, Catch.ich$lat,cex=0.2); 
points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")

par(opar)
dev.off()

########################################################################
# Fit GAM to log-transformed values, include year and month terms  (#5b)
# Spatial model by year groups (5-year groups)
# From https://www.rdocumentation.org/packages/mgcv/versions/1.8-31/topics/s
# In the GAM, 'by' is a numeric or factor variable of the same dimension as each covariate. In the numeric vector case the elements multiply the smooth, 
# evaluated at the corresponding covariate values (a `varying coefficient model' results). For the numeric by variable case the resulting 
# smooth is not usually subject to a centering constraint (so the by variable should not be added as an additional main effect). 
# In the factor by variable case a replicate of the smooth is produced for each factor level (these smooths will be centered, 
# so the factor usually needs to be added as a main effect as well). See gam.models for further details. A by variable may also be 
# a matrix if covariates are matrices: in this case implements linear functional of a smooth (see gam.models and linear.functional.terms for details).
# id = 1 forces all s(as.factor(year)) terms to have the same smoothing param, but gives nonsensical results for herring
# Also see: https://www.rdocumentation.org/packages/mgcv/versions/1.8-31/topics/gam.models

# Table of names
name.tbl <- data.frame(spp.num = seq(26,31),spp.name = c("herring","red white hake","silver hake","mackerel","butterfish","sandlance"))

# Year grouping
Catch.ich <- Catch
year.int <- 3                                                                   # Size of year group (1, 3, 5, 10)
Catch.ich$year.group <- ceiling((Catch.ich$year-1971+1)/year.int)               # create year grouping (e.g., 1973 = year.group 1 (1971-1975))
table(Catch.ich$year.group,Catch.ich$month)                                     # tabulate year group by month
min.year.group <- min(Catch.ich$year.group)
max.year.group <- max(Catch.ich$year.group)
num.year.group <- max.year.group - min.year.group + 1                           # determine the number of year groups

# Run model
pdf("gam spatial by year group.pdf")
par(mfrow = c(2,2))
for (i in 26:31) {       
  j <- i - 25   # counter for species names

  fit.spp <- gam(Catch.ich[,i] ~ s(lon, lat, by=as.factor(year.group),id=0) + s(year) + s(month,bs='cc'), data=Catch.ich, dist ="normal")
  summary(fit.spp); # gam.check(fit.spp)
  
  plot(fit.spp,pages=0,seWithMean=TRUE,select=num.year.group+1,ylim=c(-0.4,0.4),main=paste(name.tbl[j,2],"by year"))
  plot(fit.spp,pages=0,seWithMean=TRUE,select=num.year.group+2,ylim=c(-0.4,0.4),main=paste(name.tbl[j,2],"by month"))
  
  # by year group plots  
  for (k in min.year.group:max.year.group) {
    vis.gam(fit.spp, view=c("lon","lat"), plot.type="contour", too.far=0.15, cond=list(year.group=k),
            main=paste(name.tbl[j,2]," year group ",k*year.int+1971-5), xlab="Longitude",ylab="Latitude",color="topo") 
    points(Lon.SML, Lat.SML, pch=16, col = "red"); map('worldHires',fill=T,add=T, col="grey")
  }
}
par(mfrow = c(1,1))
dev.off()
####################################################################






