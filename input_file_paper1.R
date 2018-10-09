########################################loading the required packages##########################################

#### This code is specific to SMAPVEX12. If you want to use it for your own dataset, go straight to parameter esimation file


rm(list=ls())

library(fields)
library(raster)
library(GenSA)
##Package ncdf4 should be installed on the system
setwd("/Users/dhruvakathuria/OneDrive - Texas A&M University/TAMU_paper_1/data") # setup working directory

################################################################################################################
#####################################Setting up airborne data###################################################
################################################################################################################

filenames <- list.files(pattern="SV12", full.names = TRUE) # the files are available on NSIDC under SMAPVEX12 validation data
SM_files <- lapply(filenames, read.csv)  # This is a list where each element of the list has individual days of 
# airborne data

SM_extract <- function(x)
{
  ind_nan <- which(is.nan(x$VSM))
  SM <- x[-ind_nan,c(4,5,6)]
  colnames(SM) <- c("X","Y","SM")
  SM
}

SM_air <- lapply(SM_files,SM_extract) ## This list subsets the location and SM for individual days.Your subset
# should be such that each of the list is a dataframe having three columns in the order (Easting, Northing, Soil Moisture)

Dates <- read.csv("airborne data/Dates.csv") ## We name the individual lists by the Dates ## we stored the dates in a separate . csv file
names(SM_air) <- Dates$Var1

################################################################################################################
#####################################Setting up finescale data###################################################
################################################################################################################

mydata_a <- read.table("smapvex12_point scale values_sm_by_day.csv", header=TRUE, sep=",")  ## file available on NSIDC under SMAPVEX12 validation data, this is the probe soil moisture data (not the core)
mydata <- mydata_a[,-c(2,3)]  ## Fine scale soil moisture for all the days is in a single file. We subset it such
## that "mydata" has three columns in the order (Date, Soil moisture, Site-ID)

Dates <- as.character(unique(mydata$Sample_Dat, incomparables = FALSE)) ## We now subset SM values by dates 
fine_SM <- list()
for (i in 1:length(Dates))
{
  fine_SM[[i]] <- subset(mydata,Sample_Dat==Dates[i],select=c(Soil_Moist, Site_ID))
}
names(fine_SM) <- Dates  ## fine_SM is a list. The name of each element of the list is the date. Each element of fine_SM is a 
## data frame with two columns in the order (Soil moisture, Site_ID)


values <- function(SM1){ ## Since we have three replicate readings at each point, 
  Sites <- as.vector(unique(SM1$Site_ID, incomparables = FALSE)) ##we define a function "values" which averages 
  SM1_average <- data.frame()                                    ## the 3 observations at every point
  for (i in 1:length(Sites))
  {vec=0
  vec <- which(SM1$Site_ID==Sites[i])
  SM1_average[i,1] <- Sites[i]
  SM1_average[i,2] <- mean(SM1$Soil_Moist[vec])
  colnames(SM1_average) <- c("Site_ID", "SM")
  }
  return(SM1_average)
} 


resa <- lapply(fine_SM,values) ## Applying function "values" to all the files
values1 <- function(dfile){(dfile$SM)} ## defiing a function to extract only the Soil moisture values
res <- lapply(resa,values1) ## Applying function "values1" to all the files, "res" is a list, where each element 
## of res is a vector of SM observations for a date

##Putting site_id back in
SM <- list()
for (i in 1:length(resa))
{
  SM[[i]] <- cbind(resa[[i]][1],res[[i]])
  names(SM[[i]]) <- c("Site_ID", "SM")
}

names(SM) <- names(resa) ## resa is a list where each element of the list is named by the date and consists
## of a dataframe with two columns (Site_ID, SM). Now SM is the averaged of the three
## replicate readings for each point

###merging location with SM using Site-id ###### the location data is available at NSIDC
locations <- read.csv("location_values.csv")[,c(3,4,8)] ## this is the locations file. Subset such that the dataframe with three
## columns in the order (Easting, Northing, Site_ID)
D <- function(a)
{
  z <- data.frame()
  z <- merge(a,locations,by="Site_ID",sort=FALSE)
}

merged_data <- lapply(SM,D) ## This is a list where every element of the list is a dataframe with 4 columns
## in the order (Site_ID, SM, Easting, Northing)

SM_fine <- merged_data


SM_fine1 <- list(SM_fine$`6/12/2012`, SM_fine$`6/15/2012`, SM_fine$`6/17/2012`,SM_fine$`6/22/2012`,
                 SM_fine$`6/25/2012`, SM_fine$`6/27/2012`, SM_fine$`6/29/2012`, SM_fine$`7/3/2012`,
                 SM_fine$`7/5/2012`, SM_fine$`7/8/2012`, SM_fine$`7/10/2012`, SM_fine$`7/13/2012`,
                 SM_fine$`7/14/2012`, SM_fine$`7/17/2012`, SM_fine$`7/19/2012`)
names(SM_fine1) <- names(SM_air)  ## We subset SM_fine for the days as SM_fine1 for which the airborne data 
##is available and give it the elemets the name of the dates



################################################################################################################
#####################################LAI_data###################################################################
################################################################################################################

LAI1 <- brick("MCD15A3H.006_aid0001.nc", varname="Lai_500m") ## LAI file from MODIS (MCD15A3H, version 6)


m <- (LAI1$X2012.06.13 - LAI1$X2012.06.09)/4   ## For the intervening days, we do linear interpolation of LAI as
LAI_Jun12 <- m*3 + LAI1$X2012.06.09            ## mentioned in the manuscript

m <- (LAI1$X2012.06.17 - LAI1$X2012.06.13)/4
LAI_Jun15 <- m*2 + LAI1$X2012.06.13

LAI_Jun17 <- LAI1$X2012.06.17

m <- (LAI1$X2012.06.25 - LAI1$X2012.06.21)/4
LAI_Jun22 <- m*1 + LAI1$X2012.06.21

LAI_Jun25 <- LAI1$X2012.06.25

m <- (LAI1$X2012.06.29 - LAI1$X2012.06.25)/4
LAI_Jun27 <- m*2 + LAI1$X2012.06.25

LAI_Jun29 <- LAI1$X2012.06.29

LAI_Jul3 <- LAI1$X2012.07.03

m <- (LAI1$X2012.07.07 - LAI1$X2012.07.03)/4
LAI_Jul05 <- m*2 + LAI1$X2012.07.03

m <- (LAI1$X2012.07.11 - LAI1$X2012.07.07)/4
LAI_Jul08 <- m*1 + LAI1$X2012.07.07

m <- (LAI1$X2012.07.11 - LAI1$X2012.07.07)/4
LAI_Jul10 <- m*3 + LAI1$X2012.07.07

m <- (LAI1$X2012.07.15 - LAI1$X2012.07.11)/4
LAI_Jul13 <- m*2 + LAI1$X2012.07.11

LAI_Jul14 <- m*3 + LAI1$X2012.07.11

m <- (LAI1$X2012.07.19 - LAI1$X2012.07.15)/4
LAI_Jul17 <- m*2 + LAI1$X2012.07.15

LAI_Jul19 <- LAI1$X2012.07.19



LAI2 <- list(June_12=LAI_Jun12 ,June_15=LAI_Jun15, June_17=LAI_Jun17, June_22=LAI_Jun22,
             June_25=LAI_Jun25, June_27=LAI_Jun27, June_29=LAI_Jun29, July_3=LAI_Jul3,
             July_5=LAI_Jul05,July_8=LAI_Jul08,July_10=LAI_Jul10,
             July_13=LAI_Jul13, July_14=LAI_Jul14, July_17=LAI_Jul17,
             July_19=LAI_Jul19)  ## We append the LAI values into a list. Note that individual elements of LAI2
## are raster files, since we used the "brick" command above

data1 <- list()
for (i in 1:(length(SM_air)))
{
  data1[[i]] <- list(Fine=SM_fine1[[i]], Air=SM_air[[i]], name=names(SM_air[i]),LAI=LAI2[[i]])
}  ## data1 is the input list. It has 4 elements corresponding to each day 
##(Data frame of fine_Scale_SM, Data frame of Airborne SM, Date, LAI raster)



library(gdata) 
keep(data1, sure=TRUE)
detach("package:gdata", unload=TRUE)  ## We remove all the variables from memory except data1

jj = 2  ## the index jj is for individual days (use any jj from 2, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14) as we used
## these 11 days in the study
data2 <- data1[[jj]] ## data2 is the list for the particular chosen day
rm(data1)

#####

sr <- "+proj=utm +zone=14N +ellps=GRS80 +datum=NAD83" ## This is the UTM for Northing/Easting for the area

##Clay and silt rasters 
## For any day in the airborne files (since its static), extract the % clay and % sand and use rasterFromXYZ(x, y, Clay)/rasterFromXYZ(x, y, Sand) to get your
## rasters. for silt, just subtract (sand + clay from 100)
raster_clay_airborne <- raster("Clay_Smapvex12.grd")  ## We converted the soil texture data into raster files, 
raster_silt_airborne <- raster("Silt_Smapvex12.grd")  ## a shapefile would work as well

##Algorithm function

input <- data2
L <- 1500  ## support of airborne pixel in meters
name <- input$name  ## Date

f1d <- input$Fine ##fine scale data
obs.grd1 <- f1d[,3:4] ## x, y coordinate of each SM observation

obs.grd2 <- extract(raster_clay_airborne,obs.grd1)  ## There are some in situ observations, which are outside the
ind_obs <- which(is.na(obs.grd2))                   ## extent of the airborne data for some of the days, we are 
## removing those in situ points
if(length(ind_obs)==0)
{
  obs.grd_del <- obs.grd1
  z <- f1d$SM  ## z is the SM observation
} else
  
{
  
  obs.grd_del <- obs.grd1[-ind_obs,]
  z <- f1d$SM[-ind_obs] 
}

obs.grd  <- obs.grd_del  ## x,y coordinates after removing the above points

raster_LAI1 <- projectRaster(input$LAI,crs=sr) ## We project the MODIS LAI raster to UTM
raster_LAI1[raster_LAI1[]>10]=NA ## We set all values for LAI > 10 to be NA as is usually done

LAI <- extract(raster_LAI1,obs.grd)     ## extracting values for covariates for our locations
Clay <- extract(raster_clay_airborne,obs.grd)
Silt <- extract(raster_silt_airborne,obs.grd)

Covariates <- cbind(LAI,Clay,Silt)
ind_cov <- which(is.na(Covariates),TRUE)[,1]  ## Here we remove all the points for which we dont have any of the
ind_cov1 <- unique(ind_cov)                   ## covariate values. For most of the days there are no such points
if(length(ind_cov1)==0)
{
  Covariates=Covariates
}else
{
  
  Covariates <- Covariates[-ind_cov1,]
  z <- z[-ind_cov1]
  obs.grd <- obs.grd[-ind_cov1,]
  
}

Covariates1 <- Covariates
Covariates <- scale(Covariates)  ## We scale our Covariates prior to our analysis

##Defining X for mean and covariance
X_wt <- cbind(rep(1,dim(Covariates)[1]),Covariates) ## for the covariance

X_1 <- X_wt ## X_1 is the for the mean. We take it as the same as the covariance in our analysis, though this can
## be changed
