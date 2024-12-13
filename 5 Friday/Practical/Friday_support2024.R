#### Uncertainty propagation under changing spatial support
##   (Friday afternoon practical)
## 
##   Part of PE-RC course "Uncertainty Analysis and Statistical Validation 
##   of Spatial Environmental Models (UA&SV)",  9 - 13 December 2024, Wageningen.
##
##   Created by Luc Steinbuch, Gerard Heuvelink and Sytze de Bruin, Wageningen UR.
##   Largely based on earlier work by Bertin Takoutsing (CIFOR-ICRAF, Wageningen UR; 
##   see for context among others https://doi.org/10.1007/s11119-024-10200-6 ) 
##   and Sriram Jallu (Wageningen UR). 
##  
##   Please note that both the script and the provided data are meant for this 
##   educational setting. We strongly recommend you to double-check the code
##   before applying it to your own data set.
##

#### Preparation: package, files, settings ####

##  Make sure you have in the same directory as this script:
##  "df_realisations_soil_prop.Rda", "crop_parameters.csv",
##  "QUEFTS_functions.R" and the directory "divisions" containing
##  several shape files. 

## Load package
library(terra)

## Set working directory to the same directory as this script.
setwd("C:/PERC/Friday")

## Define the Coordinate Reference System we will be using for all spatial objects
target_crs <- "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +units=m +ellps=WGS84 +datum=WGS84"

## Load the shapefiles with the administrative "divisions" (in French: 
## "départements") in the West-Region of Cameroon . These divisions 
## will later be used for spatial aggregation. They are represented 
# as spatial polygon objects:
sppol_divisions <- vect("divisions/West_Region_Cameroon.shp")
sppol_divisions <- project(sppol_divisions, target_crs)
division_names <- sppol_divisions$NAME_2

## Let's have a quick look
plot(sppol_divisions, 
     col="lightblue", 
     border="darkblue",
     main = "Divisions in region West-Cameroon")
text(sppol_divisions, labels=division_names, cex=0.8, col="red")

## load the QUEFTS (crop yield simulation) functions
source("QUEFTS_functions.R")

## Load df_realizations_soil_prop containing 100 realizations of 4 soil props
## simulated using a geostatistical model for the residuals and a trend 
## predicted by a random forest model. The script can be found in 
## UA_SV_day5_Preparation_soil_property_realisations.zip

load(file= "df_realisations_soil_prop.Rda")

#### Explore soil property realizations ####

## Explore the data.frame we just loaded
dim(df_realisations_soil_prop)
df_realisations_soil_prop[1:5, 1:8]
df_realisations_soil_prop[1:5, 105:110]

npix <- nrow(df_realisations_soil_prop) ## number of pixels on raster map
nreal <- (ncol(df_realisations_soil_prop)-2)/4 ## number of realizations


## Create a terra raster object using the coordinates and the 5th column 
ra_example_map <- rast(df_realisations_soil_prop[,c(1:2, 5)], 
                       crs = target_crs)
names(df_realisations_soil_prop)[5]
plot(ra_example_map, main = "Realization #3 of soil property pH")
lines(sppol_divisions)

ra_example_map <- rast(df_realisations_soil_prop[,c(1:2, 6)], 
                       crs = target_crs)
names(df_realisations_soil_prop)[6]
plot(ra_example_map, main = "Realization #4 of soil property pH")
lines(sppol_divisions)

ra_example_map <- rast(df_realisations_soil_prop[,c(1:2, 250)], 
                       crs = target_crs)

names(df_realisations_soil_prop)[250]
plot(ra_example_map, main = "Realization #48 of soil property Olsen P")
lines(sppol_divisions, col="red")

summary(ra_example_map)

hist(ra_example_map)


## Feel free to explore any other soil property realization you are
## interested in!


#### Simulate maize yield realizations ####

## For more information about the used QUEFTS model, see for example
## https://doi.org/10.1007/s11119-024-10200-6, section 3.2. In short, 
## it predicts crop yield based on available nutrients and an input
## water limited yield prediction using empirical rules.

## First we have to set several global QUEFTS parameters

## Default fertilization  
## Feel free to change, but note that zeros will give zero yield
N_fert <- 120
P_fert <- 70
K_fert <- 70

## Prepare QUEFTS parameters
settings <- data.frame (waterlimited_yield = 5000,
                        recover_N = 0.5,
                        recover_P = 0.1,
                        recover_K = 0.5,
                        AvgTemp = 25)

crop_pars <- read.csv("crop_parameters.csv")
crop_pars <- subset(crop_pars,crop == "Maize")
crop_id <- crop_pars["crop_id"]


## Get default soil parameters
## We need some dummy variables, because this QUEFTS implementation
## expects a particular data structure
df_dummy <- data.frame(pH = NA,
                       SOC = NA, 
                       olsenP = NA,
                       exchK = NA)
tmp <- quefts_soil_parameters(df_dummy)
soil_prop <- tmp[[1]]
soil_param <- tmp[[2]]
soil_param$temp<-settings$AvgTemp 

fert_n <- settings$recover_N * N_fert  # N uptake
fert_p <- settings$recover_P * P_fert  # P uptake
fert_k <- settings$recover_K * K_fert  # N uptake

crop_param <- quefts_crop_parameters(crop_id = crop_id, 
                                     crop_pars = crop_pars)[[1]]

## >> Finished setting the global QUEFTS parameters <<


## In the file QUEFTS_functions.R, we prepared the function 
## `calculate_yield` which takes the global variables 
## and the four soil properties, and returns the yield in kg/ha.
## A few examples:

                  ### pH, SOC, P,  K
calculate_yield(x = c(7,   50, 15, 0.04))
calculate_yield(x = c(3,   10, 4, 0.01))

####
## We can apply this function also on raster level. As example
## we show here the 1st realization:

vn_index_columns <- c(1,2,3,103,203,303)
## check if we have the right columns
names(df_realisations_soil_prop)[vn_index_columns]

## Create one "spatial raster stack" of those four layers
rs_one_soil_prop_realisation <- rast(df_realisations_soil_prop[,vn_index_columns], 
                                     crs = target_crs)

## check visually
plot(rs_one_soil_prop_realisation)

## Calculate the yield as spatial raster object
ra_yield_one_realisation <- app(rs_one_soil_prop_realisation, 
                               fun = calculate_yield)

plot(ra_yield_one_realisation, main = "Yield for first realisation")


###
## Using above principle, we create a loop and fill a list with all 
## yield realizations. This will take several minutes of calculation 
## time.

## Create empty list
li_yield_realisations <- list()

## For our own convenience: Progress indication
pb = txtProgressBar(min = 0, 
                    max = nreal, 
                    initial = 0) 

## Loop over all soil property realizations
for (i in 1:nreal)
{
  vn_index_columns <- c(1,2,        # Coordinates
                        i+2+0*nreal, # pH
                        i+2+1*nreal, # SOC
                        i+2+2*nreal, # P
                        i+2+3*nreal) # K
  
  ## Create one "spatial raster stack"
  rs_one_soil_prop_realisation <- rast(df_realisations_soil_prop[,vn_index_columns], 
                                       crs = target_crs)
  
  ## Calculate the yield as spatial raster object
  ra_yield_one_realisation <- app(rs_one_soil_prop_realisation, 
                                 fun = calculate_yield)
  
  ## Add to list
  li_yield_realisations[[i]] <- ra_yield_one_realisation
  
  setTxtProgressBar(pb,i)

}  ## End loop over all soil property realisations

close(pb)


## Create one raster stack from list with single raster layers
rs_yield_realisations <- rast(li_yield_realisations)

## Assign an unique name to each yield realization
names(rs_yield_realisations) <- paste0("Yield_", 1:nreal)


#### Explore yield realizations per layer ####

summary(rs_yield_realisations)

plot(rs_yield_realisations, 
     y = 1, # select layer
     main = "Predicted maize yield, yield realisation #1")

plot(rs_yield_realisations$Yield_2, ## Another possibility to select a layer
     main = "Predicted maize yield, yield realisation #2")



#### Explore yield realizations per pixel ####

## Yield per pixel; average over all yield realizations
ra_avg_yield <- app(rs_yield_realisations, mean, na.rm=TRUE)
plot(ra_avg_yield, main = "Yield per pixel averaged over 100 yield realizations")

## Yield standard deviation per pixel; sd over all yield realizations
ra_sd_yield <- app(rs_yield_realisations, sd, na.rm=TRUE)
plot(ra_sd_yield, main = "SD yield per pixel over 100 yield realizations")

####
## Yield prediction interval width per pixel
## In this exercise, we define the prediction interval as 
## the difference between the 95th and 5th percentile, i.e., the 90% PI

## For understanding the code: we demonstrate the helper function in
## two steps

## Step 1: quantile
quantile(x = 100:300, # A vector with values
         probs = c(0.05, 0.95), # The under- and upper quantile
         names = FALSE) # Force to return plain numbers

## Step 2: the difference between the two numbers
diff(
  quantile(x = 100:300, 
              probs = c(0.05, 0.95), 
              names = FALSE)
)

## For the map of the prediction interval over all realizations, 
## we apply above step 2:
ra_yield_pred_interv <- app(rs_yield_realisations, 
                            fun = function(x) {diff( 
                              quantile(x, 
                                       probs = c(0.05, 0.95), 
                                       names = FALSE,
                                       na.rm = TRUE) 
                              )}
                            )
plot(ra_yield_pred_interv, 
     main = "Yield prediction interval width (90%) over 100 yield realizations")


###
## We can also answer more involved questions. Such as: "Which percentage 
##  of the realizations has a yield below a given threshold?
threshold <- 3200 
ra_below_threshold <- app(rs_yield_realisations, 
                       fun = function(x) {100*mean(x < threshold, na.rm=TRUE)})

plot(ra_below_threshold, 
     main = paste0("Percentage of realizations below ", threshold, " kg/ha")
)


#### Spatial aggregation: explore yield per administrative division ####

## As reminder: we have these divisions 
plot(ra_avg_yield)
lines(sppol_divisions)
text(sppol_divisions, 
     labels=division_names, 
     cex=0.5, 
     col="black")


## Calculate mean yield per division based on mean yield map
df_mean_per_division <- cbind(division_names, 
  extract(ra_avg_yield, 
          sppol_divisions, 
          fun="mean", 
          na.rm=TRUE)
)

df_mean_per_division

sppol_divisions$mean_yield <- round(df_mean_per_division$mean)
plot(sppol_divisions, col=hcl.colors(8),
     y = "mean_yield",
     main = "Mean yield over 100 realizations, per division")


##  Mean yield per realization per region 
df_mean_per_region_per_realisation <- cbind(division_names,
                            extract(rs_yield_realisations, 
                                    sppol_divisions, 
                                    fun="mean", 
                                    na.rm=TRUE)
                            )

df_mean_per_region_per_realisation[, 1:9]

## What is the yield sd over the realizations, per region?

df_yield_sd_per_region <- 
  data.frame(division = division_names, 
            yield_sd = 
    apply(X = df_mean_per_region_per_realisation[, -(1:2)], ## apply over all columns except the first two 
          MARGIN = 1, ## apply over all columns, per row
          FUN = "sd"
          )
  )

sppol_divisions$yield_sd <- round(df_yield_sd_per_region$yield_sd)
plot(sppol_divisions, 
     y = "yield_sd", col=hcl.colors(25), type = "continuous",
     main = "SD yield over 100 realizations, per division")



## What is the yield prediction interval over the realizations, per region?

diff(  quantile(x = 10:20, probs = c(0.05, 0.95), names = FALSE) )

df_prediction_interval_per_region <- 
  data.frame(division = division_names, 
             prediction_interval = 
               apply(X = df_mean_per_region_per_realisation[, -(1:2)], ## apply over all columns except the first two 
                     MARGIN = 1, ## apply over all columns, per row
                     FUN =  function(x) {diff( 
                       quantile(x, probs = c(0.05, 0.95), names = FALSE,na.rm = TRUE) 
                     )}
               )
  )


sppol_divisions$prediction_interval <- round(df_prediction_interval_per_region$prediction_interval)

plot(sppol_divisions, 
     y = "prediction_interval", col=hcl.colors(8),
     main = "90% mean yield PIW, per division")


#### Spatial aggregation: summarise yield over whole area ####


## Average per realisation (represented as a vector of numbers)
vn_average_yield_per_realisation <- unlist(global(rs_yield_realisations, mean, 
                                                  na.rm=TRUE))

## Average over whole area

mean(vn_average_yield_per_realisation)

## sd over the averaged realizations
sd(vn_average_yield_per_realisation)

## Prediction interval width for the whole area 
n_yield_pred_int_width_whole_area <-  diff(  
                                       quantile(x = vn_average_yield_per_realisation, 
                                                probs = c(0.05, 0.95), 
                                                names = FALSE)
                                       )


#### Compare yield prediction interval over different aggregation levels ####

## For a nice visualization, we show the maps with the prediction interval 
## on pixel level, on division level, and for the whole area.

## Prepare 2 rasters, just for plotting with the panel function
ra_yield_pred_interv_per_division <- 
  rasterize(sppol_divisions, ra_yield_pred_interv, field="prediction_interval")

sppol_divisions$prediction_interval_global <- 
                            n_yield_pred_int_width_whole_area
ra_yield_pred_interv_whole <- 
  rasterize(sppol_divisions, ra_yield_pred_interv, field="prediction_interval_global")

## Plot combined, using same color scale
rs_plot <- c(ra_yield_pred_interv,
       ra_yield_pred_interv_per_division,
       ra_yield_pred_interv_whole
       )
names(rs_plot) <- c("Pixel support 90% PIW", 
              "Division support 90% PIW", 
              "Global support 90% PIW")

panel(rs_plot, nc = 3, nr = 1, fun=\()lines(sppol_divisions), box=TRUE)
## Please note: the "panel" function was quite recently added to the terra package.
## If it doesn't work in your case, update terra or run the following code:

par(mfrow=c(1,3)) 
plot(ra_yield_pred_interv, main = "Pixel support 90% PIW", range = c(75, 1225))
plot(ra_yield_pred_interv_per_division, main = "Division support 90% PIW",
     range = c(75, 1225), type = "continuous")
  lines(sppol_divisions)
plot(ra_yield_pred_interv_whole, main = "Global support 90% PIW", 
     range = c(75, 1225), type = "continuous")
par(mfrow=c(1,1))
