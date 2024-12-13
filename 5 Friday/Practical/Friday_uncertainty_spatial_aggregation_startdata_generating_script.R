#### Preparation ####

## Running time of this script: ca 2 1/2 minute (with nmax = 25 and nreal = 100)
## It saves at the end the file df_realisations_soil_prop.Rda (25 MB with above settings)

## Load needed packages 
library(gstat)
library(terra)

## Set working directory to the same directory as this script. Works only in RStudio
setwd(
  dirname(
    rstudioapi::getActiveDocumentContext()$path
    )
  )

## Make the pseudo-randomness which is used in several functions repeatable
set.seed(12345)

## Define the Coordinate Reference System we will be using for all spatial objects
target_crs <- "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +units=m +ellps=WGS84 +datum=WGS84"

## Read Random Forest maps and residuals from RData file
load("trend_part_of_DSM_test_train_split_LOSCV.RData")

## Which objects have just been loaded?
ls()  
## There are several files we will not use in this script.
## But we have the prediction maps of the 4 soil prperties,
## and the residual maps of our observations vs. the same
## locations on those those 4 maps.

#### Correlation between soil property residuals ####
## To have some idea if there are spetial relationships between
## the soil properties
residuals_df <- cbind(residuals_pH_df, Predicted_Carbon = residuals_Carbon_df$Predicted, Residuals_Carbon = residuals_Carbon_df$Residuals_Carbon, 
                      Predicted_OlsenP = residuals_olsenP_df$Predicted, Residuals_OlsenP = residuals_olsenP_df$Residuals_olsenP, 
                      Predicted_ExchK = residuals_ExchK.mmolkg_df$Predicted, Residuals_ExchK = residuals_ExchK.mmolkg_df$Residuals_ExchK.mmolkg)

head(residuals_df)

correlation_data <- residuals_df[, c("Residuals_pH", "Residuals_Carbon", "Residuals_OlsenP", "Residuals_ExchK")]
colnames(correlation_data) <- c("pH", "SOC", "Olsen P", "Exch K")

correlation_matrix <- cor(correlation_data, use = "pairwise.complete.obs")

print(round(correlation_matrix, 4))

#### Spatial modelling of residuals ####
#### First fit the variogram models for the 4 soil properties separatly

##  Variogram model for pH ##
vg_pH <- variogram(Residuals_pH ~ 1, 
                   locations = ~Latitude+Longitude,
                   data = residuals_df)

vgm_model_pH <- vgm(nugget = 0.02, psill = 0.2, range = 18000, model = "Sph")
vgm_model_fit_pH <- fit.variogram(vg_pH, vgm_model_pH)
vgm_model_fit_pH
plot(vg_pH, pl = T, model = vgm_model_fit_pH, main = "Fitted variogram for pH")
 
##  Variogram model for SOC ##
vg_soc <- variogram(Residuals_Carbon ~ 1, 
                    locations = ~Latitude+Longitude,
                    data = residuals_df)

vgm_model_soc <- vgm(nugget = 23, psill = 50, range = 10000, model = "Sph")
vgm_model_fit_soc <- fit.variogram(vg_soc, vgm_model_soc)
vgm_model_fit_soc
plot(vg_soc, pl = T, model = vgm_model_fit_soc, main = "Fitted variogram for SOC") 

##  Variogram model for Olsen P ##
vg_olsenp <- variogram(Residuals_OlsenP ~ 1, 
                       locations = ~Latitude+Longitude,
                       data = residuals_df)
 
vgm_model_olsenp <- vgm(nugget = 2, psill = 8, range = 20000, model = "Sph")
vgm_model_fit_olsenp <- fit.variogram(vg_olsenp, vgm_model_olsenp)
 vgm_model_fit_olsenp
plot(vg_olsenp, 
     pl = T, 
     model = vgm_model_fit_olsenp,
     main = "Fitted variogram for P")
 
 
##  Variogram model for Exch K ##
vg_exchk <- variogram(Residuals_ExchK ~ 1, 
                      locations = ~Latitude+Longitude, 
                      data = residuals_df)
vgm_model_exchk <- vgm(nugget = 0.000014, psill = 0.00001 , range = 12000, model = "Sph")
vgm_model_fit_exchk <- fit.variogram(vg_exchk, vgm_model_exchk)
vgm_model_fit_exchk
 
plot(vg_exchk, vgm_model_fit_exchk,
     main = "Fitted variogram for K")


##  Variogram model for co-kriging ##
## We use above models as starting values
nmax <- 25


g <- gstat(NULL, id = "pH", 
           nmax=nmax, 
           form = Residuals_pH ~ 1, 
           locations = ~Latitude+Longitude,
           model = vgm_model_fit_pH,
           data = residuals_df)
g <- gstat(g, id = "SOC", 
           nmax=nmax, 
           locations = ~Latitude+Longitude,
           model = vgm_model_fit_soc,
           form = Residuals_Carbon ~ 1, 
           data = residuals_df)
g <- gstat(g, id = "OlsenP", 
           nmax=nmax, 
           locations = ~Latitude+Longitude,
           form = Residuals_OlsenP ~ 1, 
           model = vgm_model_fit_olsenp,
           data = residuals_df)
g <- gstat(g, 
           id = "ExchK", 
           nmax=nmax, 
           locations = ~Latitude+Longitude,
           form = Residuals_ExchK ~ 1, 
           model = vgm_model_fit_exchk,
           data = residuals_df)
g




# g_fitted <- gstat(g, id = c("pH"),  model = vgm_model_fit_pH, fill.all = T)
# #g_fitted <- gstat(g, id = "pH", model = g$model, fill.all = T)
# g_fitted

## Linear Model of Coregionalization
v_auto_cross <- variogram(g)
g_fitted <- fit.lmc(v_auto_cross, g, vgm_model_fit_pH, fit.method = 6, correct.diagonal = 1.01)
g_fitted

plot(variogram(g_fitted), 
     model = g_fitted$model,
     main = "Fitted variograms and co-variagrams for all soil properties")


#### Spatial prediction of residuals for all variables by co-kriging ####

## Load prediction raster from file...
pred_raster <- readRDS("original_raster_grid.rds")
## .. and turn it into a dataframe with the coordinates rather 
## than a terra raster object.
df_pred_grid <- as.data.frame(pred_raster, xy = TRUE, cells = TRUE)[2:3]

head(df_pred_grid)

names(df_pred_grid) <- c("Longitude", "Latitude")

## Set number of soil property realisations
nreal <- 100

df_realisations_soil_prop <- predict(g_fitted, 
                            newdata = df_pred_grid,
                            nsim=nreal,
                            debug.level = -1 ## For indicating progress
                            )
dim(df_realisations_soil_prop)
df_realisations_soil_prop[1:5, 1:8]

## Note that we now have modelled the difference with the RF map. 
## We need to add those differences to the RF map, per soil property.

## Add RF map pH
indices_columns_pH <- 3:(nreal+2)
names(df_realisations_soil_prop)[indices_columns_pH]

## Just to make sure that we fill the dataframe columnwise,
## we create a matrix:
ma_pH <- matrix(data = predictions_map_df_pH$Preds_pH, 
                nrow = nrow(df_realisations_soil_prop),
                ncol = nreal,
                byrow = FALSE)
ma_pH[1:5, 1:5]

df_realisations_soil_prop[,indices_columns_pH] <- 
  df_realisations_soil_prop[,indices_columns_pH] + ma_pH

## Add RF map SOC

indices_columns_SOC <- (1*nreal+3):(1*nreal+102)
names(df_realisations_soil_prop)[indices_columns_SOC]

ma_SOC <- matrix(data = predictions_map_df_Carbon$Preds_Carbon, 
                nrow = nrow(df_realisations_soil_prop),
                ncol = nreal,
                byrow = FALSE)
ma_SOC[1:5, 1:5]

df_realisations_soil_prop[,indices_columns_SOC] <- 
  df_realisations_soil_prop[,indices_columns_SOC] + ma_SOC

df_realisations_soil_prop[1:5, 101:108]

## Add RF map Olsen P
df_realisations_soil_prop[1:5, 205:208]

indices_columns_OlsenP <- (2*nreal+3):(2*nreal+102)
names(df_realisations_soil_prop)[indices_columns_OlsenP]

ma_OlsenP <- matrix(data = predictions_map_df_olsenP$Preds_olsenP, 
                 nrow = nrow(df_realisations_soil_prop),
                 ncol = nreal,
                 byrow = FALSE)
ma_OlsenP[1:5, 1:5]

df_realisations_soil_prop[,indices_columns_OlsenP] <- 
  df_realisations_soil_prop[,indices_columns_OlsenP] + ma_OlsenP

df_realisations_soil_prop[1:5, 205:208]


names(df_realisations_soil_prop)

## Add RF map Exch K

df_realisations_soil_prop[1:5, 305:308]

indices_columns_ExchK <- (3*nreal+3):(3*nreal+102)
names(df_realisations_soil_prop)[indices_columns_ExchK]

ma_ExchK <- matrix(data = predictions_map_df_ExchK.mmolkg$Preds_ExchK.mmolkg, 
                    nrow = nrow(df_realisations_soil_prop),
                    ncol = nreal,
                    byrow = FALSE)
ma_ExchK[1:5, 1:5]

df_realisations_soil_prop[,indices_columns_ExchK] <- 
  df_realisations_soil_prop[,indices_columns_ExchK] + ma_ExchK

df_realisations_soil_prop[1:5, 305:308]

## Remove negative values (because they have no meaning),
## except for the coordinates

## How many values < 0?
sum(df_realisations_soil_prop[, -(1:2)] <  0)

## Matrix of TRUE where < 0
ma_indices_smaller <- (df_realisations_soil_prop[, -(1:2)] < 0)

## Add first two columns for coordinates
ma_indices_smaller <- cbind(FALSE, FALSE, ma_indices_smaller)

## Check
dim(ma_indices_smaller)
ma_indices_smaller[1:5, 1:5]

## Set values < 0 to 0
df_realisations_soil_prop[ma_indices_smaller] <- 0

## Check: How many values < 0?
sum(df_realisations_soil_prop[, -(1:2)] <  0)

## For some reason, the rast function in the terra package wants the first column to
## be the longitude (x), and the second latitude (y), while
## gstat returs it otherwise, so we need to swap:

df_realisations_soil_prop[1:3, 1:3]

df_realisations_soil_prop[c(1,2)] <- df_realisations_soil_prop[c(2,1)]
names(df_realisations_soil_prop)[c(1,2)] <- names(df_realisations_soil_prop)[c(2,1)]

df_realisations_soil_prop[1:3, 1:3]

## Save the result for use in the script "Friday_uncertainty_spatial_aggregation.R"
save(df_realisations_soil_prop, file= "df_realisations_soil_prop.Rda")

