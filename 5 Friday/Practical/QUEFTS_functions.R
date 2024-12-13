## Tailor made function to work with Friday_uncertainty_spatial_aggregation.R
## Luc Steinbuch, December 2024
calculate_yield <- function(x)
{
  
  if(all(!is.na(x)))
  {
    soil_prop$pH <- x[1]
    soil_prop$SOC <- x[2]
    soil_prop$olsenP <- x[3]
    soil_prop$exchK <- x[4]
    
    tmp <- soil_supply(param = soil_param, soil = soil_prop, temp = settings$AvgTemp)
    soilN <- tmp[[1]]
    soilP <- tmp[[2]]
    soilK <- tmp[[3]]
    
    yield_value <- quefts(
      param = crop_param,
      supply = list(SN = soilN + fert_n, SP = soilP + fert_p, SK = soilK + fert_k),
      max_yield = settings$waterlimited_yield
    )[[2]]
    
    return(yield_value)
  } else
    return(NA) 
  
  
}














# Soil supply of nutrients as calculated by QUEFTS. The version of QUEFTS used 
# here incorporates the suggestions and additions by Sattari 2014
soil_supply <- function (param, soil, temp) 
{
  soil$Temp<-temp
  with(c(param, soil),
       {
         # pH correction factors for N
         if (pH < 4.7) 
         {
           fN <- 0.4
         } 
         else if (pH >= 4.7 && pH < 7) 
         {
           fN <- 0.25 * (pH - 3)
         } 
         else if (pH >= 7) 
         {
           fN <- 1
         }
         
         # pH correction factors for P
         if (pH < 4.7) 
         {
           fP <- 1
         } 
         else if (pH >= 4.7 && pH < 6) 
         {
           fP <- 1 - 0.5 * (pH - 6)^2 
         } 
         else if (pH >= 6 && pH < 6.7) 
         {
           fP <- 1
         } 
         else if (pH >= 6.7 && pH < 8) 
         {
           fP <- 1 - 0.5 * (pH - 6.7)^2
         } 
         else if (pH >= 8) 
         {
           fP <- 0.57
         }
         
         # pH correction factors for K
         if (pH < 6.8) 
         {
           fK <- 6.1 * pH^(-1.2)
         } 
         else if (pH >= 6.8) 
         {
           fK <- 0.6
         }
         
         # Scaling parameters for QUEFTS-equations
         if (is.na(Temp)){ 
           Temp <- 25
         }
         alfaN <- 2 * 2^((Temp - 9) / 9)
         

         # Potential N (1), P (2) and K (3) supply - QUEFTS functions
         soilN <- fN * alfaN * SOC                               #(1) -> kg N ha-1 season-1
         soilP <- fP * alfaP * SOC + betaP * olsenP              #(2) -> kg P ha-1 season-1
         soilK <- (fK * alfaK * exchK) / (gammaK + betaK * SOC)  #(3) -> kg K ha-1 season-1
         
         return(list(soilN, soilP, soilK))
       })
}

# Nutrient uptake depends on the soil supply of the nutrient and the supply of 
# other nutrients
nutrient_uptake <- function(S1=NA, S2=NA, d1=NA, a1=NA, d2=NA, a2=NA, r1=NA, 
                            r2=NA)
{
  # N, P and K uptakes based on QUEFTS
  if (S1 < r1 + ((S2 - r2) * a2 / d1))
  {
    uptakeX_givenY = S1
  }
  else if (S1 > r1 + ((S2 - r2) * (2 * d2 / a1 - a2 / d1)))
  {
    uptakeX_givenY = r1 + (S2 - r2) * (d2 / a1)
  }
  else
  {
    uptakeX_givenY = S1 - 0.25 * (S1 - r1 - (S2 - r2) * (a2 / d1))**2 / ((S2 - r2) * (d2 / a1 - a2 / d1))
  }
  # Nutrient uptake given availability of other nutrient
  return(uptakeX_givenY)
}

# Actual uptake of nutrients as determined by the availability of the other 
# nutrients
actual_uptake <- function(supply, param)
{
  with(c(supply, param),
       {
         UNP <- nutrient_uptake(S1 = SN, S2 = SP, d1 = dN, a1 = aN, d2 = dP, a2 = aP, r1 = rN, r2 = rP)
         UNK <- nutrient_uptake(S1 = SN, S2 = SK, d1 = dN, a1 = aN, d2 = dK, a2 = aK, r1 = rN, r2 = rK)
         UN <- min(UNP, UNK)
         UPN <- nutrient_uptake(S1 = SP, S2 = SN, d1 = dP, a1 = aP, d2 = dN, a2 = aN, r1 = rP, r2 = rN)
         UPK <- nutrient_uptake(S1 = SP, S2 = SK, d1 = dP, a1 = aP, d2 = dK, a2 = aK, r1 = rP, r2 = rK)
         UP <- min(UPN, UPK)
         UKN <- nutrient_uptake(S1 = SK, S2 = SN, d1 = dK, a1 = aK, d2 = dN, a2 = aN, r1 = rK, r2 = rN)
         UKP <- nutrient_uptake(S1 = SK, S2 = SP, d1 = dK, a1 = aK, d2 = dP, a2 = aP, r1 = rK, r2 = rP)
         UK <- min(UKN, UKP)
         
         return(c(UN, UP, UK))
       })
}

# Maximum and minimum yields based on the uptake of nutrients and maximum and minmum use efficiences
max_min_yields <- function(param, uptake) 
{
  with(c(param, uptake),
       {
         YNA <- (UN - rN) * aN
         YND <- (UN - rN) * dN
         YPA <- (UP - rP) * aP
         YPD <- (UP - rP) * dP
         YKA <- (UK - rK) * aK
         YKD <- (UK - rK) * dK
         
         return(c(YNA, YND, YPA, YPD, YKA, YKD))
       })
}

# Grain yield as determined by the availability of pairs of nutrients.
yield_nutrients_combined <- function(U1=NA, d1=NA, a1=NA, Y2A=NA, Y2D=NA, Y3D=NA, r1=NA)
{
  # Yield calculated based on the combined uptake of 2 nutrients, while take into account
  # the availability of the third nutrient.  
  
  # Determine which nutrient limited yield is lowest.
  YxD = min(Y2D, Y3D)
  
  # If the uptake of one of the nutrients, and therefore the yield associated with that 
  # nutrient, is zero the overall yield is also zero.
  if (U1 == 0 || YxD == 0) 
  {
    Y12 = 0
  }
  else
  {
    Y12 = Y2A + (2 * (YxD - Y2A) * (U1 - r1 - Y2A / d1)) / (YxD / a1 - Y2A / d1) -
      (YxD - Y2A) * (U1 - r1 - Y2A / d1)**2 / (YxD / a1 - Y2A / d1)**2
  }
  # Return the calculated yield based on the uptake of nutrients 1 and 2
  return(Y12)
}

# Grain yield calculated based on the availability of pairs of nutrients
final_yield <- function(param, uptake, yields, max_yield = max_yield)
{
  with(c(param, uptake, yields), 
       {
         YNP <- yield_nutrients_combined(U1 = UN, d1 = dN, a1 = aN, Y2A = YPA, Y2D = YPD, Y3D = YKD, r1 = rN)
         YNK <- yield_nutrients_combined(U1 = UN, d1 = dN, a1 = aN, Y2A = YKA, Y2D = YKD, Y3D = YPD, r1 = rN)
         YPN <- yield_nutrients_combined(U1 = UP, d1 = dP, a1 = aP, Y2A = YNA, Y2D = YND, Y3D = YKD, r1 = rP)
         YPK <- yield_nutrients_combined(U1 = UP, d1 = dP, a1 = aP, Y2A = YKA, Y2D = YKD, Y3D = YND, r1 = rP)
         YKN <- yield_nutrients_combined(U1 = UK, d1 = dK, a1 = aK, Y2A = YNA, Y2D = YND, Y3D = YPD, r1 = rK)
         YKP <- yield_nutrients_combined(U1 = UK, d1 = dK, a1 = aK, Y2A = YPA, Y2D = YPD, Y3D = YND, r1 = rK)
         #first yield estimate
         YE <- mean(c(YNP, YNK, YPN, YPK, YKN, YKP))
         
##TS: 21DEc2016 adjusted and added max_yield, otherwise a linear response curve will be found....

         #check for minimum conditions
         YNPc <- min(c(YNP,YND,YPD,YKD,max_yield))
         YNKc <- min(c(YNK,YND,YPD,YKD,max_yield))
         YPNc <- min(c(YPN,YND,YPD,YKD,max_yield))
         YPKc <- min(c(YPK,YND,YPD,YKD,max_yield))
         YKNc <- min(c(YKN,YND,YPD,YKD,max_yield))
         YKPc <- min(c(YKP,YND,YPD,YKD,max_yield))
         #Final estimate
         YEc <- mean(c(YNPc, YNKc, YPNc, YPKc, YKNc, YKPc))
         
         return(YEc)         
       })
}

# Nutrient uptake and grain yield as calculated by QUEFTS
quefts <- function(param, supply, max_yield)
{
  # Actual uptake of nutrients.
  tmp <- actual_uptake(supply=supply, param=param)
  UN <- tmp[[1]]
  UP <- tmp[[2]]
  UK <- tmp[[3]]
  uptake <- list(UN=UN, UP=UP, UK=UK)
  
  # Maximum and minimum yields, depending on maximum accumulation and dilution.
  tmp <- max_min_yields(par=param, upt=uptake)
  YNA <- tmp[[1]]
  YND <- tmp[[2]]
  YPA <- tmp[[3]]
  YPD <- tmp[[4]]
  YKA <- tmp[[5]]
  YKD <- tmp[[6]]
  yields <- list(YNA=YNA, YND=YND, YPA=YPA, YPD=YPD, YKA=YKA, YKD=YKD)
  
  # Final yield based on the combinations of nutrient uptake and minimum + maximum yields.
  yield <- final_yield(param=param, uptake=uptake, yields=yields, max_yield = max_yield)
  
  return(list(uptake, yield))
}


quefts_soil_parameters <- function(mydata)
{
  for (i in 1:nrow(mydata)){
  pH <- mydata[i,'pH']
  SOC <- mydata[i,'SOC']
  olsenP <- mydata[i,'olsenP']
  exchK <- mydata[i,'exchK']
  }
  
  soil_prop <- list(olsenP=olsenP, exchK=exchK, pH=pH, SOC=SOC)
  
  # Parameter values for the soil supply functions in QUEFTS. If temperature is given the
  # relation described in Sattari et al. 2014 is used for alfaN.
  temp <- 25
  alfaN <- 6.8 #default = 6.8
  alfaP <- 0.35 #default = 0.35
  betaP <- 0.5
  alfaK <- 400 #default = 400
  betaK <- 0.9
  gammaK <- 2
  soil_param <- list(temp=temp, alfaN=alfaN, alfaP=alfaP, betaP=betaP, alfaK=alfaK, betaK=betaK, gammaK=gammaK)
  
  return(list(soil_prop, soil_param))
}

quefts_crop_parameters <- function(target_yield=NA, crop_id=NA, crop_pars=NA)
{
  current_crop <- crop_pars[crop_pars$crop_id == crop_id, ]
  
  # Nutrient use efficiencies associated with maximum accumulation (a) and dilution (d) 
  # in kg DM per kg nutrient
  aN <- current_crop$aN
  dN <- current_crop$dN
  aP <- current_crop$aP
  dP <- current_crop$dP
  aK <- current_crop$aK
  dK <- current_crop$dK
  # Nutrient uptake at zero grain yield
  rN <- current_crop$rN
  rP <- current_crop$rP
  rK <- current_crop$rK
  tolerance <- 0.01
  crop_param <- list(aN=aN, dN=dN, aP=aP, dP=dP, aK=aK, dK=dK, rN=rN, rP=rP, rK=rK, tolerance=tolerance)
  
  # Target yield to calculate nutrient requirements for
  yield <- 0
  yields <- list(yield=yield, target_yield=target_yield)
  
  # Minimum and maximum applications of N and the mean of these two values.
  minSN <- target_yield / dN
  maxSN <- target_yield / aN
  medSN  <- mean(c(maxSN, minSN))
  
  return(list(crop_param, yields, minSN, maxSN, medSN))
}

