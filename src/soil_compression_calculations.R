

# ======================================================
# load packages

library(PeriodicTable)
library(raster)
library(rgdal)

# ======================================================
# Load input data

# set working directory
# outputs from the spatial processing script should go here
# enter the directory and uncomment line of code below
# setwd("C:/users/.../spatial_output")
setwd("C:/Users/Eric/Dropbox/C_initiative/SOC_critique/SOC_crit_CA_geospatial/June_upload")

# sampled points
spts <- read.csv("pts_dat_US_060221.csv")

# subset of the value1 table from cultivated land
ds <- read.csv("gSSURGO_US_value1_subset.csv")

# exclude map units without carbon data
ds <- ds[!is.na(ds$soc0_5) & !is.na(ds$soc0_20) & !is.na(ds$soc0_30) & !is.na(ds$soc0_100), ]

# exclude sample points not represented by map units
spts <- spts[spts$mu %in% ds$mukey, ]



# ======================================================
# function for fitting a cumulative depth curve
fit_cdp <- function(mukey, ds) {
  # subset data for the map unit
  mudata <- ds[ds$mukey == mukey, ]

  # Cumulative OC stocks for the map unit (g m-2)
  muC <- as.numeric(c(0, mudata[, c("soc0_5", "soc0_20", "soc0_30", "soc0_100", "soc0_150")]))

  # SSURGO depth intervals (cm)
  muD <- c(0, 5, 20, 30, 100, 150)

  # fit a spline function to get the cumulative depth profile (output in g m-2)
  cdp0 <- splinefun(muD, muC, method = "hyman")
  return(cdp0)
}

# ======================================================
# change_BD
# Function to simulate the effect of changing A-horizon bulk density on carbon depth profiles
# Inputs:
# D_eval = depths to be evaluated (cm)
# D_A0 = original depth of A horizon (cm)
# dBD = change in bulk density (%)
# cdp0 = fitted spline function for  carbon depth profile (returns cumulative C at D_eval)

# (1) output is only defined when the A horizon depth is less than or equal to the evaluation depth
# (2) change in bulk density must be positive (an increase in density)
# (3) the function calculates an "effective sampling depth" after an increase in density
# (4) this depth is in the reference frame of the initial carbon depth profile
# (5) the carbon stock in the additional increment is returned as an output

change_BD <- function(D_A0, dBD, cpd0, D_eval = 30) {
  if (D_eval < D_A0 | dBD < 0) {
    cat("\n\nD_eval < D_A0 or dBD < 0; returning NA")
    return(NA)
  }

  # C0 = original cumulative carbon as a function of sampling depth (g m-2)
  C0 <- cdp0(D_eval)

  # BDf = scaling factor for bulk density (unitless)
  BDf <- ((100 + dBD) / 100)

  # D_A1 = new A horizon depth (cm)
  D_A1 <- D_A0 / BDf

  # dD = difference in A horizon depths (cm)
  dD <- D_A0 * (dBD / 100) # D_A0 - D_A1

  # evaluate the spline function to get new cumulative C with depth (g m-2)
  C1 <- cdp0(D_eval + dD)

  # precompute the difference (g m-2)
  dC <- C1 - C0

  # convert units to MT CO2e ha-1

  dC <- dC / 100 * ((mass("C") + 2 * mass("O")) / mass("C"))
  return(dC = dC)
}

# ======================================================
# preallocate data frame for results
mukeys <- unique(ds$mukey)
lmukey <- length(mukeys)
result <- data.frame(
  mu = mukeys,
  dC_10_1 = rep(NA, lmukey),
  dC_10_2 = rep(NA, lmukey),
  dC_10_3 = rep(NA, lmukey),
  dC_10_4 = rep(NA, lmukey),
  dC_10_5 = rep(NA, lmukey),
  dC_20_1 = rep(NA, lmukey),
  dC_20_2 = rep(NA, lmukey),
  dC_20_3 = rep(NA, lmukey),
  dC_20_4 = rep(NA, lmukey),
  dC_20_5 = rep(NA, lmukey),
  dC_30_1 = rep(NA, lmukey),
  dC_30_2 = rep(NA, lmukey),
  dC_30_3 = rep(NA, lmukey),
  dC_30_4 = rep(NA, lmukey),
  dC_30_5 = rep(NA, lmukey)
)

# loop through map units

for (i in 1:lmukey) {
  mukey <- mukeys[i]
  cdp0 <- fit_cdp(mukey, ds)

  result$dC_10_1[i] <- change_BD(10, 1, cdp0) # MT CO2e ha-1
  result$dC_10_2[i] <- change_BD(10, 2, cdp0) # MT CO2e ha-1
  result$dC_10_3[i] <- change_BD(10, 3, cdp0) # MT CO2e ha-1
  result$dC_10_4[i] <- change_BD(10, 4, cdp0) # MT CO2e ha-1
  result$dC_10_5[i] <- change_BD(10, 5, cdp0) # MT CO2e ha-1

  result$dC_20_1[i] <- change_BD(20, 1, cdp0) # MT CO2e ha-1
  result$dC_20_2[i] <- change_BD(20, 2, cdp0) # MT CO2e ha-1
  result$dC_20_3[i] <- change_BD(20, 3, cdp0) # MT CO2e ha-1
  result$dC_20_4[i] <- change_BD(20, 4, cdp0) # MT CO2e ha-1
  result$dC_20_5[i] <- change_BD(20, 5, cdp0) # MT CO2e ha-1

  result$dC_30_1[i] <- change_BD(30, 1, cdp0) # MT CO2e ha-1
  result$dC_30_2[i] <- change_BD(30, 2, cdp0) # MT CO2e ha-1
  result$dC_30_3[i] <- change_BD(30, 3, cdp0) # MT CO2e ha-1
  result$dC_30_4[i] <- change_BD(30, 4, cdp0) # MT CO2e ha-1
  result$dC_30_5[i] <- change_BD(30, 5, cdp0) # MT CO2e ha-1

  pct <- as.character(round(100 * i / (lmukey), 1))
  pct <- paste(pct, paste(rep(" ", 4 - nchar(pct)), collapse = ""), sep = "")
  cat("\r\t", pct, "% complete\n\b")
}

# merge data by map unit
spts <- merge(spts, result, by = "mu")

# mean aggregate by 50 km cell
cellmeans <- aggregate(spts[, grepl("dC", colnames(spts))], list(cellID = spts$cellID), mean)

# read the cell summary statistics table
cells <- read.csv("pts_grid_US_060221.csv")

# calculate the hectares of cultivated land in each cell
cells$farm_hectares <- cells$f_cult * 50^2 * 100

# merge/append the cell mean bulk density effect (MT CO2e ha-1)
cells <- merge(cells, cellmeans, by = "cellID")

# save the result
write.csv(cells, "cell_data.csv", row.names = F)

# Ag sector emissions in 2018:
# https://www.epa.gov/sites/production/files/2020-04/documents/us-ghg-inventory-2020-main-text.pdf
Ag_emissions <- 618.5

# calculate the total amount of carbon susceptible to the bulk density effect, in MMT CO2e
BD_effect <- sum(cells$farm_hectares * cells$dC_20_3 * 1e-6, na.rm = T)

# the bulk density effect as a percentage of 2018 agriculture emissions
pct_BD_effect <- BD_effect / Ag_emissions * 100
