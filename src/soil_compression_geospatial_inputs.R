# Code for extracting SSURGO map units and cultivation status in a 50 km grid across the USA
# outputs:
# (1) grid cell centers with fraction cultivated
# (2) map units associated with cultivated points sampled from each cell

# ======================================================
# load packages
library(rgdal)
library(raster)

# ======================================================
# read gSSURGO data
# this must be downloaded and stored locally
# url: https://gdg.sc.egov.usda.gov/
# put these files in the same directory

# ***user input***
# substitute directory here and uncomment this line of code:
# setwd("C:/.../gSSURGO_CONUS")

# read the raster file that records map units
mu <- raster("gSSURGO_CONUS.tif")

# read the "value1" lookup table that links soil properties with map units
v1 <- readOGR("gSSURGO_CONUS.gdb", "Valu1", dropNULLGeometries = F, stringsAsFactors = F)

# subset the SOC stock data from value1
v1 <- v1[, colnames(v1) == "mukey" | grepl("soc", colnames(v1))]

# ======================================================
# read the 2020 cropland data layer (2016-2020)
# this must be downloaded and stored locally
# url: https://www.nass.usda.gov/Research_and_Science/Cropland/Release/index.php


# ***user input***
# substitute directory here and uncomment this line of code:
# setwd("C:/.../2020_Cultivated_Layer")
cr <- raster("2020_Cultivated_Layer.img")

# ======================================================
# read shapefile of US counties
# this must be downloaded and stored locally
# url: https://www.census.gov/geographies/mapping-files/time-series/geo/carto-boundary-file.html

# ***user input***
# substitute directory here and uncomment this line of code:
# setwd("C:/.../cb_2018_us_county_20m")
cts <- readOGR("cb_2018_us_county_20m.shp", "cb_2018_us_county_20m")

# reproject shapefil of counties into the CRS of the cropland data layer
cts <- spTransform(cts, crs(cr))
cts <- crop(cts, extent(cr))

# ======================================================
# create 50 km grid template for the conterminous USA
US50km_template <- raster(res = 50000, ext = extent(cts), crs = crs(cts))
US50km <- rasterize(cts, US50km_template)
USgrid <- rasterToPolygons(US50km)
USgrid$cellID <- 1:nrow(USgrid)

# attach lat/lon centroids for each grid cell
USgrid <- spTransform(USgrid, CRS("+proj=longlat +datum=WGS84"))
library(rgeos)
USpts <- gCentroid(USgrid, byid = T)
USgrid@data <- cbind(USgrid@data, USpts@coords)



# ======================================================
# preallocate fields

# commodity crops
USgrid$f_cult <- rep(NA, nrow(USgrid))

# ======================================================
# loop through 50 km cells, sampling points to extract data

n_cells <- nrow(USgrid)
n_samp <- 250

for (i in 1:n_cells) {
  cat("\n", i, "\n")
  cell <- USgrid[i, ]
  cell <- SpatialPolygons(cell@polygons, proj4string = crs(cell))
  ex <- extent(cell)

  # Take a random sample of n points from the 50 km cell
  rs <- SpatialPoints(data.frame(
    x = sample(ex[1]:ex[2], n_samp, replace = T),
    y = sample(ex[3]:ex[4], n_samp, replace = T)
  ),
  proj4string = crs(cell)
  )

  # extract the cultivated layer to get fraction of cell cultivated
  crp <- extract(cr, rs)
  USgrid$f_cult[i] <- mean(crp == 2, na.rm = T)

  # only proceed with points that are cultivated
  rs <- rs[crp == 2 & !is.na(crp)]

  # reproject the random sample to match the gSSURGO data
  rstr <- spTransform(rs, crs(mu))

  # extract gSSURGO map units
  mup <- extract(mu, rstr)

  # store data from individual sample points
  pt_dat <- data.frame(
    COUNTYNS = rep(USgrid$layer[i], length(mup)),
    cellID = rep(USgrid$cellID[i], length(mup)), mu = mup
  )
  if (!exists("pts_dat")) {
    pts_dat <- pt_dat
  } else {
    pts_dat <- rbind(pts_dat, pt_dat)
  }
}

# ======================================================
# write the outputs

# ***user input***
# create an output directory
# enter directory and uncomment this line of code
# setwd("C:/.../spatial_outputs")

# save the individual sample points (map units)
write.csv(pts_dat, "pts_dat_US_060221.csv", row.names = F)

# save the cell summary statistics (fraction cultivated)
write.csv(USgrid@data, "pts_grid_US_060221.csv", row.names = F)

# save a subset of the value1 lookup table that includes only the sampled map units
v2 <- v1[v1$mukey %in% pts_dat$mu, ]
write.csv(v2, "gSSURGO_US_value1_subset.csv", row.names = F)
