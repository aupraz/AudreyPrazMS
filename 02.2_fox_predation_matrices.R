setwd("C:/Users/dafcluster5/Box Sync/Audrey Praz")
library(raster)
library(sf)
library(parallel)
library(pbapply)
library(rnaturalearth)
library(ggplot2)
library(rgdal)
library(doSNOW)
library(patchwork)


# region from matrix
region <- readRDS("data/numbat_temporal_region.RDS")
str(region, max.level = 2) # attributes for raster 
templ_ras <- raster(res = 0.20,
                    ext = extent(111.05,155.05,-45.025,-9.025),
                    crs = crs(raster()))
values(templ_ras) <- 1
proj_ext <- projectExtent(templ_ras, crs = attr(region, "proj_wkt"))
res(proj_ext) <- attr(region, "res")
all.equal(coordinates(proj_ext)[,1],  attr(region, "X")) ## TRUE
all.equal(coordinates(proj_ext)[,2],  attr(region, "Y")) ## TRUE

# cont outlines
conts <- ne_countries(scale = 50, returnclass = "sf")
conts <- conts[conts$name == "Australia", ]
conts <- st_crop(conts, extent(111.05, 155.05, -45.025, -9.025))

TEST_REGION <- TRUE
if (TEST_REGION) {
  templ_rasprj <- projectRaster(templ_ras, proj_ext, method = "ngb")
  suppressWarnings({max_ext <- apply(region, 1, max, na.rm = TRUE)})
  templ_rasprj[] <- max_ext
  templ_rasprj[is.infinite(templ_rasprj)] <- NA
  plot(templ_rasprj,
       legend = FALSE,
       addfun = lines(as_Spatial(st_transform(conts, crs = st_crs(proj_ext)))))
}

##### Fox arrival ####
# Load the dataset representing fox arrival across time.
fox_abundance<- stack("data/animations/ensemble_fox_abundance.grd")
#temp_mask <- raster::raster("data/climate/template_full_20km.grd")
fox_abundance_stack <- stack()
for (i in 1:nlayers(fox_abundance)) {
  fox_abundance_stack <- addLayer(fox_abundance_stack, 
                                         projectRaster(fox_abundance[[i]], 
                                                       cat_arrival_int, 
                                                       method = "bilinear"))
}
dim(fox_abundance_stack)

# Define a threshold for presence/absence
threshold <- 0  # Set the threshold (e.g., >0 means presence)

# Convert each layer to binary presence/absence
fox_presence_absence <- calc(fox_abundance_stack, fun = function(x) { 
  x[x > threshold] <- 1  # Assign 1 for presence
  x[x <= threshold] <- 0 # Assign 0 for absence
  return(x)
})

# Save the resulting stack (optional)
writeRaster(fox_presence_absence, "data/fox_invasion_abundance/ensemble_fox_presence_absence.grd", format = "raster", overwrite = TRUE)

# Load the binarized dataset
fox_arrival <-  stack("data/fox_invasion_abundance/ensemble_fox_presence_absence.grd")
plot(fox_arrival)

# Copy the first layer and replicate it
first_layerfox <- fox_arrival[[1]]
values(first_layerfox)[values(first_layerfox) == 1] <- 0  # Turn 1s into 0s
prepend_stackfox <- stack(replicate(121, first_layerfox))

# Merge with existing stack
fox_arrival_int <- stack(prepend_stackfox, fox_arrival)

# convert fox to temporal arrival mask
names(fox_arrival_int) <- paste0("X_",seq(1750, by = 1, length = nlayers(fox_arrival_int)))
fox_arrival_int <- setZ(fox_arrival_int, seq(1750, by = 1, length = nlayers(fox_arrival_int)), "Years")
plot(fox_arrival_int[[seq(1, 271, length = 9)]], zlim = c(0, 1))


# Convert to matrix and save
fox_arrival_intMat <- as.matrix(fox_arrival_int)
#object.size(hv_rast_maskedMat)
mode(fox_arrival_intMat) <- "integer"
#object.size(hv_rast_maskedMat)
attr(fox_arrival_intMat, "timesteps") <- seq(1750, by = 1, length = ncol(fox_arrival_intMat))
attr(fox_arrival_intMat, "res") <- res(fox_arrival_int)[1]
attr(fox_arrival_intMat, "extent") <- "-2423364, 2436636, -2158270, 1981730"
attr(fox_arrival_intMat, "proj4") <- "+proj=laea +lat_0=-27 +lon_0=133 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
attr(fox_arrival_intMat, "proj_wkt") <- wkt(fox_arrival_int)
attr(fox_arrival_intMat, "scale.fac") <- NA
attr(fox_arrival_intMat, "X") <- coordinates(fox_arrival_int)[,1]
attr(fox_arrival_intMat, "Y") <- coordinates(fox_arrival_int)[,2]
str(fox_arrival_intMat)

sapply(names(attributes(fox_arrival_intMat)), function(x) {
  all.equal(attr(fox_arrival_intMat, x), attr(region, x))
})


saveRDS(fox_arrival_intMat, file = "data/fox_invasion_abundance/fox_arrival_temporal_mask.RDS", compress = TRUE)

#### Upper estimate ####
# original res in km
orig_res <- res(raster("data/fox_invasion_abundance/upper_bound.tif"))[1]/1000

upper <- raster("data/fox_invasion_abundance/upper_bound.tif")/orig_res
lower <- raster("data/fox_invasion_abundance/lower_bound.tif")/orig_res
mean <- raster("data/fox_invasion_abundance/estimate.tif")/orig_res
abundance <- stack(lower, mean, upper, quick = TRUE)
temp_mask <- templ_rasprj
temp_mask[is.na(temp_mask) | temp_mask == 0] <- 0
abundance <- aggregate(abundance, fact = res(fox_arrival_int)/res(abundance), 
                       fun = max, na.rm = TRUE)
abundance <- mask(projectRaster(abundance, proj_ext, method = "bilinear"), 
                  mask = temp_mask, maskvalue = 0)
abundance <- cover(abundance, templ_rasprj)
q99 <- quantile(values(abundance[[3]]), 0.95, na.rm = TRUE)
abundance[abundance > q99] <- q99
new_res <- res(abundance)[1]/1000
values(abundance) <- floor(values(abundance)*(new_res*new_res))
plot(abundance, zlim = c(0, floor(q99*(new_res*new_res))),
     col = hcl.colors(100, "Fall"),
     legend = FALSE,
     addfun = function() lines(as_Spatial(st_transform(conts, st_crs(proj_ext)))))

# paleopop requires a relative harvest rate - stretch to [0, 1]
abundance <- stretch(abundance, minv = 0, maxv = 1, smin = 0, smax = floor(q99*(new_res*new_res)))
abundance
plot(abundance, zlim = c(0, 1), 
     col = hcl.colors(100, "Fall"), legend = FALSE,
     addfun = function() lines(as_Spatial(st_transform(conts, st_crs(proj_ext)))))

par(mfrow = c(1, 3))
hist(abundance[[1]])
hist(abundance[[2]])
hist(abundance[[3]])
par(mfrow = c(1, 1))

# convert abundance to matrix
abundance_matrix <- round(as.matrix(abundance * 10000), 0)
mode(abundance_matrix) <- "integer"
attr(abundance_matrix, "timesteps") <- NA
attr(abundance_matrix, "res") <- res(abundance)[1]
attr(abundance_matrix, "extent") <- "-2423364, 2436636, -2158270, 1981730"
attr(abundance_matrix, "proj4") <- "+proj=laea +lat_0=-27 +lon_0=133 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
attr(abundance_matrix, "proj_wkt") <- wkt(abundance)
attr(abundance_matrix, "scale.fac") <- 10000
attr(abundance_matrix, "X") <- coordinates(abundance)[, 1]
attr(abundance_matrix, "Y") <- coordinates(abundance)[, 2]
str(abundance_matrix)

# dims, dimnames, timesteps, and scale.fac are different so dont test
sapply(names(attributes(abundance_matrix))[-c(1:3, 8)], function(x) {
  all.equal(attr(abundance_matrix, x), attr(region, x))
})

saveRDS(abundance_matrix, file = "data/fox_invasion_abundance/fox_abundance_estimate.RDS", compress = TRUE)
