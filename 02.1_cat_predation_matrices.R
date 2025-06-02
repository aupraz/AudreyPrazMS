library(raster)
library(sf)
library(parallel)
library(pbapply)
library(rnaturalearth)

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

##### Cat arrival ####
cat_arrival <- projectRaster(raster("data/cat_predation_layers/cat_invasion_20km.grd"),
                             proj_ext, method = "ngb")
cat_arrival; unique(cat_arrival)  # arrival times are years BP (2020 C.E.)
values(cat_arrival) <- 2020 - values(cat_arrival)  # Correct arrival times
unique(cat_arrival)
plot(cat_arrival)

# explode into unique layers based on values
cat_arrival <- layerize(cat_arrival)
plot(cat_arrival, 
     addfun = function() {
       lines(as_Spatial(st_transform(conts, crs = st_crs(proj_ext))))
     },
     legend = FALSE)

# Adjust layers to ensure cumulative inclusion
for (i in 2:nlayers(cat_arrival)) {
  s <- stack(cat_arrival[[i]], cat_arrival[[i - 1]])
  s <- calc(s, max)
  values(cat_arrival[[i]]) <- values(s)
}
options(scipen=999)
plot(cat_arrival, 
     addfun = function() {
       lines(as_Spatial(st_transform(conts, crs = st_crs(proj_ext))))
     },
     legend = FALSE)

# Set correct Z values for 1820-2020
cat_arrival <- setZ(cat_arrival, z = seq(1820, 1890, by = 10))

# Create the `cat_arrival_int` stack with adjusted years
cat_arrival_int <- stack(replicate(n = length(seq(1820, 2020, by = 1)), cat_arrival[[1]]), quick = TRUE)
names(cat_arrival_int) <- paste0("X_", 1820:2020)
values(cat_arrival_int) <- NA
cat_arrival_int <- setZ(cat_arrival_int, 1820:2020, "Year")

# Interpolate the arrival times for 1820-2020
for (time in getZ(cat_arrival)) {
  j <- which(getZ(cat_arrival_int) == time)
  values(cat_arrival_int[[j]]) <- values(cat_arrival[[which(getZ(cat_arrival) == time)]])
}

# Mask cat_arrival_int with Australia's boundaries
conts_raster <- rasterize(as_Spatial(st_transform(conts, crs = st_crs(proj_ext))),
                          cat_arrival_int[[1]], field = 1, background = NA)
cat_arrival_int <- mask(cat_arrival_int, conts_raster)

# Fill gaps in the stack using linear interpolation
cat_arrival_int <- approxNA(cat_arrival_int, method = "linear", yleft = 0, yright = 1)
cat_arrival_int[cat_arrival_int < 0.5] <- 0
cat_arrival_int[cat_arrival_int > 0] <- 1

# Verify the stack after masking
plot(cat_arrival_int[[50:62]], zlim = c(0, 1), nr = 4,
     addfun = function() {
       lines(as_Spatial(st_transform(conts, crs = st_crs(proj_ext))))
     },
     legend = FALSE)


# Copy the first layer and replicate it
first_layer <- cat_arrival_int[[1]]
values(first_layer)[values(first_layer) == 1] <- 0  # Turn 1s into 0s
prepend_stack <- stack(replicate(70, first_layer))

# Merge with existing stack
cat_arrival_int <- stack(prepend_stack, cat_arrival_int)


# convert cats to temporal arrival mask
names(cat_arrival_int) <- paste0("X_",seq(1750, by = 1, length = nlayers(cat_arrival_int)))
cat_arrival_int <- setZ(cat_arrival_int, seq(1750, by = 1, length = nlayers(cat_arrival_int)), "Years")
plot(cat_arrival_int[[seq(1, 271, length = 9)]], zlim = c(0, 1))

# Convert to matrix and save
cat_arrival_intMat <- as.matrix(cat_arrival_int)
#object.size(hv_rast_maskedMat)
mode(cat_arrival_intMat) <- "integer"
#object.size(hv_rast_maskedMat)
attr(cat_arrival_intMat, "timesteps") <- seq(1750, by = 1, length = ncol(cat_arrival_intMat))
attr(cat_arrival_intMat, "res") <- res(cat_arrival_int)[1]
attr(cat_arrival_intMat, "extent") <- "-2423364, 2436636, -2158270, 1981730"
attr(cat_arrival_intMat, "proj4") <- "+proj=laea +lat_0=-27 +lon_0=133 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
attr(cat_arrival_intMat, "proj_wkt") <- wkt(cat_arrival_int)
attr(cat_arrival_intMat, "scale.fac") <- NA
attr(cat_arrival_intMat, "X") <- coordinates(cat_arrival_int)[,1]
attr(cat_arrival_intMat, "Y") <- coordinates(cat_arrival_int)[,2]
str(cat_arrival_intMat)

sapply(names(attributes(cat_arrival_intMat)), function(x) {
  all.equal(attr(cat_arrival_intMat, x), attr(region, x))
})


saveRDS(cat_arrival_intMat, file = "data/cat_predation_layers/cat_arrival_temporal_mask.RDS", compress = TRUE)

#### Upper estimate ####
# original res in km
orig_res <- res(raster("data/cat_predation_layers/upper bound.tif"))[1]/1000

upper <- raster("data/cat_predation_layers/upper bound.tif")/orig_res
lower <- raster("data/cat_predation_layers/lower bound.tif")/orig_res
mean <- raster("data/cat_predation_layers/estimate.tif")/orig_res
abundance <- stack(lower, mean, upper, quick = TRUE)
temp_mask <- templ_rasprj
temp_mask[is.na(temp_mask) | temp_mask == 0] <- 0
abundance <- aggregate(abundance, fact = res(cat_arrival_int)/res(abundance), 
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

# paleopop requires a relative abundance rate - stretch to [0, 1]
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

saveRDS(abundance_matrix, file = "data/cat_predation_layers/cat_abundance_estimate.RDS", compress = TRUE)











