library(raster)
library(sf)
library(parallel)
library(pbapply)
library(rnaturalearth)

# Analysis region
region <- raster("paleopop_numbats/sim_inputs/numbat_region_20km.grd")
region

# cont outlines
conts <- ne_countries(scale = 50, returnclass = "sf")
conts <- conts[conts$name == "Australia", ]
conts <- st_crop(conts, extent(region))

plot(region, addfun = function() lines(as_Spatial(conts)))

# input projections
hv_unproj <- list.files("data/niche_cuts/", pattern = "\\.grd$",
                        recursive = TRUE, full.names = TRUE)
head(hv_unproj)

plot(raster(hv_unproj[1]), addfun = function() lines(as_Spatial(conts)),
     zlim = c(0, 1), col = hcl.colors(100, "Fall", rev = TRUE))

# temp file and directory for unzipping
temp <- tempfile()
tempd <- tempdir()

bras <- region
bras <- crop(bras, region)
bras_ext <- projectExtent(bras, crs = 'PROJCS["Custom_Lambert_Azimuthal",
 GEOGCS["GCS_WGS_1984",
  DATUM["D_WGS_1984",
   SPHEROID["WGS_1984",6378137.0,298.257223563]],
  PRIMEM["Greenwich",0.0],
  UNIT["Degree",0.0174532925199433]],
 PROJECTION["Lambert_Azimuthal_Equal_Area"],
 PARAMETER["False_Easting",0.0],
 PARAMETER["False_Northing",0.0],
 PARAMETER["Central_Meridian",133],
 PARAMETER["Latitude_Of_Origin",-27],
 UNIT["Meter",1.0]]')
res(bras_ext) <- c(20000,20000)
bras_ext

plot(projectRaster(region, bras_ext, method = "ngb"), addfun = function(x)
  lines(as_Spatial(st_transform(conts, st_crs(bras_ext)))), legend = FALSE)

plot(projectRaster(brick(hv_unproj[1])[[seq(1, 271, length = 4)]], bras_ext, method = "bilinear"), 
     addfun = function(x) {
  lines(as_Spatial(st_transform(conts, st_crs(bras_ext))))},
  legend = FALSE, zlim = c(0, 1), col = hcl.colors(100, "Fall", rev = TRUE))

comm_ext_prj <- projectRaster(region, bras_ext, method = "ngb")

# common extent as matrix
comm_ext_prj_temp <- stack(replicate(n = 271, comm_ext_prj), quick = TRUE)
names(comm_ext_prj_temp) <- paste0("X_",seq(1750, by = 1, length = nlayers(comm_ext_prj_temp)))
# Convert to matrix and save
comm_ext_Mat <- as.matrix(comm_ext_prj_temp)
#object.size(hv_rast_maskedMat)
mode(comm_ext_Mat) <- "integer"
#object.size(hv_rast_maskedMat)
attr(comm_ext_Mat, "timesteps") <- seq(1750, by = 1, length = ncol(comm_ext_Mat))
attr(comm_ext_Mat, "res") <- res(comm_ext_prj_temp)[1]
attr(comm_ext_Mat, "extent") <- "-2423364, 2436636, -2158270, 1981730"
attr(comm_ext_Mat, "proj4") <- "+proj=laea +lat_0=-27 +lon_0=133 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
attr(comm_ext_Mat, "proj_wkt") <- wkt(comm_ext_prj_temp)
attr(comm_ext_Mat, "scale.fac") <- NA
attr(comm_ext_Mat, "X") <- coordinates(comm_ext_prj_temp)[,1]
attr(comm_ext_Mat, "Y") <- coordinates(comm_ext_prj_temp)[,2]
str(comm_ext_Mat)
saveRDS(comm_ext_Mat, file = "first_try/numbat_temporal_region.RDS", compress = TRUE)

# Iterate through each hypervolume projection,
# mask, convert to matrix, and save as RDS

cls <- makeCluster(25L)
clusterExport(cls, c("hv_unproj", "comm_ext_prj", "bras_ext", "conts"))

pblapply(hv_unproj, function(hv) {
  require(raster); require(sf)
  #hv_name <- gsub("\\.grd$", "", paste(strsplit(hv, "_")[[1]][4:5], collapse = "_"), )
  #outfile <- file.path(
  #  "data/niche_cuts/proj_matrices/projected",paste0(hv_name,".RDS"))
  hv_name <- gsub("\\.grd$", "", basename(hv))  # Extract filename without directory and remove ".grd"
  hv_name <- gsub("cut_", "", hv_name)          # Remove "cut_" prefix
  outfile <- file.path("data/niche_cuts/proj_matrices/projected", paste0(hv_name, ".RDS"))
  
  if (file.exists(outfile)) {
    return(NULL)
  }
  hv_rast <- brick(hv)
  hv_rast <- readAll(hv_rast)
  # project
  hv_rast <- readAll(projectRaster(hv_rast, bras_ext, method = "bilinear"))
  # Multiply by spatio-temporal mask
  stopifnot(compareRaster(hv_rast, comm_ext_prj))
  hv_rast <- hv_rast * comm_ext_prj
  hv_rast[hv_rast < 0] <- 0
  hv_rast[hv_rast > 1] <- 1
  # plot(hv_rast[[seq(1, nlayers(hv_rast), length = 4)]],
  #      addfun = function(x)
  #        lines(as_Spatial(st_transform(conts, st_crs(bras_ext)))))
  values(hv_rast) <- round(values(hv_rast), 4)
  names(hv_rast) <- paste0("X_",seq(1750, by = 1, length = nlayers(hv_rast)))
  # Convert to matrix and save
  hv_rast_maskedMat <- round(as.matrix(hv_rast * 1000), 0) # multiply by 1000 for scaling
  #object.size(hv_rast_maskedMat)
  mode(hv_rast_maskedMat) <- "integer"
  #object.size(hv_rast_maskedMat)
  attr(hv_rast_maskedMat, "timesteps") <- seq(1750, by = 1, length = ncol(hv_rast_maskedMat))
  attr(hv_rast_maskedMat, "res") <- res(hv_rast)[1]
  attr(hv_rast_maskedMat, "extent") <- "-2423364, 2436636, -2158270, 1981730"
  attr(hv_rast_maskedMat, "proj4") <- "+proj=laea +lat_0=-27 +lon_0=133 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  attr(hv_rast_maskedMat, "proj_wkt") <- wkt(hv_rast)
  attr(hv_rast_maskedMat, "scale.fac") <- 1000
  attr(hv_rast_maskedMat, "X") <- coordinates(hv_rast)[,1]
  attr(hv_rast_maskedMat, "Y") <- coordinates(hv_rast)[,2]
  str(hv_rast_maskedMat)
  saveRDS(hv_rast_maskedMat, file = outfile, compress = TRUE)
  # delete temp files
  removeTmpFiles(h = (0.0083)/3) # 10 seconds
  return("DONE")
}, cl = cls)
stopCluster(cls)

# Check all correct
f_list <- list.files(path =  "data/niche_cuts/proj_matrices/projected",
                     pattern = "\\.RDS$", full.names = TRUE)
head(f_list); tail(f_list)
cls <- makeCluster(10L)
clusterExport(cls, c("f_list"))
dim_check <- pblapply(f_list,
                     function(i) {
                       message(i)
                       f <- readRDS(i)
                       return(c(dim(f)[1] == 50301, dim(f)[2] == 271))
                     }, cl = cls)
stopCluster(cls)
dim1 <- sapply(dim_check, "[", 1:2)[1,]
dim2 <- sapply(dim_check, "[", 1:2)[2,]
sum(dim1) == sum(dim2) # TRUE

