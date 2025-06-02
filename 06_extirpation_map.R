library(raster)
library(poems)
library(paleopop)
library(abc)
library(data.table)
library(sf)
library(rnaturalearth)

burn_in_steps <- 50 # 25 generations
timesteps <- 271 + burn_in_steps ## to start at 30k BP
timeseq <- seq(to = 2020, by = 1, length.out = timesteps)

# read in paleopop region
region <- readRDS("third_try/sim_inputs/numbat_paleopopregion.RDS")
plot(region$region_raster)

r <- region$region_raster
r[] <- +(!is.na(r[]))
r[r == 0] <- NA
poly_region <- rasterToPolygons(r, dissolve = TRUE)

# Continents
conts <- ne_countries(scale = 50, returnclass = "sf")
conts <- conts[conts$name == "Australia", ]
conts <- st_transform(conts, crs = st_crs(region$region_raster))
conts <- st_crop(conts, extent(region$region_raster))

# location of simulations
out_dir <- "third_try/out_sims_final"

# read in ABC object
first_pass <- readRDS("third_try/final_selection/abc_first_pass_filtered_try2.RDS")
summary(first_pass)

# read in simulation params
sim_par <- fread("third_try/ST_abc_v3/lhs2_corrected.csv")
sim_par[, UID := 1:.N]
setcolorder(sim_par, c(13, 1:12))
sim_par[, Region := first_pass$region]
sim_par <- sim_par[Region == TRUE, ]
sim_par

# read in extirpation dates
ext_dates <- gsub("UniqueID_|_results.RData", "",
                  list.files(out_dir, ".RData"))
stopifnot(sum(ext_dates %in% sim_par[["UniqueID"]]) == nrow(sim_par))

ext_files <- file.path(out_dir, paste0("UniqueID_", sim_par[["UniqueID"]], "_results.RData"))
ext_files

ext_list <- lapply(ext_files, function(i,...) {
  r <- region$raster_from_values(readRDS(i)$extirpation)+timeseq[1]
  # NA == never occupied OR never extinct. Set to 2022
  r[is.na(r)] <- 2022
  # mask back to region
  r <- mask(r, region$region_raster)
  # extinction during burn in == NA
  r[r < (timeseq[burn_in_steps] + 1)] <- NA
  return(r)
})
ext_stack <- stack(ext_list, quick = TRUE)
ext_stack[[seq(1, 50, l = 8)]]
plot(ext_stack[[seq(1, 50, l = 8)]], zlim = c(1750, 2022), legend = FALSE,
     addfun = function() lines(as_Spatial(conts)))

# weighted median extinction time
ext_time <- calc(x = ext_stack, function(x,...) {
  i <- x
  if (sum(is.na(i)) == length(i)) {
    return(NA)
  } else {
    i <- Hmisc::wtd.quantile(i, w = (1/first_pass$dist[first_pass$region == TRUE]),
                             probs = 0.5, normwt = TRUE, na.rm = TRUE)
    return(i)
  }
})
ext_time

#writeRaster(ext_time,
#            filename = "third_try/weighted_extirpation2_10000.grd",
#            datatype = "INT4U")

brks <- c(1750, seq(1760, 2022, by = 20))
brks
labs <- as.character(brks)
labs <- labs[-1]
labs[1] <- "<1760"
labs[length(labs)] <- ">2020"

p1 <- function(){
  par(mar = c(0,0,0,4))
  plot(as_Spatial(conts), axes = FALSE)
  plot(poly_region, col = "grey85", lty = 3, add = TRUE)
  plot(ext_time, zlim = c(1750,2022),
       breaks = c(brks, 2022),
       box = FALSE,
       col = hcl.colors(15, "Fall", rev = TRUE),
       axes = FALSE, legend = FALSE,
       addfun = function() {
         lines(as_Spatial(conts))
       },
       add = TRUE)
  legend(x = 2.025e6, y = 1.35e6, cex = 0.75,
         bty = "n", text.font = 2,
         legend = "Extirpation\nYear")
  legend(
    x = 2.15e6, y = 1e6, bty = "n",
    legend = labs,
    ncol = 1, cex = 0.75,
    fill = hcl.colors(14, "Fall", rev = TRUE))
}
p1()
{png(filename = "third_try/weighted_extirpation2_10000.png",
     width = 8, height = 5, units = "in", 
     bg = "white", res = 320, type = "cairo",
     antialias = "subpixel")
  p1()
  dev.off()
  if (require(magick)) {
    g1 <- image_read("third_try/weighted_extirpation2_10000.png")
    g1 <- image_trim(g1)
    image_write(image = g1, path = "third_try/weighted_extirpation2_10000.png")
  }}











library(raster)
library(poems)
library(paleopop)
library(abc)
library(data.table)
library(sf)
library(rnaturalearth)

burn_in_steps <- 30 # 25 generations
timesteps <- 271 + burn_in_steps ## to start at 30k BP
timeseq <- seq(to = 2020, by = 1, length.out = timesteps)

# read in paleopop region
region <- readRDS("third_try/sim_inputs/numbat_paleopopregion.RDS")
plot(region$region_raster)

r <- region$region_raster
r[] <- +(!is.na(r[]))
r[r == 0] <- NA
poly_region <- rasterToPolygons(r, dissolve = TRUE)

# Continents
conts <- ne_countries(scale = 50, returnclass = "sf")
conts <- conts[conts$name == "Australia", ]
conts <- st_transform(conts, crs = st_crs(region$region_raster))
conts <- st_crop(conts, extent(region$region_raster))

# location of simulations
out_dir <- "third_try/out_sims2"

# read in ABC object
first_pass <- readRDS("third_try/abc/abc_first_pass_10000_3run.RDS")
summary(first_pass)

# read in simulation params
sim_par <- fread("third_try/lhs3_10000.csv")
sim_par[, UID := 1:.N]
setcolorder(sim_par, c(13, 1:12))
sim_par[, Region := first_pass$region]
sim_par <- sim_par[Region == TRUE, ]
sim_par

# read in extirpation dates
ext_dates <- gsub("UniqueID_|_results.RData", "",
                  list.files(out_dir, ".RData"))
stopifnot(sum(ext_dates %in% sim_par[["UniqueID"]]) == nrow(sim_par))

ext_files <- file.path(out_dir, paste0("UniqueID_", sim_par[["UniqueID"]], "_results.RData"))
ext_files

ext_list <- lapply(ext_files, function(i,...) {
  r <- region$raster_from_values(readRDS(i)$extirpation)+timeseq[1]
  # NA == never occupied OR never extinct. Set to 2022
  r[is.na(r)] <- 2022
  # mask back to region
  r <- mask(r, region$region_raster)
  # extinction during burn in == NA
  r[r < (timeseq[burn_in_steps] + 1)] <- NA
  return(r)
})
ext_stack <- stack(ext_list, quick = TRUE)
ext_stack[[seq(1, 50, l = 8)]]
plot(ext_stack[[seq(1, 50, l = 8)]], zlim = c(1750, 2022), legend = FALSE,
     addfun = function() lines(as_Spatial(conts)))

# weighted median extinction time
ext_time <- calc(x = ext_stack, function(x,...) {
  i <- x
  if (sum(is.na(i)) == length(i)) {
    return(NA)
  } else {
    i <- Hmisc::wtd.quantile(i, w = (1/first_pass$dist[first_pass$region == TRUE]),
                             probs = 0.5, normwt = TRUE, na.rm = TRUE)
    return(i)
  }
})
ext_time

writeRaster(ext_time,
            filename = "third_try/weighted_extirpation3_10000.grd",
            datatype = "INT4U")

brks <- c(1750, seq(1760, 2022, by = 20))
brks
labs <- as.character(brks)
labs <- labs[-1]
labs[1] <- "<1760"
labs[length(labs)] <- ">2020"

p1 <- function(){
  par(mar = c(0,0,0,4))
  plot(as_Spatial(conts), axes = FALSE)
  plot(poly_region, col = "grey85", lty = 3, add = TRUE)
  plot(ext_time, zlim = c(1750,2022),
       breaks = c(brks, 2022),
       box = FALSE,
       col = hcl.colors(15, "Fall", rev = TRUE),
       axes = FALSE, legend = FALSE,
       addfun = function() {
         lines(as_Spatial(conts))
       },
       add = TRUE)
  legend(x = 2.025e6, y = 1.35e6, cex = 0.75,
         bty = "n", text.font = 2,
         legend = "Extirpation\nYear")
  legend(
    x = 2.15e6, y = 1e6, bty = "n",
    legend = labs,
    ncol = 1, cex = 0.75,
    fill = hcl.colors(14, "Fall", rev = TRUE))
}
p1()
{png(filename = "third_try/weighted_extirpation3_10000.png",
     width = 8, height = 5, units = "in", 
     bg = "white", res = 320, type = "cairo",
     antialias = "subpixel")
  p1()
  dev.off()
  if (require(magick)) {
    g1 <- image_read("third_try/weighted_extirpation3_10000.png")
    g1 <- image_trim(g1)
    image_write(image = g1, path = "third_try/weighted_extirpation2_10000.png")
  }}

