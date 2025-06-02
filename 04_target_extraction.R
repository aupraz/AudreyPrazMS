#### Libraries ####
library(raster)
library(poems)
library(paleopop)
library(sf)
library(rnaturalearth)
library(data.table)
library(pbapply)
library(doSNOW)

library(pbarETA, warn = FALSE)

#### Workflow options ####

# Generate a projected template regional raster
temp_mask <- readRDS("data/numbat_temporal_region.RDS")
templ_ras <- raster(res = 0.20,
                    ext = extent(111.05,155.05,-45.025,-9.025),
                    crs = crs(raster()))
values(templ_ras) <- 1
proj_ext <- projectExtent(templ_ras, crs = attr(temp_mask, "proj_wkt"))
res(proj_ext) <- attr(temp_mask, "res")
all.equal(coordinates(proj_ext)[,1],  attr(temp_mask, "X")) ## TRUE
all.equal(coordinates(proj_ext)[,2],  attr(temp_mask, "Y")) ## TRUE
proj_ext

# simulation inputs
input_dir <- "third_try/sim_inputs/"

# location of HS projections
data_dir <- "data/niche_cuts/proj_matrices/projected"

# location of validation targets
target_dir <- "third_try/validation_targets"

# location of simulations
sim_dir <- "third_try/out_sims"

# time seq not including burnin
parallel_cores <- 20L # cluster: 20 ?
nsims <- length(list.files(sim_dir, pattern = ".RData"))

str(readRDS("data/numbat_temporal_region.RDS"), max.level = 2)

burn_in_steps <- 50 # 25 generations
timesteps <- 271 + burn_in_steps ## to start at 30k BP
plot_seq <- floor(seq(burn_in_steps+1, timesteps, length.out = 6))

timeseq <- seq(to = 2020, by = 1, length.out = (timesteps - burn_in_steps))



weightedCentre <- function(x, y, z) {
  require(matrixStats)
  if (anyNA(c(x, y, z))) {
    stop("There are missing values present in x, y, or z")
  }
  if (length(z) != length(x)) {
    stop("Number of weights supplied not equal to number of coordinates")
  }
  x_y <- cbind(x, y)
  w <- z
  mean <- matrixStats::colWeightedMeans(x_y, w = w, na.rm = FALSE)
  wght_cent <- matrix(mean, ncol = 2)
  rownames(wght_cent) <- c("mean")
  return(wght_cent)
}

#### Step 1: Create targets ####

# paleopop region
region <- readRDS(file.path(input_dir, "numbat_paleopopregion.RDS"))
region$region_raster

# Aus shape
conts <- ne_countries(scale = 50, returnclass = "sf")
conts <- conts[conts$name == "Australia", ]
conts <- st_transform(conts, crs = st_crs(region$region_raster))
conts <- st_crop(conts, extent(region$region_raster))

plot(region$region_raster,
     col = hcl.colors(ncell(region$region_raster), "Lajolla"),
     legend = FALSE,
     addfun = function() lines(as_Spatial(conts)))

# load target files

## Target 1: known fossil locations
target_locations <- fread(file.path(target_dir, "numbat_targets_new.csv"), select = 1:6)
colnames(target_locations) <- c("ID", "Spp", "Year", "Latitude", "Longitude", "Type")
target_locations <- st_as_sf(target_locations,
                             coords = c("Longitude", "Latitude"),
                             crs = "EPSG:4326")
target_locations <- st_transform(target_locations, crs = st_crs(region$region_raster))

# cell indices for projected locations
target_locations$region_index <- cellFromXY(region$region_raster, st_coordinates(target_locations))
target_locations
plot(region$region_raster, col = "#a3be8c",
     legend = FALSE,
     addfun = function() {
       points(as_Spatial(target_locations[1]), pch = 19, cex = 0.25, col = "red")
       lines(as_Spatial(conts))
     })

## Target 2: regional extirpation times
final_locations <- fread(file.path(target_dir, "final_locations.csv"), select = 1:6)
colnames(final_locations) <- c("ID", "Spp", "Year", "Latitude", "Longitude", "Type")
final_locations <- st_as_sf(final_locations,
                            coords = c("Longitude", "Latitude"),
                            crs = "EPSG:4326")
final_locations <- st_transform(final_locations, crs = st_crs(region$region_raster))

# cell indices for projected locations
final_locations$region_index <- cellFromXY(region$region_raster, st_coordinates(final_locations))
final_locations

extinction_bioregions <- raster(file.path(target_dir, "extinction_layer2.grd"))
crs(extinction_bioregions) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
extinction_bioregions <- projectRaster(extinction_bioregions, proj_ext, method = "ngb")
extinction_bioregions[values(extinction_bioregions) == 0] <- NA
#extinction_bioregions[values(extinction_bioregions) == "2020"] <- 2009 # change 2020 to 2009 as per signor-lipps date
plot(extinction_bioregions)

extinction_bioregions[] <- as.factor(extinction_bioregions[])
levels(extinction_bioregions) <- levels(extinction_bioregions)[[1]]
extinction_bioregions

{plot(extinction_bioregions,
      col = hcl.colors(6, "TealRose"), legend = FALSE,
      addfun = function() {
        points(as_Spatial(target_locations[1]), pch = 19, cex = 0.25, col = "red")
        lines(as_Spatial(conts))
      },
      axes = FALSE, box = FALSE)
  legend(
    x = 2.15e6, y = 0.75, bty = "n",
    legend = levels(extinction_bioregions)[[1]][,2],
    fill = hcl.colors(6, "TealRose"))
}

# convert to polygons
subregions <- st_as_sf(rasterToPolygons(extinction_bioregions, dissolve = TRUE))
subregions$Year <- levels(extinction_bioregions)[[1]][, 2]

#subregions <- subregions[subregions$extinction_layer >1,]

plot(subregions)

regional_pops <- lapply(unique(subregions$Year), function(reg) {
  # if remnant pop, use surrounding cells
  if (reg == "2022") {
    reg_buff <- st_buffer(subregions[subregions$Year == reg, ],
                          dist = 100*1000, # 100 km
                          endCapStyle = "SQUARE", joinStyle = "MITRE")
    ras_region <- raster::mask(x = region$region_raster,
                               mask = reg_buff)
  } else {
    ras_region <- raster::mask(x = region$region_raster,
                               mask = subregions[subregions$Year == reg, ])
  }
  values(ras_region)[which(!is.na(values(ras_region)))]
})
names(regional_pops) <- unique(subregions$Year)
sum_cells <- Reduce("+", lapply(regional_pops, function(x) length(unique(x))))
sum_cells
str(regional_pops)

locs.2022 <- sapply(1:2, function(i,...) {
  cell <- cellFromXY(region$region_raster,
                     xy = st_coordinates(final_locations[i, ]))
  cells <- raster::adjacent(region$region_raster, cell, include = TRUE, directions = 8)[, 2]
  cells_ID <- region$region_raster[cells]
  return(cells_ID)
})
locs.2022 <- as.vector(locs.2022)

regional_pops <- regional_pops
#regional_pops$`2022` <- locs.2022

# quick check that cells are in correct spots
{pops_1886 <- region$coordinates[regional_pops$`1886`, ]
  pops_1938 <- region$coordinates[regional_pops$`1938`, ]
  pops_2003 <- region$coordinates[regional_pops$`2003`, ]
  pops_2009 <- region$coordinates[regional_pops$`2008`, ]
  pops_2022 <- region$coordinates[regional_pops$`2022`, ]
  plot(region$region_raster, col = "grey80", legend = FALSE, box = FALSE,
       axes = FALSE)
  points(pops_1886, pch = 16, col = "#1b9e77")
  points(pops_1938, pch = 16, col = "#d95f02")
  points(pops_2003, pch = 16, col = "#7570b3")
  points(pops_2009, pch = 16, col = "red")
  points(pops_2022, pch = 16, col = "blue")
  
  legend("right", c("1886", "1938", "2003", "2008", "2022"), bty = "n",
         pch = 16, col = c("#1b9e77", "#d95f02", "#7570b3", "red", "blue", "grey80"))}

#### Step 2: Extract targets ####


sample_data<- fread("third_try/lhs2_10000.csv")
out_sims <- file.path(sim_dir, paste0("UniqueID_", sample_data$UniqueID, "_results.RData"))
idx <- which(file.exists(out_sims))
out_sims <- out_sims[idx]
sample_data <- sample_data[idx, ]

ext_times <- data.frame(
  Region = names(regional_pops),
  extTime = as.integer(names(regional_pops))
)
ext_times

length(out_sims) == nrow(sample_data)


# loop through each simulation and calculate the targets
cl <- makeCluster(20)
registerDoSNOW(cl)
iterations <- nrow(sample_data)
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) {setTxtProgressBar(pb, n)}
opts <- list(progress = progress)


getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


ext_metrics <- foreach(sample_row = seq_len(nrow(sample_data)),
                       .packages = c("raster", "poems", "paleopop",
                                     "data.table", "sf", "matrixStats", "dplyr"),
                       .inorder = TRUE,
                       .errorhandling = "pass",
                       .verbose = FALSE) %dopar% {
                         #.options.snow = opts) %dopar% {
                         sim <- readRDS(out_sims[sample_row])$abundance[, -c(1:burn_in_steps)]
                         
                         # target_1: persistence of the whole population
                         abund_seq <- colSums(sim,na.rm = TRUE)
                         persistence <- length(abund_seq[abund_seq > 0])
                         
                         # number of cells surrounding each validation site
                         ## using 9 cells with maximum year of each focal cell
                         match_9 <- sapply(1:nrow(target_locations), function(i,...) {
                           j <- target_locations[i, ]$Year
                           yr <- which(timeseq == j)
                           r <- region$raster_from_values(sim[, yr])
                           r[] <- +(r[]>0)
                           cell <- cellFromXY(r,
                                              xy = st_coordinates(target_locations[i, ]))
                           cells <- raster::adjacent(r, cell, include = TRUE, directions = 8)[, 2]
                           occ_cells <- sum(as.vector(raster::extract(x = r, y = cells)), na.rm = TRUE)
                           occ_cells <- +(occ_cells > 0)
                           occ_cell_tab <- cbind(cell, j, occ_cells)
                           #return(occ_cells)
                           return(occ_cell_tab)
                         })
                         
                         dat <- as.data.frame(t(match_9))
                         colnames(dat) <- c("cell", "year", "occ")
                         
                         for(u in 1:length(unique(dat$cell))){
                           dat1 <- dat[dat$cell == unique(dat$cell)[u],]
                           dat1 <- dat1[dat1$year == max(dat1$year),]
                           
                           if(u == 1){match_9 <- dat1}else{match_9 <- rbind(match_9, dat1)}
                         
                           }
                        
                         match_9 <- match_9 %>% distinct(cell, .keep_all = TRUE)
                         
                         # How many validation sites were simulated correctly in time/space across the
                         # entire study region?
                         sites_9 <- sum(match_9$occ)
                         sites_9perc <- round((sum(match_9$occ)/nrow(match_9))*100)/100
                         ## using 1 cells
                         match_1 <- sapply(1:nrow(target_locations), function(i,...) {
                           j <- target_locations[i, ]$Year
                           yr <- which(timeseq == j)
                           r <- region$raster_from_values(sim[, yr])
                           r[] <- +(r[]>0)
                           cell <- cellFromXY(r,
                                              xy = st_coordinates(target_locations[i, ]))
                           occ_cells <- sum(as.vector(raster::extract(x = r, y = cell)), na.rm = TRUE)
                           occ_cells <- +(occ_cells > 0)
                           occ_cell_tab <- cbind(cell, j, occ_cells)
                           #return(occ_cells)
                           return(occ_cell_tab)
                         })
                         
                         dat <- as.data.frame(t(match_1))
                         colnames(dat) <- c("cell", "year", "occ")
                         
                         for(u in 1:length(unique(dat$cell))){
                           dat1 <- dat[dat$cell == unique(dat$cell)[u],]
                           dat1 <- dat1[dat1$year == max(dat1$year),]
                           
                           if(u == 1){match_1 <- dat1}else{match_1 <- rbind(match_1, dat1)}
                           
                         }
                         
                         match_1 <- match_1 %>% distinct(cell, .keep_all = TRUE)
                         # How many validation sites were simulated correctly in time/space across the
                         # entire study region?
                         sites_1 <- sum(match_1$occ)
                         sites_1perc <- round((sum(match_1$occ)/nrow(match_1))*100)/100
                         
                         
                         # final location
                         persist_9 <- sapply(1:2, function(i,...) {
                           #j <- target_locations[i, ]$Year
                           #yr <- which(timeseq == j)
                           r <- region$raster_from_values(sim[, 271]) 
                           r[] <- +(r[]>0)
                           
                           # cell <- cellFromXY(r,
                           #                    xy = c(-32.783143749939406, 116.97204892883563))
                           
                           
                           cell <- cellFromXY(r,
                                              xy = st_coordinates(final_locations[i, ]))
                           cells <- raster::adjacent(r, cell, include = TRUE, directions = 8)[, 2]
                           #final_abund <- sum(as.vector(raster::extract(x = r, y = cells)), na.rm = TRUE)
                           #if(final_abund>0){occ_cells <- 1}else{occ_cells <- 0}
                           occ_cells <- sum(as.vector(raster::extract(x = r, y = cells)), na.rm = TRUE)
                           occ_cells <- +(occ_cells > 0)
                           return(occ_cells)
                         })
                         
                         persist_9 <- mean(persist_9)
                         
                         # final abundance
                         for(r in c(261:271)){
                           abund_rast <- region$raster_from_values(sim[, r])
                           est <- base::sum(values(abund_rast), na.rm = TRUE)
                           
                           if(r == 261){est.dat <- est}else{est.dat <- rbind(est.dat,est)}
                         }
                         
                         #final_abund <- round(base::mean(est.dat),0)
                         final_abund <- round(est.dat[11],0)
                         
                         #distance from final locations
                         abundance_array<-sim[,persistence] #here, we take the layer of the last known populations, regardless of whether they persist to today or not
                         
                         abundance_array[abundance_array == 0] <-NA
                         abundance_array[abundance_array > 0] <-2
                         
                         abundance.rast <- region$region_raster*0
                         abundance.rast[region$region_indices] <- abundance_array
                         #extinction_location <- raster::xyFromCell(abundance.rast, which(as.logical(raster::as.matrix(abundance.rast))))
                         for(i in unique(values(abundance.rast))){
                          p1 <- rasterToPoints(abundance.rast, fun = function(xx) xx==i, spatial = T)
                         }
                         
                         target_raster <- region$region_raster*0
                         target_raster[region$region_indices][locs.2022] <- 2
                         #target_location <- raster::xyFromCell(target_raster, which(as.logical(raster::as.matrix(target_raster))))
                         
                         for(i in unique(values(target_raster))){
                           p2 <- rasterToPoints(target_raster, fun = function(xx) xx==i, spatial = T)
                         }
                         
                         # if(sum(abundance_array, na.rm = TRUE) == 0){
                         #   ext_dist <- NA
                         # }else{
                         #   ext_dist <- spDists(p1,p2)
                         #   ext_dist <- round(min(as.vector(ext_dist)),0)
                         # }
                         
                         ext_dist <- spDists(p1,p2)
                         ext_dist <- round(mean(as.vector(ext_dist)),0) # mean distance
                         #ext_dist <- getmode(as.vector(round(ext_dist,0))) # modal distance
                         
                         regional_metrics <- lapply(names(regional_pops), function(reg,...) {
                           message(reg)
                           extTime <- (ext_times[ext_times$Region == reg, ][["extTime"]])
                           if (reg != "2022") {
                             extTime.reg <- which(timeseq == extTime)
                           } else {
                             extTime.reg <- 2022
                           }
                           reg_sim <- sim[regional_pops[[reg]], ]
                           stopifnot(nrow(reg_sim) == length(regional_pops[[reg]]))
                           if (!is.null(dim(reg_sim))) {
                             ts <- colSums(reg_sim, na.rm = TRUE)
                           } else {
                             # if region is single cell, nothing to sum
                             ts <- reg_sim
                           }
                           ts_ext <- NULL
                           # time of extinction?
                           # if (all(ts == 0)) {
                           #   ts_ext <- 1
                           #   } else {
                           #     ts_ext <- base::match(0, ts)
                           # }
                           # ext_time <- timeseq[ts_ext]
                           
                           if (all(ts > 0)) {
                             ts_ext <- length(ts)
                           } else {
                             ts_ext <- base::match(0, ts)
                           }
                           if(extTime == "2022"){ext_time <- extTime}else{ext_time <- timeseq[ts_ext]}
                           
                           ext_pen <- fifelse(test =
                                                ## correct timing?
                                                ext_time == timeseq[extTime.reg],
                                              yes = 0, # no penalty
                                              ## overshoot
                                              no = fifelse(test = ext_time > timeseq[extTime.reg],
                                                           yes = ext_time - timeseq[extTime.reg],
                                                           ## undershoot
                                                           no = fifelse(test = ext_time < timeseq[extTime.reg],
                                                                        yes = ext_time - timeseq[extTime.reg],
                                                                        ## extant (or errors)
                                                                        no = NA_real_)))
                           extinct <- ifelse(ts_ext == 271, "Extant", "Extinct")
                           centroid <- lapply(extinct, function(x,...) {
                             if (x == "Extinct") {
                               # location of extinction
                               if (ts_ext == 1) {
                                 return(cbind("X" = NA_real_, "Y" = NA_real_))
                               } else {
                                 # raster one step before regional extinction time
                                 r <- region$raster_from_values(sim[, ts_ext-1])
                                 xy <- region$coordinates[regional_pops[[reg]], ]
                                 rExt <- setDT(raster::extract(r, xy, df = TRUE))
                                 colnames(rExt)[2] <- "layer"
                                 rExt[, `:=`(x = xy[, 1], y = xy[, 2])]
                                 set(x = rExt, i = which(rExt[["layer"]] == 0), j = "layer", value = NA)
                                 rExt <- na.omit(rExt, cols = "layer")
                                 wghtCnt <- weightedCentre(rExt[["x"]], rExt[["y"]], z = rExt[["layer"]])
                                 # find the final location
                                 return(cbind("X" = wghtCnt[,1], "Y" = wghtCnt[,2]))
                               }
                             } else if (x == "Extant") {
                               return(cbind("X" = NA_real_, "Y" = NA_real_))
                             }}
                           )
                           DT <- data.table("SampleID" = sample_data[sample_row, ][["UniqueID"]],
                                            "Region" = reg,
                                            "Extinct" = extinct,
                                            "ExtTS" = ts_ext,
                                            "ExtTime" = ext_time,
                                            "RegionExtTime" = timeseq[extTime.reg],
                                            "ExtPen" = ext_pen,
                                            "ExtX" = centroid[[1]][, 1],
                                            "ExtY" = centroid[[1]][, 2])
                           return(DT)
                         })
                         regional_metrics <- rbindlist(regional_metrics)
                         
                         regional_metrics[, `:=`(AllSites_9 = sites_9,
                                                 SitesProp_9 = sites_9perc,
                                                 AllSites_1 = sites_1,
                                                 SitesProp_1 = sites_1perc,
                                                 Final_abund = final_abund,
                                                 Persist_Target = persistence,
                                                 persistance_dist = ext_dist,
                                                 Persistence_locns = persist_9)]
                         removeTmpFiles(h = 0.16)
                         return(regional_metrics)
                       }
stopCluster(cl)
ext_metrics

saveRDS(ext_metrics, "third_try/validation_targets/validation_targets_10000.RDS", compress = TRUE)



# collapse

ext_metric_summ <- rbindlist(pblapply(seq_along(ext_metrics), function(sample) {
  require(data.table)
  DT <- ext_metrics[[sample]]
  if(!length(DT) == 17){
  DT2 <- data.table(SampleID = DT[1L, ][["SampleID"]],
                    AllSites_9 = NA,
                    SitesProp_9 = NA,
                    AllSites_1 = NA,
                    SitesProp_1 = NA,
                    Final_abund = NA,
                    Persistence_Target = NA,
                    Persistence_locns = NA,
                    Persistence_dist = NA,
                    Ext.1886 = NA,
                    ExtTS.1886 = NA,
                    ExtPen.1886 = NA,
                    Ext.1938 = NA,
                    ExtTS.1938 = NA,
                    ExtPen.1938 = NA,
                    Ext.2003 = NA,
                    ExtTS.2003 = NA,
                    ExtPen.2003 = NA,
                    Ext.2009 = NA,
                    ExtTS.2009 = NA,
                    ExtPen.2009 = NA,
                    Ext.2022 = NA,
                    ExtTS.2022 = NA,
                    ExtPen.2022 = NA)
  }else{
  DT2 <- data.table(SampleID = DT[1L, ][["SampleID"]],
                    AllSites_9 = DT[1L, ][["AllSites_9"]],
                    SitesProp_9 = DT[1L, ][["SitesProp_9"]],
                    AllSites_1 = DT[1L, ][["AllSites_1"]],
                    SitesProp_1 = DT[1L, ][["SitesProp_1"]],
                    Final_abund = DT[1L, ][["Final_abund"]],
                    Persistence_Target = DT[1L, ][["Persist_Target"]],
                    Persistence_locns = DT[1L, ][["Persistence_locns"]],
                    Persistence_dist = DT[1L, ][["persistance_dist"]],
                    Ext.1886 = DT[Region == "1886", ][["Extinct"]],
                    ExtTS.1886 = DT[Region == "1886", ][["ExtTS"]],
                    ExtPen.1886 = DT[Region == "1886", ][["ExtPen"]],
                    Ext.1938 = DT[Region == "1938", ][["Extinct"]],
                    ExtTS.1938 = DT[Region == "1938", ][["ExtTS"]],
                    ExtPen.1938 = DT[Region == "1938", ][["ExtPen"]],
                    Ext.2003 = DT[Region == "2003", ][["Extinct"]],
                    ExtTS.2003 = DT[Region == "2003", ][["ExtTS"]],
                    ExtPen.2003 = DT[Region == "2003", ][["ExtPen"]],
                    Ext.2009 = DT[Region == "2008", ][["Extinct"]],
                    ExtTS.2009 = DT[Region == "2008", ][["ExtTS"]],
                    ExtPen.2009 = DT[Region == "2008", ][["ExtPen"]],
                    Ext.2022 = DT[Region == "2022", ][["Extinct"]],
                    ExtTS.2022 = DT[Region == "2022", ][["ExtTS"]],
                    ExtPen.2022 = DT[Region == "2022", ][["ExtPen"]])
  }

  
  DT2[, ExtPattern := sum(c(c(Ext.1886,Ext.1938,Ext.2003, Ext.2009)  == "Extinct",
                            Ext.2022 == "Extant"))]
  return(DT2)
}), use.names = TRUE)

ext_metric_summ[,  ExtPen.2022 := 271-ExtTS.2022]
ext_metric_summ[, ExtPen := sqrt(((ExtPen.1886^2) + (ExtPen.1938^2) + (ExtPen.2003^2) + (ExtPen.2009^2) + (ExtPen.2022^2))/5)] #RMSE of extinction penalty

anyNA(ext_metric_summ)


saveRDS(ext_metric_summ,"third_try/validation_targets/validation_targets_summarised.RDS",compress = TRUE)


#ext_metric_summ <- readRDS("third_try/validation_targets/validation_targets_summarised.RDS")

ext_metric_summ

# Count simulations where Persistence_locns > 0 for 2020
persist_locns_2020 <- ext_metric_summ[Persistence_locns > 0 & Ext.2022 == "Extant", ]
n_simulations <- nrow(persist_locns_2020)
n_simulations
sample_ids <- persist_locns_2020$SampleID
sample_ids

# Initialize vectors to store abundances
location_dryandra_abundances <- c()
location_perup_abundances <- c()

# Loop through each simulation and determine abundance at each final location
for (id in sample_ids) {
  # Load the corresponding simulation results
  sim_path <- file.path(sim_dir, paste0("UniqueID_", id, "_results.RData"))
  if (file.exists(sim_path)) {
    sim <- readRDS(sim_path)$abundance[, 271]  # Extract the last timestep
    
    # Convert to raster
    sim_raster <- region$raster_from_values(sim)
    
    # Get cell indices for each final location
    dryandra_cell <- cellFromXY(sim_raster, st_coordinates(final_locations[1, ]))
    perup_cell <- cellFromXY(sim_raster, st_coordinates(final_locations[2, ]))
    
    # Extract abundance values
    dryandra_abund <- sim_raster[dryandra_cell]
    perup_abund <- sim_raster[perup_cell]
    
    # Only store non-zero abundances
    if (!is.na(dryandra_abund) && dryandra_abund > 0) {
      location_dryandra_abundances <- c(location_dryandra_abundances, dryandra_abund)
    }
    if (!is.na(perup_abund) && perup_abund > 0) {
      location_perup_abundances <- c(location_perup_abundances, perup_abund)
    }
  }
}

# Output results
cat("Number of simulations reaching Dryandra:", length(location_dryandra_abundances), "\n")
cat("Number of simulations reaching Perup:", length(location_perup_abundances), "\n")

if (length(location_dryandra_abundances) > 0) {
  cat("Mean abundance at Dryandra:", mean(location_dryandra_abundances, na.rm = TRUE), "\n")
} else {
  cat("No simulations reached Dryandra.\n")
}

if (length(location_perup_abundances) > 0) {
  cat("Mean abundance at Perup:", mean(location_perup_abundances, na.rm = TRUE), "\n")
} else {
  cat("No simulations reached Perup.\n")
}


# If you want more details, you can print the full vectors
print(location_dryandra_abundances)
print(location_perup_abundances)







