library(data.table)
library(doSNOW)
library(gtools)
library(hypervolume) # v 3.0.1 +
library(raster)
library(pbapply)
library(rnaturalearth)
library(sf)
library(ggplot2)

# Load helper functions
source("F:/Box Sync/students/Audrey Praz/code/helper_funs.R")

# Define directories
data_dir <- "F:/Box Sync/students/Audrey Praz/third_try/sim_inputs"
cut_dir <- "F:/Box Sync/students/Audrey Praz/data/niche_cuts"
output_dir <- "F:/Box Sync/students/Audrey Praz/data/hypervolumes"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# STEP 1: BUILD HYPERVOLUMES ####
# Load the numbat and climate data
matched <- fread("data/numbat_data/Exported Numbat data/full_matched_records.csv")
matched
colnames(matched) <- c("UID", "Spp", "Name", "Year", "Lat", "Lon", "Type", 
                       "WQT", "CQT", "PREC")
matched

# Subset to relevant climatic variables (Warmest quarter temp, Coolest quarter temp, Precipitation)
vars <- colnames(matched)[c(8:10)]
vars

# Center and scale measurements
df_Msub <- scale(matched[, ..vars], center = TRUE, scale = TRUE)
df_Msub.bind <- cbind(df_Msub, matched[, c(2,4:7)])

# Estimate bandwidth for kernel density estimation
bw <- estimate_bandwidth(df_Msub, method = "silverman-1d")

# Construct the hypervolume using Gaussian KDE (Kernel Density Estimation)
climate_hv.Msub <- hypervolume(data = df_Msub,
                               name = "full_hypervolume",
                               verbose = TRUE,
                               method = "gaussian",
                               kde.bandwidth = estimate_bandwidth(data = df_Msub,
                                                                  method = "fixed",
                                                                  value = bw))
# Visualize the hypervolume
plot(climate_hv.Msub, show.legend = FALSE)
# Save data to output directory  
saveRDS(climate_hv.Msub, file.path(output_dir, "ClimateHypervolume.RDS"), compress = TRUE)
saveRDS(df_Msub, file.path(output_dir, "ClimateDF.RDS"), compress = TRUE)


# STEP 2: PROJECTION AND SPATIAL CALCULATION ####
# Climatic Projections

climate_hv.Msub <- readRDS(file.path(output_dir, "ClimateHypervolume.RDS"))
df_Msub <- readRDS(file.path(output_dir, "ClimateDF.RDS"))

# Load raster mask and climatic projections
ras_mask <- raster("F:/Box Sync/students/Audrey Praz/data/template_full_20km.grd")
WQT <- brick("F:/Box Sync/students/Audrey Praz/data/climate/HARM_warmest_qtr.grd")
CQT <- brick("F:/Box Sync/students/Audrey Praz/data/climate/HARM_coolest_qtr.grd")
PREC <- brick("F:/Box Sync/students/Audrey Praz/data/climate/HARM_annual_rainfall.grd")

# Create a parallel processing cluster (4 cores)
clProj <- makeCluster(4L)
clusterExport(clProj, c("WQT", "CQT", "PREC", "ras_mask", "df_Msub", "climate_hv.Msub"))

# Project hypervolume to each raster layer (time step)
HV_seq <- pblapply(seq_len(nlayers(WQT)), function(i) {
  require("hypervolume"); require("raster")
  s <- mask(stack(WQT[[i]], CQT[[i]], PREC[[i]]), ras_mask)
  s <- scale(s, center = attr(df_Msub,"scaled:center"), scale = attr(df_Msub,"scaled:scale"))
  hv_proj <- hypervolume_project(hv = climate_hv.Msub, rasters = s,
                                 type = "probability",
                                 parallel = TRUE, n.cores = 3L,
                                 verbose = FALSE,
                                 set.edges.zero = TRUE,
                                 edges.zero.distance.factor = 25,
                                 weight.exponent = -3
                                 )
  plot(hv_proj)
  return(hv_proj)
}, cl = clProj)
stopCluster(clProj)

# Stack all projected hypervolumes into a single raster object
HV_seq <- stack(HV_seq, quick = TRUE)
# Assign names for each time steps
names(HV_seq) <- paste0("X.", seq(1750, 2020, by = 1))
plot(HV_seq[[seq(1, nlayers(HV_seq), l = 6)]])

# Calculate the 95th percentile and rescale to {0,1}
xVal <- values(HV_seq)
q95 <- quantile(xVal[xVal > 0], 0.95, na.rm = TRUE)
HV_seqR <- stretch(x = HV_seq, minv = 0, maxv = 1, smin = 0, smax = q95)
plot(HV_seqR[[seq(1, nlayers(HV_seqR), l = 6)]], zlim = c(0, 1),
     col = hcl.colors(100), legend = FALSE, nr = 3)

# Save the final raster of projected suitability
writeRaster(HV_seqR, file.path(output_dir, "Numbat_suitability.grd"), 
            overwrite = TRUE)


# STEP 3: NICHE CUTS ####

## Generating samples
## widths of the "cutting box"
widths <- seq(0.8, 1.0, by = 0.05) ## cuts from 0.8 to 1.0 by 0.05
id_cols <- names(matched)[c(1:6)]

## cutting
pb <- txtProgressBar(max = length(widths), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

cls <- length(widths)
cls <- makeCluster(cls)
registerDoSNOW(cls)


foreach(width = widths,
        .inorder = TRUE,
        .options.snow = opts,
        .packages = c("data.table")) %dopar% {
          fracwidth <- width
          fracpnt <- seq(0 + (fracwidth/2), 1 - (fracwidth/2), by = 0.05)
          loop <- 0
          ## this is for 3 climatic variables
          for(x in seq_along(fracpnt)){
            for (y in seq_along(fracpnt)){
              for(z in seq_along(fracpnt)) {
                v1 <- subset_data(as.data.frame(matched[[vars[1]]]), fracpnt[x], fracwidth)
                v2 <- subset_data(as.data.frame(matched[[vars[2]]]), fracpnt[y], fracwidth)
                v3 <- subset_data(as.data.frame(matched[[vars[3]]]), fracpnt[z], fracwidth)
                v <- v1 * v2 * v3
                if(all(v == 0)){
                  next
                } else {
                  loop <- loop + 1
                  sample <- cbind(copy(matched), v)
                  sample_out <- copy(sample)[v == 1, ]
                  sample_out[,`:=`(width = fracwidth,
                                   x = fracpnt[x],
                                   y = fracpnt[y],
                                   z = fracpnt[z],
                                   file_rows = .N,
                                   file = loop)]
                  fwrite(sample_out, 
                         file = paste0(cut_dir, sprintf("/cut_%s", fracwidth), sprintf("_%04d.csv", loop)))
                }
              }
            }
          }
        }
stopCluster(cls); registerDoSEQ(); gc()

## read files and check for duplicates
## the duplicates are defined based on the UniqueID, which is a number defining the specific set of values
## the UniqueID was assigned to each row of the dataset resulting from the pairing of fossils with the climate. Because
## there are no duplicated sets of climatic variable, each set has assigned a UniqueID.
## If two niche cuts are identical, they have the same UniqueIDs (and therefore also the same number of rows), so we check
## for UniqueIDs and remove the ones that are identical.

cuts <- list.files(cut_dir, "\\.csv$", full.names = TRUE)
idx <- lapply(cuts, function(x) {
  DT <- fread(x)
  DT[["UID"]]
})
idx <- which(!duplicated(idx))

cuts <- cuts[idx]

## read in only the unique files
nodup <- rbindlist(pblapply(cuts, fread))
setorder(nodup, file, width, UID)
nodup

# we remove cuts with less than 20 records
nodup2 <- nodup[ , `:=`(COUNT = .N , IDX = 1:.N) , by = c("width", "file")][
  , .SD[COUNT >= 20], keyby = .(width, file)]
setcolorder(nodup2, 
            neworder = c("UID", "Spp", "Year", "PREC", "WQT", "CQT", "width", 
                         "x", "y", "z", "file_rows", "file"))
setorder(nodup2, file, width, UID)
nodup2

saveRDS(nodup2, file.path(cut_dir, "unique_cuts.RDS"))

# Projecting cuts
full_hv <- readRDS(file.path(output_dir, "ClimateHypervolume.RDS"))

# process each of the unique hypervolumes
cls <- makeCluster(12L)
registerDoSNOW(cls)
pb <- txtProgressBar(max = length(cuts), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
foreach(i = seq_along(cuts),
          .options.snow = opts, .errorhandling = "remove",
          .packages = c("data.table", "hypervolume"),
          .export = c("cuts", "full_hv"),
          .noexport = c(paste(ls(), collapse = ","))) %dopar% {
            # Read in the csv
            DT <- fread(cuts[i], select = colnames(full_hv@Data))
            DT <- as.matrix(DT)
            ## scale relative to full HV
            DT <- scale(DT, 
                        center = attr(full_hv@Data, "scaled:center"), 
                        scale = attr(full_hv@Data, "scaled:scale"))
            ## build the gaussian hypervolume
            hv_name <- gsub(".csv", "", sapply(strsplit(cuts[i], "/"), tail, 1))
            hv <- hypervolume_gaussian(DT, name = hv_name, chunk.size = 1000,
                                       verbose = TRUE, sd.count = 3,
                                       kde.bandwidth = estimate_bandwidth(data = DT, method = "fixed",
                                                                          value = full_hv@Parameters$kde.bandwidth)
                                       )
            saveRDS(hv, file = paste0(file.path(output_dir, hv_name), ".RDS"),
                    compress = TRUE)
          }
registerDoSEQ(); gc()
stopCluster(cls)


# Climatic Projections for each "cut" hypervolume
cut_hv <- list.files(output_dir, "\\.RDS$", full.names = TRUE)
cut_hv <- cut_hv[grepl("cut_", cut_hv)]
length(cut_hv) == length(cuts)  

for (cut in cut_hv) {
  message(cut)
  if(file.exists(gsub("RDS", "grd", gsub("hypervolumes", "niche_cuts", cut)))){next}
  hv_cut <- readRDS(cut)
  clProj <- makeCluster(30L)
  clusterExport(clProj, c("WQT", "CQT", "PREC", "ras_mask", "hv_cut"))
  HV_seq <- pblapply(seq_len(nlayers(WQT)), function(i) {
    require("hypervolume"); require("raster")
    s <- mask(stack(WQT[[i]], CQT[[i]], PREC[[i]]), ras_mask)
    # data in the hypervolume contains the center/scaling values from the full_hv
    s <- scale(s, center = attr(hv_cut@Data,"scaled:center"), scale = attr(hv_cut@Data,"scaled:scale"))
    hv_proj <- hypervolume_project(hv = hv_cut, rasters = s,
                                   type = "probability",
                                   parallel = FALSE, # = 3L,
                                   verbose = FALSE,
                                   set.edges.zero = TRUE,
                                   edges.zero.distance.factor = 25,
                                   weight.exponent = -3)
    return(hv_proj)
  }, cl = clProj)
  stopCluster(clProj)
  HV_seq <- stack(HV_seq, quick = TRUE)
  names(HV_seq) <- paste0("X.", seq(1750, 2020, by = 1))
  # Calculate the 95th percentile and rescale to {0,1}
  xVal <- values(HV_seq)
  q95 <- quantile(xVal[xVal > 0], 0.95, na.rm = TRUE)
  HV_seqR <- stretch(x = HV_seq, minv = 0, maxv = 1, 
                     smin = min(xVal, na.rm = TRUE), smax = q95)
  ras_out <- gsub("RDS", "grd", gsub("hypervolumes", "niche_cuts", cut))
  writeRaster(HV_seqR, ras_out, overwrite = TRUE)
}

# STEP 4: PLOTTING A SUMMARY FIGURE ####

# Import a land-sea mask to account for broader projection space for Numbat habitat suitability
landsea <- raster("data/numbat_region_20km.grd")

# Import habitat suitability projections
hs_numbat <- brick("data/Numbat_suitability.grd")
hs_numbat <- mask(hs_numbat, landsea)
names(hs_numbat) <- paste("X", seq(1750, 2020, 1))
hs_numbat

plot(landsea, colNA = "blue")
plot(hs_numbat[[seq(1, nlayers(hs_numbat), l = 6)]], colNA = "blue4")

stack_df <- setDT(as.data.frame(hs_numbat, xy = TRUE, na.rm = TRUE))
stack_df <- melt(stack_df, id.vars = c("x", "y"), variable.name = "layer",
                 value.name = "hs")
stack_df

# import Australian coastline shapefile
aus <- ne_countries(scale = 50, country = "australia", returnclass = "sf")
aus <- st_crop(aus, extent(hs_numbat))
aus; plot(aus[1])

revmako <- rev(viridis::mako(20))
revvirid <- rev(viridis(20))


final_pts <- fread("D:/Box Sync/Box Sync/students/Audrey Praz/GIS stuff/translocation_sites.csv")
final_pts <- final_pts[type == "Translocation sites"]
final_pts_sf <- st_as_sf(final_pts, coords = c("x", "y"), crs = 4326)
final_pts_sf$Subpop <- "Translocation sites"


p1 <- ggplot() +
  geom_sf(data = aus, inherit.aes = FALSE, fill = 'grey80', col = "grey30") +
  geom_raster(data = stack_df[layer == "X.2020", ],
              aes(x = x, y = y, fill = hs)) +
  scale_fill_gradientn(colors=revmako, 
                       breaks=seq(0.4,1,0.2), limits=c(0.4,1),
                       na.value = "grey85",
                       guide = guide_colourbar(title = "Suitability",
                                               title.position = "top",
                                               title.theme = element_text(size = 8, colour = "black", face='bold'),
                                               title.hjust = 0.5,
                                               title.vjust = 5,
                                               label = TRUE,
                                               label.position = "left",
                                               barwidth = 0.5,
                                               #barheight = 15,
                                               ticks = FALSE,
                                               direction = "vertical",
                                               frame.colour = "black")) +
  # scale_fill_viridis_c(breaks = seq(0.4, 1, 0.1),
  #                      limits = c(0.4, 1),
  #                      na.value = "grey85",
  #                      option='G',
  #                      guide = guide_colourbar(title = "Suitability",
  #                                              title.position = "top",
  #                                              title.theme = element_text(size = 8, colour = "black", face='bold'),
  #                                              title.hjust = 0.5,
  #                                              title.vjust = 5,
  #                                              label = TRUE,
  #                                              label.position = "left",
  #                                              barwidth = 0.5,
  #                                              #barheight = 15,
  #                                              ticks = FALSE,
  #                                              direction = "vertical",
  #                                              frame.colour = "black")) +
  geom_sf(data = aus, inherit.aes = FALSE, fill = NA, col = "grey30", linewidth=1) +
  geom_sf(data = final_pts_sf, shape=21, color= "black", fill = "firebrick1", size = 4, show.legend = TRUE) +  # <-- Plot points after all else
  theme_bw() +
  theme(legend.position = "right", axis.text = element_blank()) +
  labs(x = NULL, y = NULL) +
  ggtitle(paste("a) hypervolume projections c. 2020 CE"))


p1

# define time frame
time <- seq(1750, 2020, 1)

#### create dataframe for suitability trend ####
sp1 <- hs_numbat
sp1[sp1 == 0] <- NA
aw_noz1 <- mask(area(sp1), sp1)
weighted_hs1 <- sp1 * aw_noz1

# area weighted mean and SD
ts_hs1 <- cellStats(weighted_hs1, sum, na.rm = TRUE) / cellStats(aw_noz1, sum, na.rm = TRUE)
sd1 <- sapply(1:nlayers(sp1), function(x,...) sqrt(Hmisc::wtd.var(values(sp1[[x]]), weights = values(aw_noz1[[x]]), na.rm = TRUE)))

# adding range area trend
sp1[!is.na(sp1)] <- 1 ## everything else to 1

r_hs1 <- cellStats(sp1, sum) ## sum == count of cells that are suitable

## create data.table for plotting
ts_hs <- data.table(Suitability = ts_hs1,
                    Year = time,
                    SD = sd1,
                    Range = r_hs1,
                    DB = "Numbat")
ts_hs

ts_hs$rangeSeanUnits <- ts_hs$Range * 0.4

# plot suitability trend
p3 <- ggplot() + 
  geom_ribbon(data = ts_hs, 
              aes(x = Year, ymin = fifelse(Suitability - SD < 0, 0, Suitability - SD), 
                  ymax = fifelse(Suitability + SD > 1, 1, Suitability + SD), 
                  fill = DB, col = DB), alpha = 1/3,
              show.legend = FALSE) +
  geom_line(data = ts_hs,
            aes(x = Year, y = Suitability, col = DB),
            show.legend = FALSE) +
  theme_bw() +
  xlab("Year BP") +
  scale_x_continuous(limits = range(time), 
                     breaks = seq(1750, 2020, 20)) +
  scale_y_continuous(limits = c(0, 1)) +
  ggtitle("b) Habitat suitability +/- SD")

# plot range area
p4 <- ggplot() + 
  geom_line(data = ts_hs,
            aes(x = Year, y = (Range*400)/1000, col = DB), linewidth = 1.5) + 
  xlab("Year BP") +
  ylab("Range 10^4 (sq. km)") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_x_continuous(limits = range(time), 
                     breaks = seq(1750, 2020, 20)) +
  scale_y_continuous(limits = range((ts_hs[["Range"]]*400)/1000),
                     labels = function(x) formatC(x, big.mark = ",", digits = 9)) +
  ggtitle("c) Range area")



# plot hypervolume in 3D
# Extract hypervolume data points
extract_hypervolume_points <- function(hypervolume) {
  data <- as.data.frame(hypervolume@RandomPoints)
  colnames(data) <- c("WQT", "CQT", "PREC")  # Rename columns for clarity
  return(data)
}

# Get points for hv_old and hv_new
hv_points <- extract_hypervolume_points(climate_hv.Msub)


# Add group labels
hv_points$Group <- "Numbat Data"




# Plot with plotly
fig <- plot_ly(
  data = hv_points,
  x = ~WQT, y = ~CQT, z = ~PREC,
  color = ~Group,
  colors = c("aquamarine4"),
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3, opacity = 0.7)  # Adjust size and opacity for clarity
)

# Add layout details
fig <- fig %>%
  layout(
    scene = list(
      xaxis = list(title = "WQT"),
      yaxis = list(title = "CQT"),
      zaxis = list(title = "Annual Rainfall")
    ),
    title = "3D Visualization of Hypervolume"
  )

# Display the plot
fig





