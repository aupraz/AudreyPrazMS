#### Libraries ####
library(data.table)
library(paleopop)
library(pbapply)
library(poems)
library(raster)
library(rnaturalearth)
library(sf)

#### Workflow options ####
parallel_cores <- 30L # cluster: 20 ?
nsims <- 25000

str(readRDS("numbat_temporal_region.RDS"), max.level = 2)

burn_in_steps <- 50 # 25 generations
timesteps <- 271 + burn_in_steps ## to start at 30k BP
plot_seq <- floor(seq(burn_in_steps+1, timesteps, length.out = 6))

# location of HS projections
data_dir <- "data/niche_cuts/proj_matrices/projected"

# output directory
out_dir <- "baseline/out_sims"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# output directory for test simulations
test_dir <- "baseline/out_sims/test_sims"

# time seq not including burnin
timeseq <- seq(to = 2020, by = 1, length.out = (timesteps - burn_in_steps))

#### Step 1: Create a paleopop region ####

# This step generates a spatio-temporally epxlicit study region #

# Temporal mask
## matrix of grid cells (50301) x time (271), in projected coordinates
## see str(temp_mask) for metadata
temp_mask <- readRDS("data/numbat_temporal_region.RDS")
str(temp_mask)

# index for timesteps
## we are only interested in modelling these timesteps
tidx <- which(attr(temp_mask, "timesteps") %in% timeseq)

# Generate a projected template regional raster
templ_ras <- raster(res = 0.20,
                    ext = extent(111.05,155.05,-45.025,-9.025),
                    crs = crs(raster()))
values(templ_ras) <- 1
proj_ext <- projectExtent(templ_ras, crs = attr(temp_mask, "proj_wkt"))
res(proj_ext) <- attr(temp_mask, "res")
all.equal(coordinates(proj_ext)[,1],  attr(temp_mask, "X")) ## TRUE
all.equal(coordinates(proj_ext)[,2],  attr(temp_mask, "Y")) ## TRUE
proj_ext

TEST_REGION <- TRUE
if (TEST_REGION) {
  templ_rasprj <- projectRaster(templ_ras, proj_ext, method = "ngb")
  plot(templ_rasprj, legend = FALSE)
  # values of template raster == maximum extent of land
  ## ALL HS projections are a subset of this area
  max_ext <- apply(temp_mask, 1, max, na.rm = TRUE)
  templ_rasprj[] <- max_ext
  templ_rasprj[is.infinite(templ_rasprj)] <- NA
  templ_rasprj
  plot(templ_rasprj, legend = FALSE)
}

# Produce a temporally explicit template raster using the temporal mask
templ_ras <- brick(nl = timesteps, crs = crs(proj_ext),
                   xmn = -2423364, xmx = 2436636,
                   ymn = -2158270, ymx = 1981730,
                   nrow = 207, ncol = 243)
templ_ras

# repeat the first layer of the mask to the number of burn_in_steps
temp_mask_w_burnin <- cbind(replicate(burn_in_steps, temp_mask[, tidx[1]]), temp_mask[, tidx])
values(templ_ras) <- temp_mask_w_burnin
plot(templ_ras[[plot_seq]],
     col = c("#FFFFFF", "#a3be8c"), colNA = "#bf616a",legend = FALSE,
     ext = extent(templ_ras))

# PaleoRegion from templ_ras
region <- PaleoRegion$new(template_raster = templ_ras,
                          remove_zeros = TRUE,
                          use_raster = TRUE)
saveRDS(region, "third_try/sim_inputs/numbat_paleopopregion2.RDS", compress = TRUE)
region$region_raster

# Continents
conts <- ne_countries(scale = 50, returnclass = "sf")
conts <- conts[conts$name == "Australia", ]
conts <- st_transform(conts, crs = st_crs(region$region_raster))
conts <- st_crop(conts, extent(region$region_raster))

plot(region$region_raster,
     col = hcl.colors(ncell(region$region_raster), "Lajolla"),
     legend = FALSE,
     ext = extent(templ_ras),
     addfun = function() lines(as_Spatial(conts)))

compareRaster(region$region_raster, templ_ras,
              extent = TRUE, rowcol = TRUE, crs = TRUE,
              res = TRUE, orig = TRUE, rotation = TRUE,
              values = FALSE)
region$region_cells # 9242 cells

# temporal mask
## Green = suitable; Grey = not suitable; Red = never suitable
plot(region$raster_from_values(region$temporal_mask[, plot_seq]),
     colNA = "#bf616a", col = c("#FFFFFF", "#a3be8c"),
     axes = FALSE,
     legend = FALSE,
     addfun = function() lines(as_Spatial(conts)))

#### Step 2: Distance matrix and environmental correlation ####
# Distance matrix
## pre-calculate once, save, then read in
distance_matrix <- DispersalGenerator$new(region = region,
                                          dispersal_max_distance = 50,
                                          # distance in km
                                          distance_scale = 1000)
distance_matrix <- distance_matrix$calculate_distance_matrix(use_longlat = FALSE)
dim(distance_matrix); summary(as.vector(distance_matrix))
saveRDS(distance_matrix, "third_try/sim_inputs/numbat_distance_matrix_projected2.RDS")

## read in the pre-calculated distance matrix
distance_matrix <- readRDS("third_try/sim_inputs/numbat_distance_matrix_projected2.RDS")
dim(distance_matrix); summary(as.vector(distance_matrix))

# Compact decomposition
## Distance-based environmental correlation (via a compacted Cholesky decomposition)
## as before, pre-calculate once, save, then read in
env_corr <- SpatialCorrelation$new(region = region,
                                   amplitude = 0.99,
                                   breadth = 100,
                                   # distances in km
                                   distance_scale = 1000)
env_corr$calculate_correlations(distance_matrix = distance_matrix) # all decimals
env_corr$calculate_cholesky_decomposition(decimals = 3)
compact_decomposition <- env_corr$get_compact_decomposition() # default threshold
saveRDS(compact_decomposition, "baseline/sim_inputs/numbat_compact_decomposition_projected.RDS")

## read in the pre-calculated compact decomposition
compact_decomposition <- readRDS("baseline/sim_inputs/numbat_compact_decomposition_projected.RDS")
env_corr$t_decomposition_compact_matrix <- compact_decomposition$matrix
env_corr$t_decomposition_compact_map <- compact_decomposition$map

#### Step 3a: Cat arrival times and abundance bounds ####
cat_arrival <- readRDS("data/cat_predation_layers/cat_arrival_temporal_mask.RDS")
cat_abund <- readRDS("data/cat_predation_layers/cat_abundance_estimate.RDS")
cat_abund[] <- cat_abund[]/attr(cat_abund, "scale.fac")
dim(cat_arrival)
dim(cat_arrival) == dim(temp_mask)

# replicate arrival times for the burn-in
cat_arrival <- cbind(replicate(burn_in_steps, cat_arrival[, tidx[1]]), cat_arrival[, tidx])
dim(cat_arrival) == dim(temp_mask_w_burnin)
any(is.infinite(cat_arrival[])); any(is.infinite(cat_abund[])) ## FALSE FALSE

# subset cats to the region and replace NA with 0
cat_arrival_matrix <- cat_arrival[region$region_indices, ]
cat_arrival_matrix[is.na(cat_arrival_matrix)] <- 0

cat_abund_matrix <- cat_abund[region$region_indices, ]
cat_abund_matrix[is.na(cat_abund_matrix)] <- 0

# resample cat predation mean, upper, and lower to the correct
# number of time steps
cat_abund_mean <- replicate(n = ncol(cat_arrival), cat_abund_matrix[, 2])*cat_arrival_matrix
cat_abund_lower <- replicate(n = ncol(cat_arrival), cat_abund_matrix[, 1])*cat_arrival_matrix
cat_abund_upper <- replicate(n = ncol(cat_arrival), cat_abund_matrix[, 3])*cat_arrival_matrix

dim(cat_arrival_matrix); dim(cat_abund_mean)

plot(region$raster_from_values(cat_arrival_matrix[, plot_seq]),
     colNA = "#bf616a", col = c("#FFFFFF", "#a3be8c"),
     axes = FALSE,
     legend = FALSE,
     addfun = function() lines(as_Spatial(conts)))

plot(region$raster_from_values(cat_abund_mean[, plot_seq]),
     breaks = seq(0, 1, 0.1),
     col = hcl.colors(11, "RdYlGn", rev = TRUE),
     axes = FALSE,
     legend = TRUE,
     addfun = function() lines(as_Spatial(conts)))

#### Step 3b: fox arrival times and abundance bounds ####
fox_arrival <- readRDS("data/fox_invasion_abundance/fox_arrival_temporal_mask.RDS")
fox_abund <- readRDS("data/fox_invasion_abundance/fox_abundance_estimate.RDS")
fox_abund[] <- fox_abund[]/attr(fox_abund, "scale.fac")
dim(fox_arrival)
dim(fox_arrival) == dim(temp_mask)

# replicate arrival times for the burn-in
fox_arrival <- cbind(replicate(burn_in_steps, fox_arrival[, tidx[1]]), fox_arrival[, tidx])
dim(fox_arrival) == dim(temp_mask_w_burnin)
any(is.infinite(fox_arrival[])); any(is.infinite(fox_abund[])) ## FALSE FALSE

# subset cats to the region and replace NA with 0
fox_arrival_matrix <- fox_arrival[region$region_indices, ]
fox_arrival_matrix[is.na(fox_arrival_matrix)] <- 0

fox_abund_matrix <- fox_abund[region$region_indices, ]
fox_abund_matrix[is.na(fox_abund_matrix)] <- 0

# resample fox predation mean, upper, and lower to the correct
# number of time steps
fox_abund_mean <- replicate(n = ncol(fox_arrival), fox_abund_matrix[, 2])*fox_arrival_matrix
fox_abund_lower <- replicate(n = ncol(fox_arrival), fox_abund_matrix[, 1])*fox_arrival_matrix
fox_abund_upper <- replicate(n = ncol(fox_arrival), fox_abund_matrix[, 3])*fox_arrival_matrix

dim(fox_arrival_matrix); dim(fox_abund_mean)

plot(region$raster_from_values(fox_arrival_matrix[, plot_seq]),
     colNA = "#bf616a", col = c("#FFFFFF", "#a3be8c"),
     axes = FALSE,
     legend = FALSE,
     addfun = function() lines(as_Spatial(conts)))

plot(region$raster_from_values(fox_abund_mean[, plot_seq]),
     breaks = seq(0, 1, 0.1),
     col = hcl.colors(11, "RdYlGn", rev = TRUE),
     axes = FALSE,
     legend = TRUE,
     addfun = function() lines(as_Spatial(conts)))


#### Step 4: Model template ####

# Population (simulation) model template for fixed parameters
model_template <- PaleoPopModel$new(
  region = region,
  time_steps = timesteps, # include burn-in
  years_per_step = 1, # generation length
  populations = region$region_cells,
  # initial_abundance: generated based on HS
  transition_rate = 1.0,
  # standard_deviation: sampled in LHS
  compact_decomposition = compact_decomposition,
  # carrying_capacity: generated based on HS
  density_dependence = "logistic",
  # growth_rate_max: sampled in LHS
  predation = TRUE,
  # predation_max: sampled in LHS
  predation_g = 0.4, # constant
  # predation_z: sampled in LHS
  # predation_max_n: sampled in LHS
  # human_density: generated
  dispersal_target_k = NULL,
  # dispersal_data: generated
  results_selection = c("abundance", "predated"),
  attribute_aliases = list(density_max = "predation_max_n")
)
model_template$is_complete() ## FALSE

#### Step 5: Create generators for initial abundance, carrying capacity, and dispersal ####
#FOR 1 NICHE CUT, WITHOUT LANDUSE#
niche_lookup <- list.files("data/niche_cuts/proj_matrices/projected", ".RDS")
niche_lookup <- data.table(niche_ref = gsub(".RDS", "", niche_lookup))
niche_lookup[, Width := as.numeric(sapply(strsplit(niche_ref, "_"), head, 1))]
niche_lookup

niche_lookup <- as.data.frame(niche_lookup)

capacity_gen <- Generator$new(description = "capacity",
                              region = region,
                              generate_rasters = FALSE,
                              burn_in_steps =  burn_in_steps,
                              generative_requirements = list(
                                hs_matrix = "file",
                                initial_abundance = "function",
                                carrying_capacity = "function"
                              ),
                              inputs = c("density_max", "niche_ref"),
                              ## use the generator to calculate
                              ## initial abundance and carrying capacity
                              outputs = c("initial_abundance", "carrying_capacity"))

## Here we tell the generator to import the HS file and save it as "hs_matrix"
capacity_gen$add_file_template("hs_matrix",
                               path_template = file.path(data_dir, "%s.RDS"),
                               path_params = "niche_ref",
                               file_type = "RDS")

## Here we subset the hs_matrix to only the region cells,
## and then add the burn in.
## Generate the carrying_capacity based on "density_max" and "hs_matrix".
capacity_gen$add_function_template("carrying_capacity",
                                   function_def = function(params) {
                                     hs_matrix <- params$hs_matrix[params$region$region_indices, ]
                                     # HS values stored as int. Divide by 1000.
                                     hs_matrix[] <- hs_matrix[]/1000
                                     # any infinites to 0
                                     hs_matrix[!is.finite(hs_matrix)] <- 0
                                     # append the first timestep n times as a burn in
                                     hs_matrix <- cbind(replicate(params$burn_in_steps,
                                                                  hs_matrix[, 1]),
                                                        hs_matrix[, 1:ncol(hs_matrix)])
                                     # floor the density values to integer
                                     floor(params$density_max*hs_matrix)
                                   },
                                   call_params = c("density_max", "hs_matrix", "burn_in_steps", "region"))

## Here we tell the generator what function to use to generate initial_abundance
## based on the carrying capacity of the first time step
capacity_gen$add_function_template("initial_abundance",
                                   function_def = function(params) {
                                     params$carrying_capacity[, 1]
                                   },
                                   call_params = c("carrying_capacity"))

# WITH LANDUSE (see script above to run the capacity generator without landuse)
# Carrying-capacity generator
landuse_ras <- stack("data/Landuse_change/landuse_ras.tif")
landuse <- readRDS("data/Landuse_change/lu_tempq50.RDS")
dim(landuse)

niche_lookup <- list.files("data/niche_cuts/proj_matrices/projected", ".RDS")
niche_lookup <- data.table(niche_ref = gsub(".RDS", "", niche_lookup))
niche_lookup <- as.data.table(niche_lookup)
niche_lookup[, Width := as.numeric(sapply(strsplit(niche_ref, "_"), head, 1))]
#niche_lookup[, Width := 1.0]
niche_lookup

niche_lookup <- as.data.frame(niche_lookup)

capacity_gen <- Generator$new(description = "capacity",
                              region = region,
                              lu_matrix = landuse,
                              generate_rasters = FALSE,
                              burn_in_steps =  burn_in_steps,
                              generative_requirements = list(
                                hs_matrix = "file",
                                lu_matrix = "file",
                                initial_abundance = "function",
                                carrying_capacity = "function"
                              ),
                              inputs = c("density_max", "niche_ref"),
                              ## use the generator to calculate
                              ## initial abundance and carrying capacity
                              outputs = c("initial_abundance", "carrying_capacity"))

## Here we tell the generator to import the HS file and save it as "hs_matrix"
capacity_gen$add_file_template("hs_matrix",
                               path_template = file.path("data/niche_cuts/proj_matrices/projected", "%s.RDS"),
                               path_params = "niche_ref",
                               file_type = "RDS")
capacity_gen$add_file_template("lu_matrix",
                               path_template = file.path("data/Landuse_change/lu_tempq50.RDS"),
                               file_type = "RDS")
## Here we subset the hs_matrix to only the region cells,
## and then add the burn in.
## Generate the carrying_capacity based on "density_max" and "hs_matrix".
capacity_gen$add_function_template("carrying_capacity",
                                   function_def = function(params) {
                                     hs_matrix <- params$hs_matrix[params$region$region_indices, ]
                                     lu_matrix <- params$lu_matrix[params$region$region_indices, ]
                                     # HS values stored as int. Divide by 1000.
                                     hs_matrix[] <- hs_matrix[]/1000
                                     # any infinites to 0
                                     hs_matrix[!is.finite(hs_matrix)] <- 0
                                     hs_matrix <- hs_matrix*(1-lu_matrix)
                                     # append the first timestep n times as a burn in
                                     hs_matrix <- cbind(replicate(params$burn_in_steps,
                                                                  hs_matrix[, 1]),
                                                        hs_matrix[, 1:ncol(hs_matrix)])
                                     # floor the density values to integer
                                     floor(params$density_max*hs_matrix)
                                   },
                                   call_params = c("density_max", "hs_matrix", "lu_matrix",
                                                   "burn_in_steps", "region"))

## Here we tell the generator what function to use to generate initial_abundance
## based on the carrying capacity of the first time step
capacity_gen$add_function_template("initial_abundance",
                                   function_def = function(params) {
                                     params$carrying_capacity[, 1]
                                   },
                                   call_params = c("carrying_capacity"))


system.time({test_capacity <- capacity_gen$generate(input_values = list(
  density_max = 700,
  # random sample for testing
  niche_ref = {set.seed(961);niche_lookup[sample(x = nrow(niche_lookup), size = 1), 1]}
))
})

lapply(test_capacity, anyNA)
summary(test_capacity$initial_abundance)
summary(as.vector(test_capacity$carrying_capacity))
plot(region$raster_from_values(test_capacity$initial_abundance),
     col = rev(hcl.colors(100, "YlGn")),
     addfun = function() lines(as_Spatial(conts)))
dim(test_capacity$carrying_capacity) # [1] 9242 321
plot(region$raster_from_values(test_capacity$carrying_capacity[, plot_seq]),
     col = rev(hcl.colors(700, "Terrain2")),
     zlim = c(0, 700),
     legend = TRUE,
     addfun = function() lines(as_Spatial(conts)))

plot(region$raster_from_values(test_capacity$carrying_capacity[, 321]),
     col = rev(hcl.colors(100, "Terrain2")),
     zlim = c(0, 100),
     legend = FALSE,
     addfun = function() lines(as_Spatial(conts)))

# Dispersal generator
## Need to run once to calculate distance based dispersal
## will take some time to run
## run once, save to RDS, then read in
dispersal_gen <- DispersalGenerator$new(region = region,
                                        dispersal_max_distance = 50, # km
                                        distance_classes = seq(10, 50, 5), # 10 to 50 km
                                        distance_scale = 1000, # km
                                        inputs = c("dispersal_p", ## dispersal_proportion
                                                   "dispersal_b"), ## dispersal_breadth
                                        decimals = 3)
dispersal_gen$calculate_distance_data(distance_matrix = distance_matrix) # pre-calculate
saveRDS(dispersal_gen$distance_data, "baseline/sim_inputs/numbat_dispersal_distance_data_projected.RDS")
#dispersal_gen$distance_data <- readRDS("baseline/sim_inputs/numbat_dispersal_distance_data_projected.RDS")
head(dispersal_gen$distance_data$base)
summary(dispersal_gen$distance_data$base) ## distance class max == 5 == 50 km

## test the dispersal generator
system.time({test_dispersal <- dispersal_gen$generate(input_values =
                                                        list(
                                                          # proportion of pop. disp
                                                          dispersal_p = 0.5,
                                                          # shape of curve (average disp. distance)
                                                          dispersal_b = 20
                                                        ))$dispersal_data})
sum(sapply(test_dispersal, anyNA)) # 0
lapply(test_dispersal, dim) # [1] 177444 5
summary(test_dispersal[[1]])


#### Step 6: Create generators for metapredator abundance ####

## metapred density generator
metapred_abundance_gen <- Generator$new(description = "metapred density generator",
                                        fox_mean = fox_abund_mean,
                                        fox_lower = fox_abund_lower,
                                        fox_upper = fox_abund_upper,
                                        cat_mean = cat_abund_mean,
                                        cat_lower = cat_abund_lower,
                                        cat_upper = cat_abund_upper,
                                        spatial_correlation = env_corr,
                                        generate_rasters = FALSE,
                                        generative_requirements = list(
                                          # distrib_var = "function",
                                          p_window = "function",
                                          meta_lower = "function",
                                          meta_upper = "function",
                                          meta_mean = "function",
                                          human_density = "distribution"),
                                        inputs = c("p", "cat_weight"),
                                        outputs = c("human_density"))

metapred_abundance_gen$add_function_template("p_window",
                                             function_def = function(params) {
                                               w <- params$p*10/100
                                               p_lower <- params$p - w
                                               p_upper <- params$p + w
                                               p_lower <- ifelse(p_lower < 0, 0, p_lower)
                                               p_upper <- ifelse(p_upper > 1, 1, p_upper)
                                               return(c(p_lower, p_upper))
                                             },
                                             call_params = c("p"))

metapred_abundance_gen$add_function_template("meta_lower",
                                             function_def = function(params) {
                                               meta_lower <- ((1-params$cat_weight)*params$fox_lower + params$cat_weight*params$cat_lower)/1
                                               return(meta_lower)
                                             },
                                             call_params = c("cat_weight", "fox_lower", "cat_lower"))
metapred_abundance_gen$add_function_template("meta_mean",
                                             function_def = function(params) {
                                               meta_mean <- ((1-params$cat_weight)*params$fox_mean + params$cat_weight*params$cat_mean)/1
                                               return(meta_mean)
                                             },
                                             call_params = c("cat_weight", "fox_mean", "cat_mean"))
metapred_abundance_gen$add_function_template("meta_upper",
                                             function_def = function(params) {
                                               meta_upper <- ((1-params$cat_weight)*params$fox_upper + params$cat_weight*params$cat_upper)/1
                                               return(meta_upper)
                                             },
                                             call_params = c("cat_weight", "fox_upper", "cat_upper"))

metapred_abundance_gen$add_distribution_template("human_density",
                                                 distr_type = "triangular",
                                                 distr_params = list(
                                                   lower = "meta_lower",
                                                   mode = "meta_mean",
                                                   upper = "meta_upper"),
                                                 sample = c("p_window"),
                                                 normalize_threshold = NULL)

# tests at 3 values (1, 0.01, 0.5) with cat_weight=0.2
{system.time({test_mp <- metapred_abundance_gen$generate(
  input_values = list(p = 0.1, cat_weight = 0.6))}) ## < 1 second
  region$raster_from_values(test_mp$human_density[, plot_seq])
  print(region$raster_from_values(test_mp$human_density[, plot_seq]))
  # quick plot to check
  plot(region$raster_from_values(test_mp$human_density[, plot_seq]),
       col = hcl.colors(100, "Zissou"),
       legend = TRUE,
       zlim = c(0, 1),
       addfun = function() lines(as_Spatial(conts)))
}

{system.time({test_mp <- metapred_abundance_gen$generate(
  input_values = list(p = 0.01, cat_weight = 0.2))}) ## < 1 second
  region$raster_from_values(test_mp$human_density[, plot_seq])
  print(region$raster_from_values(test_mp$human_density[, plot_seq]))
  # quick plot to check
  plot(region$raster_from_values(test_mp$human_density[, plot_seq]),
       col = hcl.colors(100, "Zissou"),
       legend = TRUE,
       zlim = c(0, 1),
       addfun = function() lines(as_Spatial(conts)))
}

{system.time({test_mp <- metapred_abundance_gen$generate(
  input_values = list(p = 0.5, cat_weight = 0))}) ## < 1 second
  region$raster_from_values(test_mp$human_density[, plot_seq])
  print(region$raster_from_values(test_mp$human_density[, plot_seq]))
  # quick plot to check
  plot(region$raster_from_values(test_mp$human_density[, plot_seq]),
       col = hcl.colors(100, "Zissou"),
       legend = TRUE,
       zlim = c(0, 1),
       addfun = function() lines(as_Spatial(conts)))
}

#### Step 7: Test simulation model ####
test_model <- model_template$clone()

test_model$set_attributes(carrying_capacity = test_capacity$carrying_capacity, # generated
                          initial_abundance = test_capacity$initial_abundance, # generated
                          growth_rate_max = log(1.30), # LHS
                          standard_deviation = 0.1, # LHS
                          density_max = 1380, # LHS e.g. 4.6/km2*(0.75*(20*20)) == 1380 per grid cell
                          dispersal_data = test_dispersal, # generated
                          dispersal_target_k = NULL, # NULL
                          occupancy_threshold = 10, # fixed or LHS. Suggest fixed. e.g. <= 10 cells == end simulation
                          abundance_threshold = 10, # LHS (allee effect) e.g. <= 10 individuals == kill cell
                          predation_max = 0.75, # maximum of 50% predated # LHS
                          predation_z = 1.25, # LHS
                          human_density = metapred_abundance_gen$generate(input_values = list(p = 0.5, cat_weight = 0.2))$human_density,  #we can also drop cat_weight and p and write human_density = test_mp$human_density 
                          results_selection = c("abundance", "predated"),
                          random_seed = 20210507)
test_model$list_completeness()
test_model$is_complete()
system.time({test_result <- paleopop_simulator(test_model)})
dim(test_result$abundance) #[1] 9242  321

dev.off()
op <- par()
par(mfrow = c(2, 1), mar = c(1,4,1,1))
plot(x = timeseq,
     y = colSums(test_result$abundance[, -c(1:burn_in_steps)], na.rm = FALSE),
     type = "l", xaxt = "n", xlab = "", ylab = "")
legend("topright", "Abundance", bty = "n", lty = 1, col = "black")
par(mar = c(2,4,0,1))
plot(x = timeseq,
     y = colSums(test_result$predated[,-c(1:burn_in_steps)], na.rm = FALSE),
     col = "red", type = "l", ylab = "", xlab = "")
legend("topright", "predation", bty = "n", lty = 1, col = "red")
par(op)
ab <- region$raster_from_values(test_result$abundance[, plot_seq])
abh <- region$raster_from_values(test_result$predated[, plot_seq])
max(cellStats(ab, max))
max(cellStats(abh, max))
plot(ab,
     # colNA = "#bf616a",
     legend = FALSE, nc = 3,
     col = hcl.colors(100, "Rocket"),
     zlim = c(0, max(cellStats(ab, max))),
     addfun = function() lines(as_Spatial(conts)))
plot(abh,
     # colNA = "#bf616a",
     legend = FALSE,
     col = hcl.colors(100, "Zissou"),
     zlim = c(0, max(cellStats(abh, max))),
     addfun = function() lines(as_Spatial(conts)))


# Extract human_density from the test_model
human_density_time <- test_model$human_density

# Sum human density across all grid cells for each timestep
human_density_sum <- colSums(human_density_time, na.rm = TRUE)

# Plot human_density over time
par(mfrow = c(1, 1), mar = c(4, 4, 2, 1))
plot(x = timeseq,
     y = human_density_sum[-c(1:burn_in_steps)], 
     type = "l",
     col = "blue",
     lwd = 2,
     xlab = "Year",
     ylab = "Human Density",
     main = "Human Density Over Time")


#### Step 8: Latin hypercube generation ####

# Generate the parameters for the latin-hypercube sampler
lhs_gen <- LatinHypercubeSampler$new()
lhs_gen$set_uniform_parameter("standard_deviation", lower = 0.00, upper = 0.35, decimals = 3)
lhs_gen$set_uniform_parameter("growth_rate_max", lower = log(1.20), upper = log(1.40), decimals = 3)
# Density max values
## Density: animals/km2 needs to be scaled by grid size (20km x 20km)
## e.g. an estimate of 7/km2 translates to 2800
## However, we assume only ~75% of the cell is suitable, so densities have
## to be scaled by 75%. E.g. 7/km2 == 7*(0.75*(20*20)) == 2,100 / grid cell
## range = 300 : 1500, corresponds to ~ 1/km2:5/km2
lhs_gen$set_uniform_parameter("density_max", lower = 300, upper = 1500, decimals = 0)
lhs_gen$set_uniform_parameter("dispersal_p", lower = 0.001, upper = 0.600, decimals = 3)
lhs_gen$set_uniform_parameter("dispersal_b", lower = 0, upper = 25, decimals = 0)
lhs_gen$set_uniform_parameter("abundance_threshold", lower = 15, upper = 60, decimals = 0) # allee effect
lhs_gen$set_uniform_parameter("predation_max", lower = 0.05, upper = 0.85, decimals = 3)
lhs_gen$set_uniform_parameter("predation_z", lower = 1, upper = 2)
lhs_gen$set_uniform_parameter("p", lower = 0.1, upper = 0.7, decimals = 2) # sample range for cats
sample_data <- lhs_gen$generate_samples(number = nsims, random_seed = 42)
head(sample_data)
dim(sample_data)

# Sample the HS projections
sample_rows <- {set.seed(17559326); setDT(niche_lookup)[sample(.N, nsims, replace = TRUE), ][, -2]}
sample_rows <- as.data.frame(sample_rows)
head(sample_rows)

round(prop.table(table(as.numeric(sapply(strsplit(split = "_", x = niche_lookup$niche_ref), "[", 1)))), 3)
round(prop.table(table(as.numeric(sapply(strsplit(split = "_", x = sample_rows$niche_ref), "[", 1)))), 3)

# bind the HS samples, back to the LHS demographic parameters
sample_data <- as.data.frame(cbind(sample_rows, sample_data))
head(sample_data)

plot_samples <- TRUE
if (plot_samples) {
   par(mfrow = n2mfrow(ncol(sample_data)-1))
   for (j in 2:ncol(sample_data)) {
     d <- sample_data[, j]
     plot(density(d, kernel = "gaussian",
                  from = min(d), to = max(d)),
          main = "", ylab = "",
          xlim = range(d),
          xlab = colnames(sample_data)[j])
 }
}
par(mfrow = c(1, 1))

# Make unique row names for saving files
{set.seed(54612)
  sample_data$UniqueID <- paste0(stringi::stri_rand_strings(nsims, 4, "[A-Z]"),
                                 stringi::stri_rand_strings(nsims, 4, "[0-9]"),
                                 "_", sample_data$niche_ref)
}
any(duplicated(sample_data$UniqueID))
head(sample_data$UniqueID)
sample_data <- sample_data[, c(12, 1:11)]
head(sample_data)

file_check <- pbsapply(unique(sample_data$niche_ref), function(x) {
 file.exists(file.path(data_dir, sprintf("%s.RDS", x)))
})
sum(file_check) == length(unique(sample_data$niche_ref))

# write the LHS samples to file
fwrite(sample_data, "baseline/lhs_25000.csv")
#sample_data <- fread("baseline/lhs_25000.csv"")

#### Step 9: Run the simulations ####
# Create a simulation manager and run the sampled model simulations
sim_model <- model_template$clone()
sim_model$random_seed <- NULL
sim_model$set_attributes(occupancy_threshold = 10, # end sim if 10 cells occupied
                         # save abundance, total predation, and extirpation time
                         results_selection = c("abundance", "predated", "extirpation"))
sim_model$attribute_aliases
sim_model$is_complete() # FALSE. Other params set by LHS

test_sims <- FALSE # if TRUE will run 100 LHS simulations

if (test_sims) {
  if (!dir.exists(test_dir)) {
    dir.create(test_dir, recursive = TRUE)
  }
  ## generate a new 100 sample LHS
  lhs_gen <- LatinHypercubeSampler$new()
  lhs_gen$set_uniform_parameter("standard_deviation", lower = 0.00, upper = 0.35, decimals = 3)
  lhs_gen$set_uniform_parameter("growth_rate_max", lower = log(1.20), upper = log(1.40), decimals = 3)
  lhs_gen$set_uniform_parameter("density_max", lower = 300, upper = 1500, decimals = 0)
  lhs_gen$set_uniform_parameter("dispersal_p", lower = 0.001, upper = 0.600, decimals = 3)
  lhs_gen$set_uniform_parameter("dispersal_b", lower = 0, upper = 25, decimals = 0)
  lhs_gen$set_uniform_parameter("abundance_threshold", lower = 15, upper = 60, decimals = 0)
  lhs_gen$set_uniform_parameter("predation_max", lower = 0.001, upper = 0.60, decimals = 3)
  lhs_gen$set_uniform_parameter("predation_z", lower = 1, upper = 2)
  lhs_gen$set_uniform_parameter("p", lower = 0, upper = 1, decimals = 2)
  lhs_gen$set_uniform_parameter("cat_weight", lower = 0, upper = 1, decimals = 2)
  test_data <- lhs_gen$generate_samples(number = 100, random_seed = 42)
  # Sample the HS projections
  sample_rows <- {set.seed(943); setDT(niche_lookup)[sample(.N, 100, replace = TRUE), ][, -2]}
  sample_rows <- as.data.frame(sample_rows)
  # bind the HS samples, back to the LHS demographic parameters
  test_data <- as.data.frame(cbind(sample_rows, test_data))
  head(test_data)
  plot_samples <- TRUE
  if (plot_samples) {
    par(mfrow = n2mfrow(ncol(test_data)-1))
    for (j in 2:ncol(test_data)) {
      d <- test_data[, j]
      plot(density(d, kernel = "gaussian", from = min(d), to = max(d)),
           main = "", ylab = "",
           xlim = range(d),
           xlab = colnames(test_data)[j])
    }
  }
  par(mfrow = c(1,1))
  # Make unique row names for saving test simulation files
  {set.seed(54612); test_data$UniqueID <- paste0(stringi::stri_rand_strings(100, 3, "[A-Z]"),
                                                   stringi::stri_rand_strings(100, 3, "[0-9]"),
                                                   "_", test_data$niche_ref)}
  stopifnot(!any(duplicated(test_data$UniqueID)))
  test_data <- test_data[, c(12, 1:11)]
  file_check <- pbsapply(unique(test_data$niche_ref), function(x) {
    file.exists(file.path(data_dir, sprintf("%s.RDS", x)))
  })
  stopifnot(sum(file_check) == length(unique(test_data$niche_ref)))
  head(test_data)
}

if (test_sims) {
  message(sprintf("RUNNING %s TEST SIMULATIONS!!!", nrow(test_data)))
  sim_manager <- SimulationManager$new(sample_data = test_data,
                                       model_template = sim_model,
                                       generators = list(capacity_gen,
                                                         dispersal_gen,
                                                         metapred_abundance_gen),
                                       parallel_cores = parallel_cores,
                                       results_filename_attributes =
                                         c(NULL, "UniqueID", "results"),
                                       results_dir = test_dir)
  sim_manager$model_template$list_completeness()
  run_output <- sim_manager$run()
  run_output$summary
  } else {
  message(sprintf("RUNNING ALL %.0f SIMULATIONS!!!", nsims))
    
  # Start timing
  start_time <- proc.time()  
  ## run all simulations!
  sim_manager <- SimulationManager$new(sample_data = sample_data,
                                       model_template = sim_model,
                                       generators = list(capacity_gen,
                                                         dispersal_gen,
                                                         metapred_abundance_gen),
                                       parallel_cores = parallel_cores,
                                       results_filename_attributes =
                                         c(NULL, "UniqueID", "results"),
                                       results_dir = out_dir)
  run_output <- sim_manager$run()
  run_output$summary
  # End timing and calculate duration
  end_time <- proc.time()
  run_duration <- end_time - start_time
  cat("Total time taken for simulations (in seconds):\n")
  print(run_duration)
}

out_files <- list.files(out_dir, ".RData", recursive = FALSE)
length(out_files)

if (length(out_files) < nsims) {
  message(sprintf("There are %s simulations to re-run", nsims - length(out_files)))
  out_sims <- gsub("UniqueID_|_results.RData", "", out_files)
  re_runs <- which(!sample_data$UniqueID %in% out_sims)
  re_runs <- sample_data[re_runs, ]
  stopifnot(nrow(re_runs) == nsims - length(out_files))
  # run the extra models!
  message(sprintf("RUNNING ALL %.0f REMAINING SIMULATIONS!!!", nrow(re_runs)))
  sim_manager <- SimulationManager$new(sample_data = re_runs,
                                       model_template = sim_model,
                                       generators = list(capacity_gen,
                                                         dispersal_gen,
                                                         metapred_abundance_gen),
                                       parallel_cores = parallel_cores,
                                       results_filename_attributes =
                                         c(NULL, "UniqueID", "results"),
                                       results_dir = out_dir)
  run_output <- sim_manager$run()
  run_output$summary
}

# quick summaries of results
results <- PaleoPopResults$new(burn_in_steps = burn_in_steps,
                               region = region,
                               time_steps = (timesteps - burn_in_steps))
res_manager <- ResultsManager$new(simulation_manager = sim_manager,
                                  simulation_results = results,
                                  generators = NULL,
                                  # total abundance & total predation
                                  summary_metrics = c("total_n",
                                                      "total_h"),
                                  summary_matrices = c("n", "h"),
                                  summary_functions = list(
                                    total_n = function(results) {
                                      sum(results$abundance)
                                    },
                                    total_h = function(results) {
                                      sum(results$predated)
                                    },
                                    n = "all$abundance",
                                    h = "all$predated"),
                                  parallel_cores = parallel_cores)

# Generate the results outputs
gen_log <- res_manager$generate()
gen_log$summary
gen_log$full_log[[1]]
lapply(res_manager$summary_matrix_list, dim)

res_ext <- data.frame(Abund = res_manager$summary_matrix_list$n[1, ],
                      Harv = res_manager$summary_matrix_list$h[1, ],
                      Time = attr(temp_mask, "timesteps"))
head(res_ext)

{with(res_ext, plot(x = Time, y = Abund, type = "l", xlab = "Year",
                   ylab = "Total Abundance", ylim = c(0, max(Abund))))
with(res_ext, lines(x = Time, y = Harv, col = "red"))}

options(scipen = 999)  # Turn off scientific notation globally
dev.off()
 op <- par()
 png(filename = "baseline/simulation_plots_final.png", height = 5, width = 8, res = 320, units = "in", bg = "white")
par(mfrow = c(2, 1), mar = c(1,4,1,1))
matplot(t(res_manager$summary_matrix_list$n), type = "l",
        col = "#B8B8B86D", lty = 1, xaxt = "n",
        x = timeseq, ylab = "Total abundance")
lines(apply(res_manager$summary_matrix_list$n, 2, mean, na.rm = TRUE), 
      col = "#F513AAC7", x = timeseq)
lines(apply(res_manager$summary_matrix_list$n, 2, median, na.rm = TRUE), 
      col = "#F513AAC7", lty = 2, x = timeseq)
legend("topright", "Abundance", bty = "n", lty = 1, col = "black")
par(mar = c(2,4,0,1))
matplot(t(res_manager$summary_matrix_list$h), type = "l",
        col = "#B8B8B86D", lty = 1, x = timeseq, xlab = "Year",
        ylab = "Total predation")
lines(apply(res_manager$summary_matrix_list$h, 2, mean, na.rm = TRUE), 
      col = "#F513AAC7", x = timeseq)
lines(apply(res_manager$summary_matrix_list$h, 2, median, na.rm = TRUE), 
      col = "#F513AAC7", lty = 2, x = timeseq)
legend("topright", "predation", bty = "n", lty = 1, col = "black")
 dev.off()
 par(op)

