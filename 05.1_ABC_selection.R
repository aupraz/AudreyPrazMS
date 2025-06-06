library(data.table)
library(abc)
library(raster)
library(poems)
library(paleopop)
library(parallel)
library(pbapply)
library(pbarETA, warn = FALSE)

# number of cores for parallel computation
parallel_cores <- 10L

# source updated ABC and extra plotting functions
source("code/helper_funs.R")

# simulation parameters

sim_par <- fread("baseline/sim_inputs/lhs_25000.csv")

#sim_par[, UID := 1:.N]
#setcolorder(sim_par, c(12, 1:11))

# Simulated targets
sim_tar <- setDT(readRDS("baseline/validation_targets/validation_targets_summarised.RDS"))

sim_tar[, ExtPatternP := round(ExtPattern/5, 1)]
sim_tar[, Remnant := Persistence_locns]

# Ensure columns are characters to avoid mismatch due to factors
sim_par[, UniqueID := as.character(UniqueID)]
sim_tar[, SampleID := as.character(SampleID)]

# Filter sim_par to only keep rows whose UniqueID is in sim_tar
sim_par <- sim_par[UniqueID %in% sim_tar$SampleID]

# Optionally, also sort both tables by ID to ensure alignment
setkey(sim_par, UniqueID)
setkey(sim_tar, SampleID)

na_idx <- apply(sim_tar, 1, anyNA)
sum(na_idx) # Some rows have NA for extirpation targets 
## NA either means no extirpation or extirpation during burn-in period
valid_rows <- !na_idx & complete.cases(sim_par[, -c(1:2)])

# Actual targets
act_tar <- data.table(
  #AllSites_1 = 1.0, 
  SitesProp_9 = 1.0,
  Persistence_Target = 271,
  Persistence_dist = 0,
  ExtPen = 0
  # Final_abund = 5000
  # ExtPatternP = 1,
  # ExtPen.1886 = 0,
  # ExtPen.1938 = 0,
  # ExtPen.2003 = 0,
  # ExtPen.2009 = 0,
  # #ExtPen.2008 = 0,
  # ExtPattern = 5,
  #ExtPatternP = 1,
  #ExtAbsPen = 0
  # Remnant = 1
)
apply(sim_tar, 2, mad, na.rm = TRUE) # 'Remnant' target will fail. MAD of 0.
apply(sim_tar, 2, sd, na.rm = TRUE) # use SD approach instead
act_tar_replicated <- act_tar[rep(1, nrow(sim_tar)), ]  # Replicate act_tar rows to match the row of sim_tar

# replace normalise function in abc
## following https://doi.org/10.1016/j.ecolmodel.2015.05.020
assignInNamespace("normalise", normalise, ns = "abc")
(abc:::normalise) # now using sd rather than MAD

# ABC

tars <- c("SitesProp_9","Persistence_Target", "Persistence_dist","ExtPen")      
          
first_pass <- abc(target = as.matrix(act_tar[, ..tars]),
                  param = as.matrix(sim_par[, -c(1:2)]),
                  sumstat = as.matrix(sim_tar[, ..tars]),
                  subset = !na_idx, # ignore sims with NA for targets
                  tol = 0.005, # best 100 sims 
                  #tol = 100/nrow(sim_par), # best 100 sims 
                  method = "rejection")
first_pass
summary(first_pass$ss)
summary(first_pass)

summary(sim_par[first_pass$region == TRUE, ]) #this won't run properly until all the simulations in the sim_par file have been run out

summary(sim_tar[first_pass$region == TRUE, ])

# plot
plot.rej(first_pass, sim_par[, -c(1:2)], "gaussian")

plot.targets(x = first_pass, 
             model_targets = sim_tar[!na_idx, ..tars], 
             obs_targets = unlist(act_tar[,..tars]))+
  scale_y_continuous(trans = "log1p", breaks = c(0, 100, 500, 2500, 5000, 10000), fill=TRUE)

saveRDS(first_pass, "baseline/abc/abc_first_pass.RDS", compress = TRUE)

# make a plot of weighted abundances etc
## weights are defined as euc.dist from idealised targets.
# Quick ensemble raster from all sims
out_dir <- "baseline/sims_25000"
files <- sim_par[first_pass$region == TRUE, ][["UniqueID"]]
files <- paste0("UniqueID_", files, "_results.RData")
files <- file.path(out_dir, files)
stopifnot(sum(sapply(files, function(x) file.exists(x))) == sum(first_pass$region))

burn_in_steps <- 50 # 25 generations
timesteps <- 271 + burn_in_steps
plot_seq <- floor(seq(burn_in_steps+1, timesteps, length.out = 6))
region <- readRDS("baseline/sim_inputs/numbat_paleopopregion.RDS")

# abundances
ras_list <- lapply(seq_len(timesteps)[-c(1:burn_in_steps)], function(time) {
  cls <- makeCluster(parallel_cores)
  clusterExport(cls, c("files", "time"))
  tempT <- pbsapply(1:length(files), function(x) {
    readRDS(files[x])$abundance[, time]
  }, cl = cls)
  stopCluster(cls)
  mat <- matrix(nrow = nrow(tempT), ncol = 1)
  mat[, 1] <- apply(tempT, 1, function(x) {floor(weighted.mean(x = x, w = 1/first_pass$dist[first_pass$region == TRUE], na.rm = TRUE))})
  rs <- region$raster_from_values(mat)
  return(rs)
})
ras_list[plot_seq-burn_in_steps]
wtdAbund <- stack(ras_list, quick = TRUE)
wtdAbund[[plot_seq-50]]
plot(wtdAbund[[plot_seq-50]], zlim = c(0, 800),
     col = hcl.colors(100, "Zissou"),
     legend = FALSE)

writeRaster(wtdAbund, "baseline/abc/wtdAbund_25000.grd")
#wtdAbund <- stack("baseline/abc/wtdAbund_25000.grd")

# predation rates
ras_list <- lapply(seq_len(timesteps)[-c(1:burn_in_steps)], function(time) {
  cls <- makeCluster(50L)
  clusterExport(cls, c("files", "time"))
  tempT <- pbsapply(1:length(files), function(x) {
    readRDS(files[x])$predated[, time]
  }, cl = cls)
  stopCluster(cls)
  mat <- matrix(nrow = nrow(tempT), ncol = 1)
  mat[, 1] <- apply(tempT, 1, function(x) {floor(weighted.mean(x = x, w = 1/first_pass$dist[first_pass$region == TRUE], na.rm = TRUE))})
  rs <- region$raster_from_values(mat)
  return(rs)
})
ras_list[plot_seq-burn_in_steps]
wtdHarv <- stack(ras_list, quick = TRUE)
wtdHarv[[plot_seq-50]]
plot(wtdHarv[[plot_seq-50]], zlim = c(0, 60),
     col = hcl.colors(100, "Zissou"),
     legend = FALSE)

#writeRaster(wtdHarv, "baseline/abc/wtdHarv_25000.grd")
#wtdHarv<- stack("baseline/abc/wtdHarv_25000.grd")


# Try to create the plot to assess the abc selection 

# Bayesian analyses
feral.gfit=abc::gfit(target=unlist(act_tar[,..tars]), sumstat=first_pass$ss, 10, statistic=mean)
summary(feral.gfit)

My_first_BF <- pp_check(sim_tar[, ..tars], unlist(act_tar[,..tars]), distance_function = "euclidean",
                        test = "t.test")
# Extract t-test results
My_first_BF_results <- data.table(
  statistic = My_first_BF$statistic,
  p_value = My_first_BF$p.value
)

# Save as CSV
write.csv(My_first_BF_results, "baseline/abc/bayes_factors_run25000.csv", row.names = FALSE)

summary_metric_data <- readRDS("baseline/validation_targets/validation_targets_summarised_new.RDS")
summary_metric_data <- cbind(sim_par, summary_metric_data)

summary_metric_data_2 <- data.table(
  Model_ID = which(first_pass$region == TRUE),  # Index of accepted models
  Distance = first_pass$dist[first_pass$region == TRUE]  # Distances from target
)

summary_metric_data$ABC_distance <- first_pass$dist  # Distances from target

summary_metric_data[, Weights := 1/(ABC_distance + .Machine$double.eps)]
print(head(summary_metric_data))
write.csv(summary_metric_data, "baseline/validation_targets/validation_targets_summarised.RDS", row.names = FALSE)



# Filter sim_par to keep only the selected models
selected_sims <- summary_metric_data[first_pass$region == TRUE, ]

selected_sims[, `:=`(
  Distance = first_pass$dist[first_pass$region == TRUE],  
  Weights = 1 / (ABC_distance + .Machine$double.eps))]

head(selected_sims)

write.csv(selected_sims, "baseline/abc/selected_sims.csv", row.names = FALSE)

