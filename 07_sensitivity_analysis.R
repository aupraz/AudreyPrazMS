# Load stacks
wtdAbund      <- stack("third_try/final_selection/wtdAbund_threshfinal.tif")
abund_cat     <- stack("counterfactuals/abc/Abund_cat.tif")
abund_foxes   <- stack("counterfactuals/abc/Abund_fox.tif")
abund_nolanduse <- stack("counterfactuals/abc/Abund_landuse.tif")
abund_nopred  <- stack("counterfactuals/abc/Abund_metapred.tif")
abund_off     <- stack("counterfactuals/abc/Abund_off.tif")
abund_foxnoland <- stack("counterfactuals/abc/Abund_landuse_fox.tif")
abund_catnoland <- stack("counterfactuals/abc/Abund_landuse_cat.tif")


stacks <- list(
  Baseline                = wtdAbund,
  No_landuse              = abund_nolanduse,
  Without_foxes           = abund_cat,
  Without_cats            = abund_foxes,
  No_landuse_and_foxes    = abund_catnoland,
  No_landuse_and_cats     = abund_foxnoland,
  Without_foxes_and_cats  = abund_nopred,        
  All_threats_removed     = abund_off
)

years <- 1750:(1750 + nlayers(wtdAbund) - 1)

# Define a function to get the total abundance for each layer
get_ts <- function(x) sapply(1:nlayers(x), function(i) sum(values(x[[i]]), na.rm=TRUE))


ts_list <- lapply(stacks, get_ts)
ts_dt <- as.data.table(ts_list)
ts_dt[, Year := years]  
setcolorder(ts_dt, c("Year", names(ts_list)))
abund_ts <- ts_dt


library(stats)

# define factor levels
fact_df <- data.table(
  land = c(1,0,1,1,0,0,1,0),
  fox  = c(1,1,0,1,0,1,0,0),
  cat  = c(1,1,1,0,1,0,0,0)
)

# Prepare storage for Sobol' main-effect indices
S_land <- S_fox <- S_cat <- numeric(length(years))

# Loop over years, fit full-factorial ANOVA, extract SS ratios
for(i in seq_along(years)) {
  Yt <- as.numeric(abund_ts[i, -1, with=FALSE])   # eight values
  df <- cbind(fact_df, Y = Yt)
  
  # full-factorial model
  mod <- lm(Y ~ land*fox*cat, data = df)
  a   <- anova(mod)
  SS  <- setNames(a$`Sum Sq`, rownames(a))
  
  totalSS <- sum(SS)
  S_land[i] <- SS["land"]
  S_fox[i]  <- SS["fox"]
  S_cat[i]  <- SS["cat"]
  
  # convert to proportions
  S_land[i] <- S_land[i] / totalSS
  S_fox[i]  <- S_fox[i]  / totalSS
  S_cat[i]  <- S_cat[i]  / totalSS
}

# Create a df for the results
sobol_df <- data.table(
  Year    = years,
  LandUse = S_land,
  Fox     = S_fox,
  Cat     = S_cat
)

# Plot 
m <- melt(sobol_df, id.vars="Year", variable.name="Threat", value.name="SobolIndex")

ggplot(m, aes(x=Year, y=SobolIndex, color=Threat)) +
  geom_line(size=1) +
  labs(
    title = "Time-resolved Sobol' First-Order Indices",
    y     = "Proportion of Variance Explained",
    x     = "Year"
  ) +
  theme_minimal()

write.csv(sobol_df, "counterfactuals/sobol_df.csv")
