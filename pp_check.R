pp_check <- function(model_metrics, obs_targets, distance_function = "euclidean",
                     test = "t.test") {
  if (missing(obs_targets))
    stop("'obs_targets' is missing")
  if (missing(model_metrics))
    stop("'model_metrics' is missing")
  if (!is.list(obs_targets) && !is.vector(obs_targets))
    stop("'obs_targets' has to be a list or vector.")
  if (!is.matrix(model_metrics) && !is.data.frame(model_metrics))
    stop("'model_metrics' has to be a matrix or data.frame.")
  if (is.data.frame(obs_targets))
    obs_targets <- as.matrix(obs_targets)
  if (is.data.frame(model_metrics))
    model_metrics <- as.matrix(model_metrics)
  if (is.list(obs_targets))
    obs_targets <- unlist(obs_targets)
  if (!(distance_function %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")))
    stop("'distance_function' must be one of the options available in dist()")
  if (!(test %in% c("t.test", "wilcox.test")))
    stop("You must choose either 't.test' (for large samples) or 'wilcox.test' (for smaller samples)")
  if (length(obs_targets) != dim(model_metrics)[2])
    stop("Number of summary statistics in 'obs_targets' has to be the same as in 'model_metrics'.")
  m <- nrow(model_metrics)
  if (m >= 100000)
    stop("Too many rows in 'model_metrics'; your computer will very likely crash.")
  if (m > 10000)
    warning("Large number of rows in 'model_metrics'; your computer may crash from computational strain.")
  nss <- length(model_metrics[1,])
  cond1 <- !any(as.logical(apply(model_metrics, 2, function(x)
    length(unique(x)) - 1)))
  if (cond1)
    stop("Zero variance in the summary statistics.")
  gwt <- rep(TRUE, length(model_metrics[, 1]))
  gwt[attributes(na.omit(model_metrics))$na.action] <- FALSE
  if (!length(colnames(model_metrics))) {
    warning("No summary statistics names are given, using S1, S2, ...")
    statnames <- paste("S", 1:nss, sep = "")
  } else {
    statnames <- colnames(model_metrics)
  }
  scaled.sumstat <- model_metrics
  sdScale <- function (x, y) {
    if (sd(y) == 0) {
      return(x)
    } else {
      return(x / sd(y))
    }
  }
  for (j in 1:nss) {
    scaled.sumstat[, j] <- sdScale(x = model_metrics[, j], y = model_metrics[, j][gwt])
  }
  target <- obs_targets
  for (j in 1:nss) {
    target[j] <- sdScale(obs_targets[j], model_metrics[, j][gwt])
  }
  d_obs <- rep(NA, m)
  for (j in 1:m) {
    d_obs[j] <- dist(
          matrix(c(scaled.sumstat[j,], target), ncol = dim(scaled.sumstat)[2], byrow = T),
          method = distance_function)
  }
  d_rep <- as.vector(dist(scaled.sumstat, method = distance_function))
  if (test == "t.test")
    output <- t.test(d_obs, d_rep, alternative = "greater")
  if (test == "wilcox.test")
    output <- wilcox.test(d_obs, d_rep, alternative = "greater")
  return(output)
}
