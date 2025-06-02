# function for subsetting/performing the cuts
subset_data <- function(data, position_fraction, width_fraction) {
  if (position_fraction - width_fraction/2 < 0 || position_fraction + width_fraction/2 > 1) {
    stop("Bounds exceed range of data")
  }
  if (width_fraction > 1 || width_fraction <= 0) {
    stop("Width fraction must be in (0,1]")
  }
  if (width_fraction > 1 || width_fraction < 0) {
    stop("Position fraction must be in [0,1]")
  }
  data_scaled <- data
  for (i in 1:ncol(data)) {
    data_scaled[,i] <- (data_scaled[,i] - min(data_scaled[,i],na.rm=T)) / (max(data_scaled[,i],na.rm=T) - min(data_scaled[,i],na.rm=T))#as.numeric(cut(indexdf[,i], breaks=ncuts, ordered=T))
  }
  data_scaled_index <- (data_scaled >= (position_fraction - width_fraction/2) & data_scaled <= (position_fraction + width_fraction/2))
  data_scaled_index  <- apply(data_scaled_index, 1, function(x) { all(x>0) })	
  return(data_scaled_index)
}

# function for normalising by SD rather than mad
normalise <- function(x, y) {
  if (sd(y) == 0) {
    return(x)
  } else {
    return(x/sd(y))
  }
}

# Function for plotting ABC results with the rejection method posteriors
plot.rej <- function(x, param, dens_kernel = "epanechnikov", minq = 0.025, maxq = 0.975, ...) {
  opar <- par()
  if (!inherits(x, "abc"))
    stop("Use only with objects of class \"abc\".", call. = FALSE)
  abc.out <- x
  mymethod <- abc.out$method
  if (mymethod != "rejection") {
    stop("plot.rej can only be used when method is 'rejection'",
         call. = FALSE)
  }
  if (!is.matrix(param) && !is.data.frame(param) && !is.vector(param)) {
    stop("'param' has to be a matrix, data.frame or vector.",
         call. = F)
  }
  if (is.null(dim(param))) {
    param <- matrix(param, ncol = 1)
  }
  if (is.data.frame(param)) {
    param <- as.matrix(param)
  }
  np <- abc.out$numparam
  numsim <- length(param)/np
  parnames <- abc.out$names$parameter.names
  cond <- isTRUE(c(match(colnames(param), parnames),
                   match(parnames, colnames(param))))
  if (cond) {
    stop("'abc.out' and 'param' are not compatible; paramater names are different.",
         call. = FALSE)
  }
  # selected with rejection or not?
  rej <- abc.out$unadj.values
  # make sure matrices are in same order
  if (is.vector(param)) {
    np.orig <- 1
    nsim <- length(param)
  } else if (is.matrix(param)) {
    np.orig <- dim(param)[2]
    nsim <- dim(param)[1]
    myorder <- match(parnames, colnames(param))
    if (isTRUE(myorder - 1:np)) {
      param <- param[, myorder]
      warning("'param' is being re-ordered according to 'abc.out'...",
              call. = FALL, immediate. = TRUE)
    }
  }
  if (np.orig != np) {
    stop("The number parameters supplied in \"param\" is different from that in \"x\".",
         call. = FALSE)
  }
  par(mfrow = n2mfrow(nr.plots = np))
  for (i in 1:np) {
    prior.d <- density(param[, i],
                       from = min(param[,i]),
                       to = max(param[, i]),
                       kernel = dens_kernel,...)
    rej.d <- density(rej[, i],
                     from = min(rej[, i], na.rm = T),
                     to = max(rej[, i], na.rm = TRUE),
                     kernel = dens_kernel,...)
    rej.dl <- density(rej[, i],
                      from = quantile(rej[, i], minq, na.rm = TRUE),
                      to = quantile(rej[, i], maxq, na.rm = TRUE),
                      kernel = dens_kernel, ...)
    myxlim <- range(c(prior.d$x, rej.d$x))
    myylim <- range(c(prior.d$y, rej.d$y))
    myxlab <- parnames[i]
    plot(prior.d, main = myxlab, col = "#1b9e77", lty = 1,
         lwd = 1, xlab = "", ylim = myylim, xlim = myxlim)
    title(sub = paste0(paste("N =", prior.d$n, "  Bandwidth =", formatC(prior.d$bw)), "; ",
                       sprintf("mean = %0.3f; sd = %0.3f",
                               mean(param[, i],na.rm = T), sd(param[, i],na.rm = T)), "\n",
                       sprintf("meanSel = %0.3f; sdSel = %0.3f",
                               mean(rej[, i],na.rm = T), sd(rej[, i],na.rm = T))))
    lines(rej.d, col = "#d95f02", lty = 1, lwd = 1)
    abline(v = quantile(rej[, i], 0.5), col = "#d95f02", lty = 2, lwd = 1)
    abline(v = quantile(param[, i], 0.5), col = "#1b9e77", lty = 2, lwd = 1)
    rug(rej[, i], side = 1, ticksize = 0.02, lwd = 1, col = "#d95f02")
    abline(v = range(rej.dl$x), col = "grey90", lty = 1, lwd = 1)
    if (i == np) {
      legend("bottomright",
             c("Prior", "Rej."),
             lty = c(1,1),
             col = c("#1b9e77", "#d95f02"),
             inset = c(0,1), xpd = TRUE, horiz = TRUE,
             bty="n"
      )
    }
  }
  par(opar)
  invisible()
}

# plot regression adjusted methods with rejection sims
plot.adj <- function(x, param) {
  if (!inherits(x, "abc"))
    stop("Use only with objects of class \"abc\".", call. = FALSE)
  abc.out <- x
  mymethod <- abc.out$method
  if (mymethod == "rejection")
    stop("plot.rej can only be used when method is \"loclinear\", \"neuralnet\" or \"ridge\".",
         call. = FALSE)
  if (!is.matrix(param) && !is.data.frame(param) && !is.vector(param))
    stop("'param' has to be a matrix, data.frame or vector.",
         call. = F)
  if (is.null(dim(param)))
    param <- matrix(param, ncol = 1)
  if (is.data.frame(param))
    param <- as.matrix(param)
  np <- abc.out$numparam
  numsim <- length(param)/np
  parnames <- abc.out$names$parameter.names
  cond <- isTRUE(c(match(colnames(param), parnames), match(parnames,
                                                           colnames(param))))
  if (cond)
    stop("'abc.out' and 'param' are not compatible; paramater names are different.",
         call. = FALSE)
  rej <- abc.out$unadj.values
  post <- abc.out$adj.values
  if (np == 1)
    rej <- matrix(rej, ncol = 1)
  if (is.vector(param)) {
    np.orig <- 1
    nsim <- length(param)
  } else if (is.matrix(param)) {
    np.orig <- dim(param)[2]
    nsim <- dim(param)[1]
    myorder <- match(parnames, colnames(param))
    if (isTRUE(myorder - 1:np)) {
      param <- param[, myorder]
      warning("'param' is being re-ordered according to 'abc.out'...",
              call. = FALL, immediate. = TRUE)
    }
  }
  if (np.orig != np) {
    stop("The number parameters supplied in \"param\" is different from that in \"x\".",
         call. = FALSE)
  }
  par(mfrow = n2mfrow(nr.plots = np))
  for (i in 1:np) {
    prior.d <- density(param[, i], from = min(param[,i]),
                       to = max(param[, i]))
    rej.d <- density(rej[, i], from = min(rej[, i], na.rm = T),
                     to = max(rej[, i], na.rm = TRUE))
    post.d <- density(post[, i], from = min(post[, i], na.rm = T),
                      to = max(post[, i], na.rm = TRUE))
    myxlim <- range(c(prior.d$x, rej.d$x, post.d$x))
    myylim <- range(c(prior.d$y, rej.d$y, post.d$y))
    myxlab <- parnames[i]
    plot(prior.d, main = myxlab, col = "#1b9e77", lty = 2,
         lwd = 1, xlab = "", ylim = myylim, xlim = myxlim)
    title(sub = paste0(paste("N =", prior.d$n, "  Bandwidth =", formatC(prior.d$bw)), "; ",
                       sprintf("mean = %0.3f; sd = %0.3f",
                               mean(param[, i],na.rm = T), sd(param[, i],na.rm = T)), "\n",
                       sprintf("meanRej = %0.3f; sdRej = %0.3f",
                               mean(rej[, i],na.rm = T), sd(rej[, i],na.rm = T)), "\n",
                       sprintf("meanAdj = %0.3f; sdAdj = %0.3f",
                               mean(post[, i],na.rm = T), sd(post[, i],na.rm = T))))
    lines(rej.d, col = "#d95f02", lty = 1, lwd = 1)
    lines(post.d, col = "#7570b3", lty = 1, lwd = 1)
    rug(rej[, i], side = 3, ticksize = 0.02, lwd = 1,
        col = "#d95f02")
    rug(post[, i], side = 1, ticksize = 0.02, lwd = 1,
        col = "#7570b3")
    if (i == np) {
      legend("bottomright",
             c("Prior", "Rej.", "Adj."),
             lty=c(1,1,1),
             col = c("#1b9e77", "#d95f02", "#7570b3"),
             inset = c(0,1), xpd = TRUE, horiz = TRUE,
             bty="n"
      )
    }
  }
  par(mfcol = c(1, 1))
}

# plot function for abc cross-validation
plot.abc.cv <- function(x, exclude = NULL, ...) {
  invisible(opar <- par())
  require(colorspace)
  if (!inherits(x, "cv4abc")){
    stop("Use only with objects of class \"cv4abc\".", call. = F)}
  cv4abc.out <- x
  tols <- cv4abc.out$tols
  numtols <- length(tols)
  np <- length(cv4abc.out$names$parameter.names)
  true <- cv4abc.out$true
  estim <- cv4abc.out$estim
  parnames <- cv4abc.out$names$parameter.names
  cv4abc.out$estim <- as.data.frame(cv4abc.out$estim)
  nval <- length(true)/np
  caption <- as.graphicsAnnot(parnames)
  cols <- hcl.colors(n = numtols, palette = "Batlow")
  cols <- colorspace::adjust_transparency(cols, alpha = 1/length(cols))
  transparent.colors = scales::alpha(cols, alpha = 1/length(cols))
  pch <- 19
  figma <- matrix(nrow = n2mfrow(np)[1],
                  data = 1:np,
                  byrow = TRUE)
  n1 <- layout(figma, widths = 1, heights = 1)
  par(mar = c(3,3,3,3))
  for (par in 1:np) {
    columns <- seq(par, numtols * np, by = np)
    if (!is.null(exclude)) {
      plot(rep(cv4abc.out$true[-exclude, par],
               times = numtols),
           unlist(cv4abc.out$estim[-exclude, columns]),
           xlab = NA, ylab = NA,
           col = rep(cols, each = nval), pch = pch)
      mtext("True value", side = 1, line = 1.75, cex = 0.75)
      mtext("Est. value", side = 2, line = 1.75, cex = 0.75)
      mtext(parnames[par], side = 3, line = 0, cex = 0.75, font = 2)
      abline(0, 1)
      if(par == 1) {
        legend("topleft",
               legend = sprintf(tols, fmt = '%#.3f'),
               cex = 1, pt.cex = 1.5,
               col = colorspace::adjust_transparency(cols, alpha = FALSE),
               pch = 19, bty = "n", horiz = TRUE)
      }
    } else {
      plot(rep(cv4abc.out$true[, par], times = numtols),
           unlist(cv4abc.out$estim[, columns]),
           xlab = NA, ylab = NA,
           col = rep(cols, each = nval), pch = pch)
      mtext("True value", side = 1, line = 1.75, cex = 0.75)
      mtext("Est. value", side = 2, line = 1.75, cex = 0.75)
      mtext(parnames[par], side = 3, line = 0, cex = 0.75, font = 2)
      abline(0, 1)
      if(par == 1) {
        legend("topleft",
               legend = sprintf(tols, fmt = '%#.3f'),
               cex = 1, pt.cex = 1.5,
               col = colorspace::adjust_transparency(cols, alpha = FALSE),
               pch = 19, bty = "o",
               box.lwd = 0, bg = "white", horiz = TRUE)
      }
    }
  }
  invisible(par(opar))
}

# plot accepted sims against targets
plot.targets <- function(x, model_targets, obs_targets, ...) {
  require(ggplot2); require(Hmisc)
  #if (!inherits(x, "abc")){
    #stop("Use only with objects of class \"abc\".", call. = FALSE)}
  if (!is.matrix(obs_targets) && !is.data.frame(obs_targets) && !is.vector(obs_targets)){
    stop("'param' has to be a matrix, data.frame or vector.",
         call. = F)}
  if (is.null(dim(obs_targets))){
    obs_targets <- matrix(obs_targets, ncol = 1, dimnames = list(colnames(x$ss), "Target"))}
  if (is.data.frame(obs_targets)){
    obs_targets<- as.matrix(obs_targets)}
  nt <- x$numstat
  numsim <- nrow(model_targets)
  tarnames <- x$names$statistics.names
  cond <- isTRUE(c(match(colnames(obs_targets), tarnames),
                   match(tarnames, colnames(model_targets))))
  if (cond) {
    stop("'targets' and 'model_targets' are not compatible; paramater names are different.",
         call. = FALSE)}
  sel_sims <- x$ss
  df <- rbindlist(list(
    as.data.table(sel_sims)[, Scen := "Selected"],
    as.data.table(model_targets)[, Scen := "AllSims"]))
  suppressWarnings({mnDF <- melt(
    df[Scen == "Selected",
       lapply(.SD, weighted.mean, w = 1/x$dist[x$region]),
       .SDcols = tarnames])})
  mnDF[, Scen := "Selected"]
  suppressWarnings({allDF <- melt(
    df[Scen == "AllSims",
       lapply(.SD, weighted.mean, w = 1/x$dist[!na_idx]),
       .SDcols = tarnames])})
  allDF[, Scen := "AllSims"]
  mnDF <- rbindlist(
    list(mnDF,allDF))
  df <- melt(df, id.vars = "Scen")
  p1 <- ggplot(data = df, aes(x = value, fill = Scen, group = Scen)) +
    facet_wrap(~variable, scales = "free") +
    geom_histogram(col = "grey70", bins = 50) +
    geom_vline(data = mnDF,
               aes(xintercept = value, colour = Scen)) +
    scale_y_continuous(expand = c(0.01, 0.01), ...) +
    geom_vline(data = as.data.table(obs_targets)[,variable := unique(unique(df$variable))],
               aes(xintercept = Target)) +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          legend.title = element_blank()) +
    labs(x = NULL, y = NULL)
  return(print(p1))
}
