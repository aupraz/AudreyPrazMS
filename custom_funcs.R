# Function for plotting ABC results with the rejection method posteriors
plot.rej <- function(x, param, dens_kernel = "epanechnikov", ...) {
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
                       kernel = dens_kernel)
    rej.d <- density(rej[, i],
                     from = min(rej[, i], na.rm = T),
                     to = max(rej[, i], na.rm = TRUE),
                     kernel = dens_kernel)
    rej.dl <- density(rej[, i],
                      from = quantile(rej[, i], 0.05, na.rm = TRUE),
                      to = quantile(rej[, i], 0.95, na.rm = TRUE),
                      kernel = dens_kernel)
    myxlim <- range(c(prior.d$x, rej.d$x))
    myylim <- range(c(prior.d$y, rej.d$y))
    myxlab <- parnames[i]
    plot(prior.d, main = myxlab, col = "#1b9e77", lty = 2,
         lwd = 1, xlab = "", ylim = myylim, xlim = myxlim)
    title(sub = paste0(paste("N =", prior.d$n, "  Bandwidth =", formatC(prior.d$bw)), "; ",
                       sprintf("mean = %0.3f; sd = %0.3f",
                               mean(param[, i],na.rm = T), sd(param[, i],na.rm = T)), "\n",
                       sprintf("meanSel = %0.3f; sdSel = %0.3f",
                               mean(rej[, i],na.rm = T), sd(rej[, i],na.rm = T))))
    lines(rej.d, col = "#d95f02", lty = 1, lwd = 1)
    rug(rej[, i], side = 1, ticksize = 0.02, lwd = 1, col = "#d95f02")
    abline(v = range(rej.dl$x), col = "grey90", lty = 1, lwd = 1)
    if (i == np) {
      legend("bottomright",
             c("Prior", "Rej."),
             lty=c(1,1),
             col = c("#1b9e77", "#d95f02"),
             inset = c(0,1), xpd = TRUE, horiz = TRUE,
             bty="n"
      )
    }
  }
  par(mfcol = c(1, 1))
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

# Default ABC function with some extra warnings/context for errors
abc.stu <- function (target, param, sumstat, tol, method, hcorr = TRUE, scaling = "sd",
                     transf = "none", logit.bounds = c(0, 0), subset = NULL, kernel = "epanechnikov",
                     numnet = 10, sizenet = 5, lambda = c(1e-04, 0.001, 0.01),
                     trace = FALSE, maxit = 500, ...) {
  call <- match.call()
  if (missing(target))
    stop("'target' is missing")
  if (missing(param))
    stop("'param' is missing")
  if (missing(sumstat))
    stop("'sumstat' is missing")
  if (!is.matrix(param) && !is.data.frame(param) && !is.vector(param))
    stop("'param' has to be a matrix, data.frame or vector.")
  if (!is.matrix(sumstat) && !is.data.frame(sumstat) && !is.vector(sumstat))
    stop("'sumstat' has to be a matrix, data.frame or vector.")
  if (missing(tol))
    stop("'tol' is missing")
  if (missing(method))
    stop("'method' is missing with no default")
  if (!any(method == c("rejection", "loclinear", "neuralnet",
                       "ridge"))) {
    stop("Method must be rejection, loclinear, or neuralnet or ridge")
  }
  if (!any(scaling == c("sd", "mad"))) {
    stop("Scaling must be sd or mad")
  }
  if (method == "rejection"){
    rejmethod <- TRUE
  }else {rejmethod <- FALSE}
  if (!any(kernel == c("gaussian", "epanechnikov", "rectangular",
                       "triangular", "biweight", "cosine"))) {
    kernel <- "epanechnikov"
    warning("Kernel is incorrectly defined. Setting to default (Epanechnikov)")
  }
  if (is.data.frame(param))
    param <- as.matrix(param)
  if (is.data.frame(sumstat))
    sumstat <- as.matrix(sumstat)
  if (is.list(target))
    target <- unlist(target)
  if (is.vector(sumstat))
    sumstat <- matrix(sumstat, ncol = 1)
  if (length(target) != dim(sumstat)[2])
    stop("Number of summary statistics in 'target' has to be the same as in 'sumstat'.")
  nss <- length(sumstat[1, ])
  cond1 <- !any(as.logical(apply(sumstat, 2, function(x) length(unique(x)) - 1)))
  if (cond1)
    stop("Zero variance in the summary statistics.")
  ltransf <- length(transf)
  if (is.vector(param)) {
    numparam <- 1
    param <- matrix(param, ncol = 1)
  } else {
    numparam <- dim(param)[2]}
  for (i in 1:ltransf) {
    if (sum(transf[i] == c("none", "log", "logit")) == 0) {
      stop("Transformations must be none, log, or logit.")
    }
    if (transf[i] == "logit") {
      if (logit.bounds[i, 1] >= logit.bounds[i, 2]) {
        stop("Logit bounds are incorrect.")
      }
    }
  }
  if (rejmethod) {
    if (!all(transf == "none")) {
      warning("No transformation is applied when the simple rejection is used.",
              call. = F)
    }
    transf[1:numparam] <- "none"
  } else {
    if (numparam != ltransf) {
      if (length(transf) == 1) {
        transf <- rep(transf[1], numparam)
        warning("All parameters are \"", transf[1], "\" transformed.",
                sep = "", call. = F)
      }
      else stop("Number of parameters is not the same as number of transformations.",
                sep = "", call. = F)
    }
  }
  gwt <- rep(TRUE, length(sumstat[, 1]))
  gwt[attributes(na.omit(sumstat))$na.action] <- FALSE
  if (missing(subset)){
    subset <- rep(TRUE, length(sumstat[, 1]))}
  gwt <- as.logical(gwt * subset)
  if (!length(colnames(param))) {
    warning("No parameter names are given, using P1, P2, ...")
    paramnames <- paste("P", 1:numparam, sep = "")
  } else {
    paramnames <- colnames(param)}
  if (!length(colnames(sumstat))) {
    warning("No summary statistics names are given, using S1, S2, ...")
    statnames <- paste("S", 1:nss, sep = "")
  } else {
    statnames <- colnames(sumstat)
  }
  scaled.sumstat <- sumstat
  for (j in 1:nss) {
    if (scaling == "mad") {
      scaled.sumstat[, j] <- abc:::normalise(x = sumstat[, j],
                                             y = sumstat[,j][gwt])
    } else {
      sdScale <- function (x, y) {
        if (sd(y) == 0) {
          return(x)
        } else {
          return(x/sd(y))
        }}
      scaled.sumstat[, j] <- sdScale(x = sumstat[, j], y = sumstat[,j][gwt])
    }
  }
  for (j in 1:nss) {
    if (scaling == "mad") {
      target[j] <- abc:::normalise(target[j], sumstat[, j][gwt])
    } else {
      target[j] <- sdScale(target[j], sumstat[, j][gwt])
    }
  }
  sum1 <- 0
  for (j in 1:nss) {
    sum1 <- sum1 + (scaled.sumstat[, j] - target[j])^2
  }
  dist <- sqrt(sum1)
  dist[!gwt] <- floor(max(dist[gwt]) + 10)
  nacc <- ceiling(length(dist) * tol)
  ds <- sort(dist)[nacc]
  wt1 <- (dist <= ds)
  aux <- cumsum(wt1)
  wt1 <- wt1 & (aux <= nacc)
  if (kernel == "gaussian") {
    wt1 <- rep(TRUE, length(dist))
  }
  for (i in 1:numparam) {
    if (transf[i] == "log") {
      if (min(param[, i]) <= 0) {
        cat("log transform: values out of bounds - correcting...")
        x.tmp <- ifelse(param[, i] <= 0, max(param[, i]), param[, i])
        x.tmp.min <- min(x.tmp)
        param[, i] <- ifelse(param[, i] <= 0, x.tmp.min,
                             param[, i])
      }
      param[, i] <- log(param[, i])
    } else if (transf[i] == "logit") {
      if (min(param[, i]) <= logit.bounds[i, 1]) {
        x.tmp <- ifelse(param[, i] <= logit.bounds[i, 1], max(param[, i]), param[, i])
        x.tmp.min <- min(x.tmp)
        param[, i] <- ifelse(param[, i] <= logit.bounds[i,1], x.tmp.min, param[, i])
      }
      if (max(param[, i]) >= logit.bounds[i, 2]) {
        x.tmp <- ifelse(param[, i] >= logit.bounds[i,
                                                   2], min(param[, i]), param[, i])
        x.tmp.max <- max(x.tmp)
        param[, i] <- ifelse(param[, i] >= logit.bounds[i,
                                                        2], x.tmp.max, param[, i])
      }
      param[, i] <- (param[, i] - logit.bounds[i, 1])/(logit.bounds[i,
                                                                    2] - logit.bounds[i, 1])
      param[, i] <- log(param[, i]/(1 - param[, i]))
    }
  }
  ss <- sumstat[wt1, ]
  unadj.values <- param[wt1, ]
  statvar <- as.logical(apply(cbind(sumstat[wt1, ]), 2, function(x) length(unique(x)) - 1))
  cond2 <- any(!statvar)
  if (cond2 && !rejmethod)
    stop(paste0("Zero variance in the summary statistics in the selected region. Try: checking summary statistics, choosing larger tolerance, or rejection method.\n",
                paste(colnames(ss)[!statvar], collapse = ", ")))
  if (rejmethod) {
    if (cond2)
      warning(paste0("Zero variance in the summary statistics in the selected region. Try: checking summary statistics, choosing larger tolerance.\n",
                     paste(colnames(ss)[!statvar], collapse = ", ")))
    weights <- rep(1, length = sum(wt1))
    adj.values <- NULL
    residuals <- NULL
    lambda <- NULL
  } else {
    if (cond2)
      cat("Warning messages:\nStatistic(s)", statnames[!statvar],
          "has/have zero variance in the selected region.\nConsider using larger tolerance or the rejection method or discard this/these statistics, which might solve the collinearity problem in 'lsfit'.\n",
          sep = ", ")
    if (kernel == "epanechnikov")
      weights <- 1 - (dist[wt1]/ds)^2
    if (kernel == "rectangular")
      weights <- dist[wt1]/ds
    if (kernel == "gaussian")
      weights <- 1/sqrt(2 * pi) * exp(-0.5 * (dist/(ds/2))^2)
    if (kernel == "triangular")
      weights <- 1 - abs(dist[wt1]/ds)
    if (kernel == "biweight")
      weights <- (1 - (dist[wt1]/ds)^2)^2
    if (kernel == "cosine")
      weights <- cos(pi/2 * dist[wt1]/ds)
    if (method == "loclinear") {
      fit1 <- lsfit(scaled.sumstat[wt1, ], param[wt1, ],
                    wt = weights)
      pred <- t(structure(cbind(fit1$coefficients)[fit1$qr$pivot,
      ], names = names(fit1$coefficients))) %*% c(1,
                                                  target)
      pred <- matrix(pred, ncol = numparam, nrow = sum(wt1),
                     byrow = TRUE)
      residuals <- param[wt1, ] - t(t(structure(cbind(fit1$coefficients)[fit1$qr$pivot,
      ], names = names(fit1$coefficients))) %*% t(cbind(1,
                                                        scaled.sumstat[wt1, ])))
      residuals <- cbind(residuals)
      the_m <- apply(residuals, FUN = mean, 2)
      residuals <- sapply(1:numparam, FUN = function(x) {
        residuals[, x] - the_m[x]
      })
      pred <- sapply(1:numparam, FUN = function(x) {
        pred[, x] + the_m[x]
      })
      sigma2 <- apply(as.matrix(residuals), FUN = function(x) {
        sum((x)^2 * weights)/sum(weights)
      }, MARGIN = 2)
      aic <- sum(wt1) * sum(log(sigma2)) + 2 * (nss + 1) *
        numparam
      bic <- sum(wt1) * sum(log(sigma2)) + log(sum(wt1)) *
        (nss + 1) * numparam
      if (hcorr == TRUE) {
        fit2 <- lsfit(scaled.sumstat[wt1, ], log(residuals^2),
                      wt = weights)
        auxaux <- t(structure(cbind(fit2$coefficients)[fit2$qr$pivot,
        ], names = names(fit2$coefficients))) %*% c(1,
                                                    target)
        pred.sd <- sqrt(exp(auxaux))
        pred.sd <- matrix(pred.sd, nrow = sum(wt1), ncol = numparam,
                          byrow = T)
        pred.si <- t(t(structure(cbind(fit2$coefficients)[fit2$qr$pivot,
        ], names = names(fit2$coefficients))) %*% t(cbind(1,
                                                          scaled.sumstat[wt1, ])))
        pred.si <- sqrt(exp(pred.si))
        adj.values <- pred + (pred.sd * residuals)/pred.si
        residuals <- (pred.sd * residuals)/pred.si
      }
      else {
        adj.values <- pred + residuals
      }
      colnames(adj.values) <- colnames(unadj.values)
      lambda <- NULL
    }
    if (method == "neuralnet") {
      linout <- TRUE
      param.mad <- c()
      for (i in 1:numparam) {
        param.mad[i] <- stats::mad(param[, i][gwt])
        param[, i] <- abc:::normalise(param[, i], param[, i][gwt])
      }
      lambda <- sample(lambda, numnet, replace = T)
      fv <- array(dim = c(sum(wt1), numparam, numnet))
      pred <- matrix(nrow = numparam, ncol = numnet)
      for (i in 1:numnet) {
        fit1 <- nnet(scaled.sumstat[wt1, ], param[wt1,
        ], weights = weights, decay = lambda[i], size = sizenet,
        trace = trace, linout = linout, maxit = maxit,
        ...)
        cat(i)
        fv[, , i] <- fit1$fitted.values
        pred[, i] <- predict(fit1, data.frame(rbind(target)))
      }
      cat("\n")
      pred.med <- apply(pred, 1, median)
      pred.med <- matrix(pred.med, nrow = sum(wt1), ncol = numparam,
                         byrow = T)
      fitted.values <- apply(fv, c(1, 2), median)
      residuals <- param[wt1, ] - fitted.values
      if (hcorr == TRUE) {
        pred2 <- matrix(nrow = numparam, ncol = numnet)
        fv2 <- array(dim = c(sum(wt1), numparam, numnet))
        for (i in 1:numnet) {
          fit2 <- nnet(scaled.sumstat[wt1, ], log(residuals^2),
                       weights = weights, decay = lambda[i], size = sizenet,
                       trace = trace, linout = linout, ...)
          cat(i)
          fv2[, , i] <- fit2$fitted.values
          pred2[, i] <- predict(fit2, data.frame(rbind(target)))
        }
        cat("\n")
        pred.sd <- sqrt(exp(apply(pred2, 1, median)))
        pred.sd <- matrix(pred.sd, nrow = sum(wt1), ncol = numparam,
                          byrow = T)
        fv.sd <- sqrt(exp(apply(fv2, c(1, 2), median)))
        adj.values <- pred.med + (pred.sd * residuals)/fv.sd
        residuals <- (pred.sd * residuals)/fv.sd
      }
      else {
        adj.values <- pred.med + residuals
      }
      colnames(adj.values) <- colnames(unadj.values)
      for (i in 1:numparam) {
        adj.values[, i] <- adj.values[, i] * param.mad[i]
      }
    }
    if (method == "ridge") {
      param.mad <- c()
      for (i in 1:numparam) {
        param.mad[i] <- stats::mad(param[, i][gwt])
        param[, i] <- abc:::normalise(param[, i], param[, i][gwt])
      }
      numnet <- length(lambda)
      fv <- array(dim = c(sum(wt1), numparam, numnet))
      pred <- matrix(nrow = numparam, ncol = numnet)
      mataux <- sqrt(diag(weights))
      paramaux <- as.matrix(mataux %*% param[wt1, ])
      scaledaux <- mataux %*% scaled.sumstat[wt1, ]
      for (parcur in (1:numparam)) {
        #message(parcur)
        fit1 <- MASS::lm.ridge(paramaux[, parcur] ~ scaledaux,
                               lambda = lambda)
        for (i in 1:numnet) {
          fv[, parcur, i] <- drop(cbind(1, scaled.sumstat[wt1,
          ]) %*% (rbind(coef(fit1))[i, ]))
          pred[parcur, i] <- drop(c(1, target) %*% (rbind(coef(fit1))[i,
          ]))
        }
      }
      pred.med <- apply(pred, 1, median)
      pred.med <- matrix(pred.med, nrow = sum(wt1), ncol = numparam,
                         byrow = T)
      fitted.values <- apply(fv, c(1, 2), median)
      residuals <- param[wt1, ] - fitted.values
      if (hcorr == TRUE) {
        pred2 <- matrix(nrow = numparam, ncol = numnet)
        fv2 <- array(dim = c(sum(wt1), numparam, numnet))
        for (parcur in (1:numparam)) {
          lresidaux <- (mataux %*% (log(residuals[, parcur]^2)))
          fit2 <- lm.ridge(lresidaux ~ scaledaux, lambda = lambda)
          for (i in 1:numnet) {
            fv2[, parcur, i] <- drop(cbind(1, scaled.sumstat[wt1,
            ]) %*% (rbind(coef(fit2))[i, ]))
            pred2[parcur, i] <- drop(c(1, target) %*%
                                       (rbind(coef(fit2))[i, ]))
          }
        }
        cat("\n")
        pred.sd <- sqrt(exp(apply(pred2, 1, median)))
        pred.sd <- matrix(pred.sd, nrow = sum(wt1), ncol = numparam,
                          byrow = T)
        fv.sd <- sqrt(exp(apply(fv2, c(1, 2), median)))
        adj.values <- pred.med + (pred.sd * residuals)/fv.sd
        residuals <- (pred.sd * residuals)/fv.sd
      }
      else {
        adj.values <- pred.med + residuals
      }
      colnames(adj.values) <- colnames(unadj.values)
      for (i in 1:numparam) {
        adj.values[, i] <- adj.values[, i] * param.mad[i]
      }
    }
  }
  if (numparam == 1) {
    unadj.values <- matrix(unadj.values, ncol = 1)
    if (method != "rejection") {
      adj.values <- matrix(adj.values, ncol = 1)
      residuals <- matrix(residuals, ncol = 1)
    }
  }
  for (i in 1:numparam) {
    if (transf[i] == "log") {
      unadj.values[, i] <- exp(unadj.values[, i])
      adj.values[, i] <- exp(adj.values[, i])
    }
    else if (transf[i] == "logit") {
      unadj.values[, i] <- exp(unadj.values[, i])/(1 +
                                                     exp(unadj.values[, i]))
      unadj.values[, i] <- unadj.values[, i] * (logit.bounds[i,
                                                             2] - logit.bounds[i, 1]) + logit.bounds[i, 1]
      adj.values[, i] <- exp(adj.values[, i])/(1 + exp(adj.values[,
                                                                  i]))
      adj.values[, i] <- adj.values[, i] * (logit.bounds[i,
                                                         2] - logit.bounds[i, 1]) + logit.bounds[i, 1]
    }
  }
  abc:::abc.return(transf, logit.bounds, method, call, numparam,
                   nss, paramnames, statnames, unadj.values, adj.values,
                   ss, weights, residuals, dist, wt1, gwt, lambda, hcorr,
                   aic, bic)
}

# plot accepted sims against targets
plot.targets <- function(x, model_targets, obs_targets, ...) {
  require(ggplot2); require(Hmisc); require(data.table)
  if (!inherits(x, "abc")){
    stop("Use only with objects of class \"abc\".", call. = FALSE)}
  if (!is.matrix(obs_targets) && !is.data.frame(obs_targets) && !is.vector(obs_targets)){
    stop("'obs_targets' has to be a matrix, data.frame or vector.",
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
    stop("'obs_targets' and 'model_targets' are not compatible; paramater names are different.",
         call. = FALSE)}
  sel_sims <- x$ss
  df <- rbindlist(list(
    as.data.table(sel_sims)[, Scen := "Selected"],
    as.data.table(model_targets)[, Scen := "AllSims"][!na_idx, ]))
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
               aes(xintercept = V1)) +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          legend.title = element_blank()) +
    labs(x = NULL, y = NULL)
  return(print(p1))
}

# Truncated distribution functions

## normal
rnormt <- function(n, range, mean, sd, ...) {
  F.a <- pnorm(min(range), mean = mean, sd = sd, ...)
  F.b <- pnorm(max(range), mean = mean, sd = sd, ...)
  u <- runif(n, min = F.a, max = F.b)
  qnorm(p = u, mean = mean, sd = sd, ...)
}

## log-normal
rlnormt <- function(n, range, meanlog, sdlog, ...) {
  F.a <- plnorm(min(range), meanlog = meanlog, sdlog = sdlog, ...)
  F.b <- plnorm(max(range), meanlog = meanlog, sdlog = sdlog, ...)
  u <- runif(n, min = F.a, max = F.b)
  qlnorm(u, meanlog = meanlog, sdlog = sdlog, ...)
}

## gamma
rgammat <- function(n, range, shape, rate, ...) {
  F.a <- pgamma(min(range), shape = shape, rate = rate,...)
  F.b <- pgamma(max(range), shape = shape, rate = rate, ...)
  u <- runif(n, min = F.a, max = F.b)
  qgamma(u, shape = shape,rate = rate, ...)
}

## beta
rbetat <- function(n, range, shape1, shape2, ...) {
  F.a <- pbeta(min(range), shape1, shape2,...)
  F.b <- pbeta(max(range), shape1, shape2,...)
  u <- runif(n, min = F.a, max = F.b)
  qbeta(u, shape1, shape2,...)
}

## triangular
rtrianglet <- function(n, range, min, max, ...) {
  F.a <- ptriang(min(range), min = min, max = max, ...)
  F.b <- ptriang(max(range), min = min, max = max, ...)
  u <- runif(n, min = F.a, max = F.b)
  qtriang(u, min = min, max = max, ...)
}

## poisson
rpoist <- function(n, range, lambda, ...) {
  F.a <- ppois(min(range), lambda, ...)
  F.b <- ppois(max(range), lambda, ...)
  u <- runif(n, min = F.a, max = F.b)
  qpois(p = u, lambda, ...)
}

## negbin
rnbt <- function(n, range, size, mu, ...) {
  F.a <- pnbinom(min(range), size = size, mu = mu, ...)
  F.b <- pnbinom(max(range), size = size, mu = mu, ...)
  u <- runif(n, min = F.a, max = F.b)
  qnbinom(u, size = size, mu = mu, ...)
}

