graceT = function (argvals, y, outcome, id, group = NULL, tolerance = 0.001, 
          maxiter = 20, norder = 6, ngrid = 5, Lfdobj = 3, lambda = 10^(4), 
          verbose = FALSE, plots = FALSE, figdir = file.path(getwd(), 
                                                             paste("Test", Sys.time(), sep = "_")), shift = TRUE) 
{
  require(lme4)
  require(fda)
  if (plots) {
    safe.dir.create <- function(path) {
      dirTest <- function(x) !is.na(isdir <- file.info(x)$isdir) & 
        isdir
      if (!dirTest(path) && !dir.create(path)) 
        stop(gettextf("cannot create directory '%s'", 
                      path), domain = NA)
    }
    if (!file.exists(figdir)) 
      safe.dir.create(figdir)
  }
  if (is.null(group)) 
    group <- rep(1, length(y))
  if (any(is.na(y))) 
    stop("NA y values not allowed.")
  if (any(is.na(argvals))) 
    stop("NA argvals not allowd.")
  if (any(is.na(id))) 
    stop("NA ids not allowed.")
  if (!file.exists(figdir)) 
    dir.create(figdir)
  dd <- data.frame(id = id, group = group, outcome = outcome, 
                   argvals = as.numeric(argvals), y = as.numeric(y))
  dd$alpha0 <- 0
  dd$alpha1 <- 0
  dd$gamma0 <- 0
  sigma <- c()
  rss.0 <- with(dd, sum((y - 0)^2))
  rss.1 <- with(dd, sum((y - mean(y))^2))
  i <- 1
  while (abs(rss.0 - rss.1)/rss.0 > tolerance & i <= maxiter) {
    dd <- dd[with(dd, order(argvals + gamma0)), ]
    fits <- lapply(unique(outcome), function(oc) {
      subs <- subset(dd, outcome == oc)
      rng <- with(subs, range(argvals + gamma0))
      grid <- seq(rng[1], rng[2], length = ngrid)
      nbasis <- ngrid + norder - 2
      wbasis <- create.bspline.basis(rng, nbasis, norder, 
                                     grid)
      cvec0 <- matrix(0, nbasis, 1)
      Wfd0 <- fd(cvec0, wbasis)
      WfdParobj <- fdPar(Wfd0, Lfdobj, lambda)
      subs$alpha0 <- 0
      subs$alpha1 <- 0
      rss.subs.0 <- with(subs, sum((y - 0)^2))
      rss.subs.1 <- with(subs, sum((y - mean(y))^2))
      j <- 1
      while (abs(rss.subs.0 - rss.subs.1)/rss.subs.0 > 
             tolerance & j <= maxiter) {
        trash <- try(capture.output(mfit <- with(subs, 
                                                 smooth.monotone(argvals + gamma0, y - alpha0 - 
                                                                   alpha1 * argvals, WfdParobj))), silent = TRUE)
        if (class(trash) == "try-error") 
          trash <- try(capture.output(mfit <- with(subs, 
                                                   smooth.monotone(argvals + gamma0, y - alpha0 - 
                                                                     alpha1 * argvals, WfdParobj, iterlim = 1))), 
                       silent = TRUE)
        if (class(trash) == "try-error") {
          mfit <- lm(y - alpha0 - alpha1 * argvals ~ 
                       argvals + gamma0, subs)
          subs$ghat <- predict(mfit)
        }
        else {
          Wfd <- mfit$Wfdobj
          beta <- mfit$beta
          subs$ghat <- beta[1] + beta[2] * eval.monfd(with(subs, 
                                                           argvals + gamma0), Wfd)
        }
        subs$smooth.resids <- with(subs, y - ghat)
        lmefit <- try(lmer(smooth.resids ~ 1 + (argvals | 
                                                  id), data = subs), silent = TRUE)
        if (class(lmefit) == "try-error") {
          lmefit <- try(lmer(smooth.resids ~ 1 + (1 | 
                                                    id), data = subs), silent = TRUE)
          if (class(lmefit) == "try-error") {
            ranefs <- subs[, c("id", "smooth.resids")]
            ranefs$alpha0 <- 0
            ranefs$alpha1 <- 0
            ranefs <- subset(ranefs, !duplicated(id))
          }
          else {
            ranefs <- ranef(lmefit)$id
            ranefs$id <- rownames(ranefs)
            ranefs$alpha0 <- ranefs[, "(Intercept)"] - 
              mean(ranefs[, "(Intercept)"])
            ranefs$alpha1 <- 0
          }
        }
        else {
          ranefs <- ranef(lmefit)$id
          ranefs$id <- rownames(ranefs)
          ranefs$alpha0 <- ranefs[, "(Intercept)"] - 
            mean(ranefs[, "(Intercept)"])
          ranefs$alpha1 <- ranefs[, "argvals"] - mean(ranefs[, 
                                                             "argvals"])
        }
        subs <- merge(subs[, !colnames(subs) %in% c("alpha0", 
                                                    "alpha1")], ranefs[, c("id", "alpha0", "alpha1")], 
                      by = "id", all.x = TRUE)
        rss.subs.0 <- rss.subs.1
        rss.subs.1 <- with(subs, sum((y - ghat)^2))
        j <- j + 1
      }
      subs$yihat <- with(subs, ghat + alpha0 + alpha1 * 
                           argvals)
      subs$sigma <- sd(as.numeric(with(subs, yihat - y)))
      gkinv <- with(subs, approxfun(ghat, argvals + gamma0, 
                                    rule = 2))
      subs$gamma0new <- with(subs, gkinv(y) - argvals)
      list(monotone = mfit, lme = lmefit, subset = subs)
    })
    dd <- do.call(rbind, lapply(fits, function(x) x$subs))
    if (plots) {
      ggplot(dd, aes(argvals + gamma0, y, group = id)) + 
        geom_line(aes(colour = group), alpha = 1/4) + 
        geom_point(aes(colour = group), size = 1, alpha = 1) + 
        ylim(min(dd$y), max(dd$y)) + geom_line(aes(argvals + 
                                                     gamma0, ghat, group = NULL), size = 1.5) + facet_wrap(group~outcome) + 
        xlab(expression(t + hat(gamma))) + labs(title = paste("Iteration", 
                                                              i, sep = " "))
      ggsave(file = file.path(figdir, paste("iteration_", 
                                            i, ".pdf", sep = "")),width = 10,height = 10)
      ggplot(dd, aes(argvals, smooth.resids, group = id)) + 
        geom_line(aes(colour = group), alpha = 1/4) + 
        geom_point(aes(colour = group), size = 1, alpha = 1) + 
        ylim(min(dd$smooth.resids), max(dd$smooth.resids)) + 
        facet_wrap(~outcome) + xlab("t") + labs(title = paste("Iteration", 
                                                              i, sep = " "))
      ggsave(file = file.path(figdir, paste("resids_iteration_", 
                                            i, ".pdf", sep = "")))
    }
    sigma <- rbind(sigma, cbind(Iteration = i, dd[!duplicated(dd$outcome), 
                                                  c("outcome", "sigma")]))
    if (shift) {
      shifts <- with(dd, aggregate(gamma0new, by = list(id = id), 
                                   FUN = mean, na.rm = TRUE))
      shifts$gamma0 <- shifts$x - mean(shifts$x)
      dd <- merge(dd[, colnames(dd) != "gamma0"], shifts[, 
                                                         c("id", "gamma0")], by = "id", all.x = TRUE)
    }
    else {
      i <- maxiter + 1
    }
    if (plots) {
      ggplot(subset(dd, !duplicated(id)), aes(y = gamma0, 
                                              x = group)) + geom_boxplot() + labs(title = paste("Iteration", 
                                                                                                i, sep = " "))
      ggsave(file = file.path(figdir, paste("gamma0_iteration_", 
                                            i, ".pdf", sep = "")))
    }
    rss.0 <- rss.1
    rss.1 <- with(dd, sum((y - ghat)^2))
    i <- i + 1
  }
  names(fits) <- unique(outcome)
  fit <- list(fits = fits, sigma = sigma)
  class(fit) <- "grace"
  fit
}
