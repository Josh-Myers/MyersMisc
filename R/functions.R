# package functions

#' MyCalPlot function
#' my hack of validate.plot.default in rms for my 2 panel plots to match with val.prob
MyCalPlot <-
  function (p, y, logit, group, weights = rep(1, length(y)), normwt = FALSE,
            pl = TRUE, smooth = TRUE, logistic.cal = TRUE, xlab = "Predicted Probability",
            ylab = "Actual Probability", lim = c(0, 1), m, g, cuts, emax.lim = c(0,
                                                                                 1), legendloc = lim[1] + c(0.55 * diff(lim), 0.27 * diff(lim)),
            statloc = c(0, 0.99), riskdist = "calibrated", cex = 0.7,
            mkh = 0.02, connect.group = FALSE, connect.smooth = TRUE,
            g.group = 4, evaluate = 100, nmin = 0)
  {
    if (missing(p))
      p <- plogis(logit)
    else logit <- qlogis(p)
    if (length(p) != length(y))
      stop("lengths of p or logit and y do not agree")
    names(p) <- names(y) <- names(logit) <- NULL
    Spi <- function(p, y) {
      z <- sum((y - p) * (1 - 2 * p))/sqrt(sum((1 - 2 * p) *
                                                 (1 - 2 * p) * p * (1 - p)))
      P <- 2 * pnorm(-abs(z))
      c(Z = z, P = P)
    }
    if (!missing(group)) {
      if (length(group) == 1 && is.logical(group) && group)
        group <- rep("", length(y))
      if (!is.factor(group))
        group <- if (is.logical(group) || is.character(group))
          as.factor(group)
      else cut2(group, g = g.group)
      names(group) <- NULL
      nma <- !(is.na(p + y + weights) | is.na(group))
      ng <- length(levels(group))
    }
    else {
      nma <- !is.na(p + y + weights)
      ng <- 0
    }
    logit <- logit[nma]
    y <- y[nma]
    p <- p[nma]
    if (ng > 0) {
      group <- group[nma]
      weights <- weights[nma]
      return(val.probg(p, y, group, evaluate, weights, normwt,
                       nmin))
    }
    if (length(unique(p)) == 1) {
      P <- mean(y)
      Intc <- qlogis(P)
      n <- length(y)
      D <- -1/n
      L01 <- -2 * sum(y * logit - logb(1 + exp(logit)), na.rm = TRUE)
      L.cal <- -2 * sum(y * Intc - logb(1 + exp(Intc)), na.rm = TRUE)
      U.chisq <- L01 - L.cal
      U.p <- 1 - pchisq(U.chisq, 1)
      U <- (U.chisq - 1)/n
      Q <- D - U
      spi <- unname(Spi(p, y))
      stats <- c(0, 0.5, 0, D, 0, 1, U, U.chisq, U.p, Q, mean((y -
                                                                 p[1])^2), Intc, 0, rep(abs(p[1] - P), 2), spi)
      names(stats) <- c("Dxy", "C (ROC)", "R2", "D", "D:Chi-sq",
                        "D:p", "U", "U:Chi-sq", "U:p", "Q", "Brier", "Intercept",
                        "Slope", "Emax", "Eavg", "S:z", "S:p")
      return(stats)
    }
    i <- !is.infinite(logit)
    nm <- sum(!i)
    if (nm > 0)
      warning(paste(nm, "observations deleted from logistic calibration due to probs. of 0 or 1"))
    f.fixed <- lrm.fit(logit[i], y[i], initial = c(0, 1), maxit = 1L)
    f.recal <- lrm.fit(logit[i], y[i])
    stats <- f.fixed$stats
    n <- stats["Obs"]
    predprob <- seq(emax.lim[1], emax.lim[2], by = 5e-04)
    lt <- f.recal$coef[1] + f.recal$coef[2] * qlogis(predprob)
    calp <- plogis(lt)
    emax <- max(abs(predprob - calp))
    Sm <- lowess(p, y, iter = 0)
    cal.smooth <- approx(Sm, xout = p, ties = mean)$y
    eavg <- mean(abs(p - cal.smooth))
    if (pl) {
      plot(0.5, 0.5, xlim = lim, ylim = lim, type = "n", xlab = xlab,
           ylab = ylab)
      abline(0, 1, lty = 2)
      lt <- 2
      leg <- "Ideal"
      marks <- -1
      if (logistic.cal) {
        lt <- c(lt, 1)
        leg <- c(leg, "Logistic calibration")
        marks <- c(marks, -1)
      }
      if (smooth) {
        if (connect.smooth) {
          lines(Sm, lty = 1)
          lt <- c(lt, 3)
          marks <- c(marks, -1)
        }
        else {
          points(Sm)
          lt <- c(lt, 0)
          marks <- c(marks, 1)
        }
        leg <- c(leg, "Actual")
      }
      if (!missing(m) | !missing(g) | !missing(cuts)) {
        if (!missing(m))
          q <- cut2(p, m = m, levels.mean = TRUE, digits = 7)
        else if (!missing(g))
          q <- cut2(p, g = g, levels.mean = TRUE, digits = 7)
        else if (!missing(cuts))
          q <- cut2(p, cuts = cuts, levels.mean = TRUE,
                    digits = 7)
        means <- as.numeric(levels(q))
        prop <- tapply(y, q, function(x) mean(x, na.rm = TRUE))
        points(means, prop, pch = 2)
        if (connect.group) {
          lines(means, prop)
          lt <- c(lt, 1)
        }
        else lt <- c(lt, 0)
        leg <- c(leg, "Grouped observations")
        marks <- c(marks, 2)
      }
    }
    lr <- stats["Model L.R."]
    p.lr <- stats["P"]
    D <- (lr - 1)/n
    L01 <- -2 * sum(y * logit - logb(1 + exp(logit)), na.rm = TRUE)
    U.chisq <- L01 - f.recal$deviance[2]
    p.U <- 1 - pchisq(U.chisq, 2)
    U <- (U.chisq - 2)/n
    Q <- D - U
    Dxy <- stats["Dxy"]
    C <- stats["C"]
    R2 <- stats["R2"]
    B <- mean((p - y)^2)
    spi <- unname(Spi(p, y))
    stats <- c(Dxy, C, R2, D, lr, p.lr, U, U.chisq, p.U, Q, B,
               f.recal$coef, emax, spi)
    names(stats) <- c("Dxy", "C (ROC)", "R2", "D", "D:Chi-sq",
                      "D:p", "U", "U:Chi-sq", "U:p", "Q", "Brier", "Intercept",
                      "Slope", "Emax", "S:z", "S:p")
    stats <- c(stats, c(Eavg = eavg))
    if (pl) {
      logit <- seq(-7, 7, length = 200)
      prob <- plogis(logit)
      pred.prob <- f.recal$coef[1] + f.recal$coef[2] * logit
      pred.prob <- plogis(pred.prob)
      if (logistic.cal)
        lines(prob, pred.prob, lty = 1)
      lp <- legendloc
      if (!is.logical(lp)) {
        if (!is.list(lp))
          lp <- list(x = lp[1], y = lp[2])
        legend(lp, leg, lty = c(2,1), pch = marks, cex = cex,
               bty = "n")
      }
      if (!is.logical(statloc)) {
        dostats <- c(1, 2, 3, 4, 7, 10, 11, 12, 13, 14, 15,
                     16)
        leg <- format(names(stats)[dostats])
        leg <- paste(leg, ":", format(stats[dostats]), sep = "")
        if (!is.list(statloc))
          statloc <- list(x = statloc[1], y = statloc[2])
        text(statloc, paste(format(names(stats[dostats])),
                            collapse = "\n"), adj = c(0, 1), cex = cex)
        text(statloc$x + 0.225 * diff(lim), statloc$y, paste(format(round(stats[dostats],
                                                                          3)), collapse = "\n"), adj = c(1, 1), cex = cex)
      }
      if (is.character(riskdist)) {
        if (riskdist == "calibrated") {
          x <- f.recal$coef[1] + f.recal$coef[2] * qlogis(p)
          x <- plogis(x)
          x[p == 0] <- 0
          x[p == 1] <- 1
        }
        else x <- p
        bins <- seq(lim[1], lim[2], length = 101)
        x <- x[x >= lim[1] & x <= lim[2]]
        f <- table(cut(x, bins))
        j <- f > 0
        bins <- (bins[-101])[j]
        f <- f[j]
        f <- lim[1] + 0.15 * diff(lim) * f/max(f)
        segments(bins, 0, bins, f)
      }
    }
    stats
  }
#' MyValPlot function
#' My hack val.prob for my plot - set line type and legend etc
MyValPlot <-
  function (x, xlab, ylab, xlim, ylim, legend = TRUE, subtitles = TRUE,
            scat1d.opts = NULL, ...)
  {
    at <- attributes(x)
    if (missing(ylab))
      ylab <- if (at$model == "lr")
        "Actual Probability"
    else paste("Observed", at$yvar.name)
    if (missing(xlab)) {
      if (at$model == "lr") {
        xlab <- paste("Predicted Pr{", at$yvar.name, sep = "")
        if (at$non.slopes == 1) {
          xlab <- if (at$lev.name == "TRUE")
            paste(xlab, "}", sep = "")
          else paste(xlab, "=", at$lev.name, "}", sep = "")
        }
        else xlab <- paste(xlab, ">=", at$lev.name, "}",
                           sep = "")
      }
      else xlab <- paste("Predicted", at$yvar.name)
    }
    p <- x[, "predy"]
    p.app <- x[, "calibrated.orig"]
    p.cal <- x[, "calibrated.corrected"]
    if (missing(xlim) & missing(ylim))
      xlim <- ylim <- range(c(p, p.app, p.cal), na.rm = TRUE)
    else {
      if (missing(xlim))
        xlim <- range(p)
      if (missing(ylim))
        ylim <- range(c(p.app, p.cal, na.rm = TRUE))
    }
    plot(p, p.app, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
         type = "n", ...)
    predicted <- at$predicted
    err <- NULL
    if (length(predicted)) {
      s <- !is.na(p + p.cal)
      err <- predicted - approx(p[s], p.cal[s], xout = predicted,
                                ties = mean)$y
      cat("\nn=", n <- length(err), "   Mean absolute error=",
          round(mae <- mean(abs(err), na.rm = TRUE), 3), "   Mean squared error=",
          round(mean(err^2, na.rm = TRUE), 5), "\n0.9 Quantile of absolute error=",
          round(quantile(abs(err), 0.9, na.rm = TRUE), 3),
          "\n\n", sep = "")
      if (subtitles)
        title(sub = paste("Mean absolute error=", round(mae,
                                                        3), " n=", n, sep = ""), cex = 0.65, adj = 1)
      do.call("scat1d", c(list(x = predicted), scat1d.opts))
    }
    lines(p, p.app, lty = 3)
    lines(p, p.cal, lty = 1)
    abline(a = 0, b = 1, lty = 2)
    if (subtitles)
      title(sub = paste("B=", at$B, "repetitions,", at$method),
            adj = 0)
    if (!(is.logical(legend) && !legend)) {
      if (is.logical(legend))
        legend <- list(x = xlim[1] + 0.55 * diff(xlim), y = ylim[1] +
                         0.32 * diff(ylim))
      legend(legend, c("Ideal", "Apparent", "Bias-corrected"),
             lty = c(2, 3, 1), bty = "n")
    }
    invisible(err)
  }
#' my val.prb ci.2 function
#' My hack of steyerberg's cal.plot function - which is an adapataion of the val.prob function from rms
#' The main thing was to get all hist bins poining up to match plot(calibrate) - I had 2 panel cal plot - needed to match the bottom histograms
#' couldn't do that with val.prob because no way to control space under hist plot - could do that with this function
#' changed degree to 1 to match val.prob which is what I really wanted to use in the first place
#' how to install the package
#' library(githubinstall)
#' githubinstall("CalibrationCurves")
#' library(CalibrationCurves)

My.val.prob.ci.2 <-
  function (p, y, logit, group, weights = rep(1, length(y)), normwt = F,
            pl = T, smooth = c("loess", "rcs", F), CL.smooth = "fill",
            CL.BT = F, lty.smooth = 1, col.smooth = "black", lwd.smooth = 1,
            nr.knots = 5, logistic.cal = F, lty.log = 1, col.log = "black",
            lwd.log = 1, xlab = "Predicted Probability", ylab = "Observed proportion",
            xlim = c(-0.02, 1), ylim = c(-0.15, 1), m, g, cuts, emax.lim = c(0,
                                                                             1), legendloc = c(0.5, 0.27), statloc = c(0, 0.85), dostats = T,
            cl.level = 0.95, method.ci = "pepe", roundstats = 2, riskdist = "predicted",
            cex = 0.75, cex.leg = 0.75, connect.group = F, connect.smooth = T,
            g.group = 4, evaluate = 100, nmin = 0, d0lab = "0", d1lab = "1",
            cex.d01 = 0.7, dist.label = 0.04, line.bins = -0.05, dist.label2 = 0.03,
            cutoff, las = 1, length.seg = 1, y.intersp = 1, lty.ideal = 1,
            col.ideal = "grey", lwd.ideal = 1, ...)
  {
    if (smooth[1] == F) {
      smooth <- "F"
    }
    smooth <- match.arg(smooth)
    if (!missing(p))
      if (any(!(p >= 0 | p <= 1))) {
        stop("Probabilities can not be > 1 or < 0.")
      }
    if (missing(p))
      p <- 1/(1 + exp(-logit))
    else logit <- log(p/(1 - p))
    if (!all(c(0, 1) %in% y)) {
      stop("The vector with the binary outcome can only contain the values 0 and 1.")
    }
    if (length(p) != length(y))
      stop("lengths of p or logit and y do not agree")
    names(p) <- names(y) <- names(logit) <- NULL
    if (!missing(group)) {
      if (length(group) == 1 && is.logical(group) && group)
        group <- rep("", length(y))
      if (!is.factor(group))
        group <- if (is.logical(group) || is.character(group))
          as.factor(group)
      else cut2(group, g = g.group)
      names(group) <- NULL
      nma <- !(is.na(p + y + weights) | is.na(group))
      ng <- length(levels(group))
    }
    else {
      nma <- !is.na(p + y + weights)
      ng <- 0
    }
    logit <- logit[nma]
    y <- y[nma]
    p <- p[nma]
    if (ng > 0) {
      group <- group[nma]
      weights <- weights[nma]
      return(val.probg(p, y, group, evaluate, weights, normwt,
                       nmin))
    }
    y <- y[order(p)]
    logit <- logit[order(p)]
    p <- p[order(p)]
    if (length(p) > 5000 & smooth == "loess") {
      warning("Number of observations > 5000, RCS is recommended.",
              immediate. = T)
    }
    if (length(p) > 1000 & CL.BT == T) {
      warning("Number of observations is > 1000, this could take a while...",
              immediate. = T)
    }
    if (length(unique(p)) == 1) {
      P <- mean(y)
      Intc <- log(P/(1 - P))
      n <- length(y)
      D <- -1/n
      L01 <- -2 * sum(y * logit - log(1 + exp(logit)), na.rm = T)
      L.cal <- -2 * sum(y * Intc - log(1 + exp(Intc)), na.rm = T)
      U.chisq <- L01 - L.cal
      U.p <- 1 - pchisq(U.chisq, 1)
      U <- (U.chisq - 1)/n
      Q <- D - U
      stats <- c(0, 0.5, 0, D, 0, 1, U, U.chisq, U.p, Q, mean((y -
                                                                 p[1])^2), Intc, 0, rep(abs(p[1] - P), 2))
      names(stats) <- c("Dxy", "C (ROC)", "R2", "D", "D:Chi-sq",
                        "D:p", "U", "U:Chi-sq", "U:p", "Q", "Brier", "Intercept",
                        "Slope", "Emax", "Eavg", "ECI")
      return(stats)
    }
    i <- !is.infinite(logit)
    nm <- sum(!i)
    if (nm > 0)
      warning(paste(nm, "observations deleted from logistic calibration due to probs. of 0 or 1"))
    i.2 <- i
    f.or <- lrm(y[i] ~ logit[i])
    f <- lrm.fit(logit[i], y[i])
    cl.slope <- confint(f, level = cl.level)[2, ]
    f2 <- lrm.fit(offset = logit[i], y = y[i])
    cl.interc <- confint(f2, level = cl.level)
    stats <- f$stats
    cl.auc <- ci.auc(y, p, cl.level, method.ci)
    n <- stats["Obs"]
    predprob <- seq(emax.lim[1], emax.lim[2], by = 5e-04)
    lt <- f$coef[1] + f$coef[2] * log(predprob/(1 - predprob))
    calp <- 1/(1 + exp(-lt))
    emax <- max(abs(predprob - calp))
    if (pl) {
      plot(0.5, 0.5, xlim = xlim, ylim = ylim, type = "n",
           xlab = xlab, ylab = ylab, las = las, ...)
      clip(0, 1, 0, 1)
      abline(0, 1, lty = lty.ideal, col = col.ideal, lwd = lwd.ideal)
      do.call("clip", as.list(par()$usr))
      lt <- lty.ideal
      lw.d <- lwd.ideal
      all.col <- col.ideal
      leg <- "Ideal"
      marks <- -1
      if (logistic.cal) {
        lt <- c(lt, lty.log)
        lw.d <- c(lw.d, lwd.log)
        all.col <- c(all.col, col.log)
        leg <- c(leg, "Logistic calibration")
        marks <- c(marks, -1)
      }
      if (smooth != "F") {
        all.col <- c(all.col, col.smooth)
      }
      if (smooth == "loess") {
        Sm <- loess(y ~ p, degree = 1)
        Sm <- data.frame(Sm$x, Sm$fitted)
        Sm.01 <- Sm
        if (connect.smooth == T & CL.smooth != "fill") {
          clip(0, 1, 0, 1)
          lines(Sm, lty = lty.smooth, lwd = lwd.smooth,
                col = col.smooth)
          do.call("clip", as.list(par()$usr))
          lt <- c(lt, lty.smooth)
          lw.d <- c(lw.d, lwd.smooth)
          marks <- c(marks, -1)
        }
        else if (connect.smooth == F & CL.smooth != "fill") {
          clip(0, 1, 0, 1)
          points(Sm, col = col.smooth)
          do.call("clip", as.list(par()$usr))
          lt <- c(lt, 0)
          lw.d <- c(lw.d, 1)
          marks <- c(marks, 1)
        }
        if (CL.smooth == T | CL.smooth == "fill") {
          to.pred <- seq(min(p), max(p), length = 200)
          if (CL.BT == T) {
            cat("Bootstrap samples are being generated.\n\n\n")
            res.BT <- replicate(2000, BT.samples(y, p,
                                                 to.pred))
            CL.BT <- apply(res.BT, 1, quantile, c(0.025,
                                                  0.975))
            colnames(CL.BT) <- to.pred
            if (CL.smooth == "fill") {
              clip(0, 1, 0, 1)
              polygon(x = c(to.pred, rev(to.pred)), y = c(CL.BT[2,
                                                                ], rev(CL.BT[1, ])), col = rgb(177, 177,
                                                                                               177, 177, maxColorValue = 255), border = NA)
              if (connect.smooth == T) {
                lines(Sm, lty = lty.smooth, lwd = lwd.smooth,
                      col = col.smooth)
                lt <- c(lt, lty.smooth)
                lw.d <- c(lw.d, lwd.smooth)
                marks <- c(marks, -1)
              }
              else if (connect.smooth == F) {
                points(Sm, col = col.smooth)
                lt <- c(lt, 0)
                lw.d <- c(lw.d, 1)
                marks <- c(marks, 1)
              }
              do.call("clip", as.list(par()$usr))
              leg <- c(leg, "Actual")
            }
            else {
              clip(0, 1, 0, 1)
              lines(to.pred, CL.BT[1, ], lty = 2, lwd = 1,
                    col = col.smooth)
              clip(0, 1, 0, 1)
              lines(to.pred, CL.BT[2, ], lty = 2, lwd = 1,
                    col = col.smooth)
              do.call("clip", as.list(par()$usr))
              leg <- c(leg, "Actual",
                       "CL flexible")
              lt <- c(lt, 2)
              lw.d <- c(lw.d, 1)
              all.col <- c(all.col, col.smooth)
              marks <- c(marks, -1)
            }
          }
          else {
            Sm.0 <- loess(y ~ p, degree = 2)
            cl.loess <- predict(Sm.0, type = "fitted",
                                se = T)
            clip(0, 1, 0, 1)
            if (CL.smooth == "fill") {
              polygon(x = c(Sm.0$x, rev(Sm.0$x)), y = c(cl.loess$fit +
                                                          cl.loess$se.fit * 1.96, rev(cl.loess$fit -
                                                                                        cl.loess$se.fit * 1.96)), col = rgb(177,
                                                                                                                            177, 177, 177, maxColorValue = 255), border = NA)
              if (connect.smooth == T) {
                lines(Sm, lty = lty.smooth, lwd = lwd.smooth,
                      col = col.smooth)
                lt <- c(lt, lty.smooth)
                lw.d <- c(lw.d, lwd.smooth)
                marks <- c(marks, -1)
              }
              else if (connect.smooth == F) {
                points(Sm, col = col.smooth)
                lt <- c(lt, 0)
                lw.d <- c(lw.d, 1)
                marks <- c(marks, 1)
              }
              do.call("clip", as.list(par()$usr))
              leg <- c(leg, "Actual")
            }
            else {
              lines(Sm.0$x, cl.loess$fit + cl.loess$se.fit *
                      1.96, lty = 2, lwd = 1, col = col.smooth)
              lines(Sm.0$x, cl.loess$fit - cl.loess$se.fit *
                      1.96, lty = 2, lwd = 1, col = col.smooth)
              do.call("clip", as.list(par()$usr))
              leg <- c(leg, "Actual",
                       "CL flexible")
              lt <- c(lt, 2)
              lw.d <- c(lw.d, 1)
              all.col <- c(all.col, col.smooth)
              marks <- c(marks, -1)
            }
          }
        }
        else {
          leg <- c(leg, "Actual")
        }
        cal.smooth <- approx(Sm.01, xout = p)$y
        eavg <- mean(abs(p - cal.smooth))
        ECI <- mean((p - cal.smooth)^2) * 100
      }
      if (smooth == "rcs") {
        par(lwd = lwd.smooth, bty = "n", col = col.smooth)
        if (!is.numeric(nr.knots)) {
          stop("Nr.knots must be numeric.")
        }
        if (nr.knots == 5) {
          tryCatch(rcspline.plot(p, y, model = "logistic",
                                 nk = 5, show = "prob", statloc = "none", add = T,
                                 showknots = F, xrange = c(min(na.omit(p)),
                                                           max(na.omit(p))), lty = lty.smooth), error = function(e) {
                                                             warning("The number of knots led to estimation problems, nk will be set to 4.",
                                                                     immediate. = T)
                                                             tryCatch(rcspline.plot(p, y, model = "logistic",
                                                                                    nk = 4, show = "prob", statloc = "none",
                                                                                    add = T, showknots = F, xrange = c(min(na.omit(p)),
                                                                                                                       max(na.omit(p))), lty = lty.smooth), error = function(e) {
                                                                                                                         warning("Nk 4 also led to estimation problems, nk will be set to 3.",
                                                                                                                                 immediate. = T)
                                                                                                                         rcspline.plot(p, y, model = "logistic", nk = 3,
                                                                                                                                       show = "prob", statloc = "none", add = T,
                                                                                                                                       showknots = F, xrange = c(min(na.omit(p)),
                                                                                                                                                                 max(na.omit(p))), lty = lty.smooth)
                                                                                                                       })
                                                           })
        }
        else if (nr.knots == 4) {
          tryCatch(rcspline.plot(p, y, model = "logistic",
                                 nk = 4, show = "prob", statloc = "none", add = T,
                                 showknots = F, xrange = c(min(na.omit(p)),
                                                           max(na.omit(p))), lty = lty.smooth), error = function(e) {
                                                             warning("The number of knots led to estimation problems, nk will be set to 3.",
                                                                     immediate. = T)
                                                             rcspline.plot(p, y, model = "logistic", nk = 3,
                                                                           show = "prob", statloc = "none", add = T,
                                                                           showknots = F, xrange = c(min(na.omit(p)),
                                                                                                     max(na.omit(p))), lty = lty.smooth)
                                                           })
        }
        else if (nr.knots == 3) {
          tryCatch(rcspline.plot(p, y, model = "logistic",
                                 nk = 3, show = "prob", statloc = "none", add = T,
                                 showknots = F, xrange = c(min(na.omit(p)),
                                                           max(na.omit(p))), lty = lty.smooth), error = function(e) {
                                                             stop("Nk = 3 led to estimation problems.")
                                                           })
        }
        else {
          stop(paste("Number of knots = ", nr.knots, sep = "",
                     ", only 5 >= nk >=3 is allowed."))
        }
        par(lwd = 1, bty = "o", col = "black")
        leg <- c(leg, "Flexible calibration (RCS)", "CL flexible")
        lt <- c(lt, lty.smooth, 2)
        lw.d <- c(lw.d, rep(lwd.smooth, 2))
        all.col <- c(all.col, col.smooth)
        marks <- c(marks, -1, -1)
      }
      if (!missing(m) | !missing(g) | !missing(cuts)) {
        if (!missing(m))
          q <- cut2(p, m = m, levels.mean = T, digits = 7)
        else if (!missing(g))
          q <- cut2(p, g = g, levels.mean = T, digits = 7)
        else if (!missing(cuts))
          q <- cut2(p, cuts = cuts, levels.mean = T, digits = 7)
        means <- as.single(levels(q))
        prop <- tapply(y, q, function(x) mean(x, na.rm = T))
        points(means, prop, pch = 2, cex = 1)
        ng <- tapply(y, q, length)
        og <- tapply(y, q, sum)
        ob <- og/ng
        se.ob <- sqrt(ob * (1 - ob)/ng)
        g <- length(as.single(levels(q)))
        for (i in 1:g) lines(c(means[i], means[i]), c(prop[i],
                                                      min(1, prop[i] + 1.96 * se.ob[i])), type = "l")
        for (i in 1:g) lines(c(means[i], means[i]), c(prop[i],
                                                      max(0, prop[i] - 1.96 * se.ob[i])), type = "l")
        if (connect.group) {
          lines(means, prop)
          lt <- c(lt, 1)
          lw.d <- c(lw.d, 1)
        }
        else lt <- c(lt, 0)
        lw.d <- c(lw.d, 0)
        leg <- c(leg, "Grouped observations")
        marks <- c(marks, 2)
      }
    }
    lr <- stats["Model L.R."]
    p.lr <- stats["P"]
    D <- (lr - 1)/n
    L01 <- -2 * sum(y * logit - logb(1 + exp(logit)), na.rm = TRUE)
    U.chisq <- L01 - f$deviance[2]
    p.U <- 1 - pchisq(U.chisq, 2)
    U <- (U.chisq - 2)/n
    Q <- D - U
    Dxy <- stats["Dxy"]
    C <- stats["C"]
    R2 <- stats["R2"]
    B <- sum((p - y)^2)/n
    Bmax <- mean(y) * (1 - mean(y))^2 + (1 - mean(y)) * mean(y)^2
    Bscaled <- 1 - B/Bmax
    stats <- c(Dxy, C, R2, D, lr, p.lr, U, U.chisq, p.U, Q, B,
               f2$coef[1], f$coef[2], emax, Bscaled)
    names(stats) <- c("Dxy", "C (ROC)", "R2", "D", "D:Chi-sq",
                      "D:p", "U", "U:Chi-sq", "U:p", "Q", "Brier", "Intercept",
                      "Slope", "Emax", "Brier scaled")
    if (smooth == "loess")
      stats <- c(stats, c(Eavg = eavg), c(ECI = ECI))
    if (!missing(cutoff)) {
      arrows(x0 = cutoff, y0 = 0.1, x1 = cutoff, y1 = -0.025,
             length = 0.15)
    }
    if (pl) {
      if (min(p) > plogis(-7) | max(p) < plogis(7)) {
        lrm.fit.1 <- lrm(y[i.2] ~ qlogis(p[i.2]))
        if (logistic.cal)
          lines(p[i.2], plogis(lrm.fit.1$linear.predictors),
                lwd = lwd.log, lty = lty.log, col = col.log)
      }
      else {
        logit <- seq(-7, 7, length = 200)
        prob <- 1/(1 + exp(-logit))
        pred.prob <- f$coef[1] + f$coef[2] * logit
        pred.prob <- 1/(1 + exp(-pred.prob))
        if (logistic.cal)
          lines(prob, pred.prob, lty = lty.log, lwd = lwd.log,
                col = col.log)
      }
      lp <- legendloc
      if (!is.logical(lp)) {
        if (!is.list(lp))
          lp <- list(x = lp[1], y = lp[2])
        legend(lp, leg, lty = lt, pch = marks, cex = cex.leg,
               bty = "n", lwd = lw.d, col = all.col, y.intersp = y.intersp)
      }
      if (!is.logical(statloc)) {
        if (dostats[1] == T) {
          stats.2 <- paste("Calibration\n", "...intercept: ",
                           sprintf(paste("%.", roundstats, "f", sep = ""),
                                   stats["Intercept"]), " (", sprintf(paste("%.",
                                                                            roundstats, "f", sep = ""), cl.interc[1]),
                           " to ", sprintf(paste("%.", roundstats, "f",
                                                 sep = ""), cl.interc[2]), ")", "\n", "...slope: ",
                           sprintf(paste("%.", roundstats, "f", sep = ""),
                                   stats["Slope"]), " (", sprintf(paste("%.",
                                                                        roundstats, "f", sep = ""), cl.slope[1]),
                           " to ", sprintf(paste("%.", roundstats, "f",
                                                 sep = ""), cl.slope[2]), ")", "\n", "Discrimination\n",
                           "...c-statistic: ", sprintf(paste("%.", roundstats,
                                                             "f", sep = ""), stats["C (ROC)"]), " (",
                           sprintf(paste("%.", roundstats, "f", sep = ""),
                                   cl.auc[2]), " to ", sprintf(paste("%.", roundstats,
                                                                     "f", sep = ""), cl.auc[3]), ")", sep = "")
          text(statloc[1], statloc[2], stats.2, pos = 4,
               cex = cex)
        }
        else {
          dostats <- dostats
          leg <- format(names(stats)[dostats])
          leg <- paste(leg, ":", format(stats[dostats],
                                        digits = roundstats), sep = "")
          if (!is.list(statloc))
            statloc <- list(x = statloc[1], y = statloc[2])
          text(statloc, paste(format(names(stats[dostats])),
                              collapse = "\n"), adj = 0, cex = cex)
          text(statloc$x + (xlim[2] - xlim[1])/3, statloc$y,
               paste(format(round(stats[dostats], digits = roundstats)),
                     collapse = "\n"), adj = 1, cex = cex)
        }
      }
      if (is.character(riskdist)) {
        if (riskdist == "calibrated") {
          x <- f$coef[1] + f$coef[2] * log(p/(1 - p))
          x <- 1/(1 + exp(-x))
          x[p == 0] <- 0
          x[p == 1] <- 1
        }
        else x <- p
        bins <- seq(0, min(1, max(xlim)), length = 101)
        x <- x[x >= 0 & x <= 1]
        f0 <- table(cut(x[y == 0], bins))
        f1 <- table(cut(x[y == 1], bins))
        j0 <- f0 > 0
        j1 <- f1 > 0
        bins0 <- (bins[-101])[j0]
        bins1 <- (bins[-101])[j1]
        f0 <- f0[j0]
        f1 <- f1[j1]
        maxf <- max(f0, f1)
        f0 <- (0.1 * f0)/maxf
        f1 <- (0.1 * f1)/maxf
        segments(bins1, line.bins, bins1, length.seg * f1 +
                   line.bins)
        segments(bins0, line.bins, bins0, length.seg * f0 +
                   line.bins)
        lines(c(min(bins0, bins1) - 0.01, max(bins0, bins1) +
                  0.01), c(line.bins, line.bins))
        text(max(bins0, bins1) + dist.label, line.bins +
               dist.label2, d1lab, cex = cex.d01)
        text(max(bins0, bins1) + dist.label, line.bins -
               dist.label2, d0lab, cex = cex.d01)
      }
    }
    if (dostats == T) {
      cat(paste("\n\n A ", cl.level * 100, "% confidence interval is given for the calibration intercept, calibration slope and c-statistic. \n\n",
                sep = ""))
    }
    stats
  }

#' My plot.xmean.ordinaly from rms
#' removed the xlab and ylab defaults from plot.xmean.ordinaly - so I could set my own
My.plot.xmean.ordinaly =
  function (x, data, subset, na.action, subn = TRUE, cr = FALSE, 
            topcats = 1, cex.points = 0.75, ...) 
  {
    X <- match.call(expand.dots = FALSE)
    X$subn <- X$cr <- X$topcats <- X$cex.points <- X$... <- NULL
    if (missing(na.action)) 
      X$na.action <- na.keep
    Terms <- if (missing(data)) 
      terms(x)
    else terms(x, data = data)
    X$formula <- Terms
    X[[1]] <- as.name("model.frame")
    X <- eval.parent(X)
    resp <- attr(Terms, "response")
    if (resp == 0) 
      stop("must have a response variable")
    nx <- ncol(X) - 1
    Y <- X[[resp]]
    nam <- as.character(attr(Terms, "variables"))
    nam <- nam[-1]
    dopl <- function(x, y, cr, xname, yname) {
      s <- !is.na(unclass(Y) + x)
      y <- y[s]
      x <- x[s]
      n <- length(x)
      f <- lrm.fit(x, y)
      fy <- f$freq/n
      ns <- length(fy) - 1
      k <- ns + 1
      intcept <- f$coef[1:ns]
      xb <- f$linear.predictors - intcept[1]
      xb <- sapply(intcept, "+", xb)
      P <- 1/(1 + exp(-xb))
      P <- matrix(P, ncol = ns)
      P <- cbind(1, P) - cbind(P, 0)
      xmean.y <- tapply(x, y, mean)
      xp <- x * P/n
      xmean.y.po <- apply(xp, 2, sum)/fy
      yy <- 1:length(fy)
      rr <- c(xmean.y, xmean.y.po)
      if (cr) {
        u <- cr.setup(y)
        s <- u$subs
        yc <- u$y
        xc <- x[s]
        cohort <- u$cohort
        xcohort <- matrix(0, nrow = length(xc), ncol = length(levels(cohort)) - 
                            1)
        xcohort[col(xcohort) == unclass(cohort) - 1] <- 1
        cof <- lrm.fit(cbind(xcohort, xc), yc)$coefficients
        cumprob <- rep(1, n)
        for (j in 1:k) {
          P[, j] <- cumprob * (if (j == k) 
            1
            else plogis(cof[1] + (if (j > 1) 
              cof[j]
              else 0) + cof[k] * x))
          cumprob <- cumprob - P[, j]
        }
        xp <- x * P/n
        xmean.y.cr <- apply(xp, 2, sum)/fy
        rr <- c(rr, xmean.y.cr)
      }
      plot(yy, xmean.y, type = "b", ylim = range(rr), axes = FALSE, 
           ...)
      mgp.axis(1, at = yy, labels = names(fy))
      mgp.axis(2)
      lines(yy, xmean.y.po, lty = 2, ...)
      if (cr) 
        points(yy, xmean.y.cr, pch = "C", cex = cex.points)
      if (subn) 
        title(sub = paste("n=", n, sep = ""), adj = 0)
    }
    for (i in 1:nx) {
      x <- X[[resp + i]]
      if (is.factor(x)) {
        f <- table(x)
        ncat <- length(f)
        if (ncat < 2) {
          warning(paste("predictor", nam[resp + i], "only has one level and is ignored"))
          next
        }
        nc <- min(ncat - 1, topcats)
        cats <- (names(f)[order(-f)])[1:nc]
        for (wcat in cats) {
          xx <- 1 * (x == wcat)
          xname <- paste(nam[resp + i], wcat, sep = "=")
          dopl(xx, Y, cr, xname, nam[resp])
        }
      }
      else dopl(x, Y, cr, nam[resp + i], nam[resp])
    }
    invisible()
  }



