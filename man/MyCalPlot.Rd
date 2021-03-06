% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{MyCalPlot}
\alias{MyCalPlot}
\title{MyCalPlot function
my hack of validate.plot.default in rms for my 2 panel plots to match with val.prob}
\usage{
MyCalPlot(p, y, logit, group, weights = rep(1, length(y)),
  normwt = FALSE, pl = TRUE, smooth = TRUE, logistic.cal = TRUE,
  xlab = "Predicted Probability", ylab = "Actual Probability",
  lim = c(0, 1), m, g, cuts, emax.lim = c(0, 1), legendloc = lim[1] +
  c(0.55 * diff(lim), 0.27 * diff(lim)), statloc = c(0, 0.99),
  riskdist = "calibrated", cex = 0.7, mkh = 0.02,
  connect.group = FALSE, connect.smooth = TRUE, g.group = 4,
  evaluate = 100, nmin = 0)
}
\description{
MyCalPlot function
my hack of validate.plot.default in rms for my 2 panel plots to match with val.prob
}
\examples{
MyCalPlot()
}
\keyword{Calibration}
\keyword{curves}
