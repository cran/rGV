#' Plot glucose values over time
#'
#' @param x vector of glucose readings
#' @param times vector of corresponding times, in minutes
#' @param unit "mg" if the units are mg/dL or "mmol" if the units are mmol/L. Null value is "mg".
#' @return A plot of glucose values over time.
#' @export
#'
#' @examples
#' cgm_plot(x=c(rep(100, 10), rep(120, 10), 105, 85), times=seq(0, 1260, 60), unit='mg')
cgm_plot <- function(x, times, unit="mg") {
  times <- times / (60 * 24)
  graphics::plot(x ~ times, type="l",
       ylim=c(0, max(x) + 50), xlab="Time, in days",
       ylab=paste0("Blood glucose, in ", ifelse(unit=="mg", "mg/dL", "mmol/L")))
}

#' Plots glucose changes over time
#'
#' @param x vector of glucose readings
#' @param times vector of corresponding times, in minutes
#' @param n number of hours between "partner" observations. Null value is 1.
#' @param s number of minutes of slack used when searching for partners. Null value is 1.
#' @param unit "mg" if the units are mg/dL or "mmol" if the units are mmol/L. Null value is "mg".
#'
#' @return A plot of n-hour glucose differences over time.
#'
#' @export
#'
#' @examples
#' diff_plot(x=c(rep(100, 10), rep(120, 10), 105, 85), times=seq(0, 1260, 60), n=1, s=1, unit='mg')
diff_plot <- function(x, times, n=1, s=1, unit="mg") {
  dt <- vector() # vector to store differences
  for (ti in times) {
    index <- which(times == ti)
    # identify partner(s) for time t
    partner <- which(times < ti - n * 60 + s & times > ti - n * 60 - s)
    if (length(partner) > 1)
      stop('Multiple glucose value pairs found for a single value. Reduce slack.')
    if (length(partner) == 1) dt[index] <-  x[index] - x[partner] # calculate differences
    if (length(partner) == 0) dt[index] <-  NA # calculate differences
  }

  times <- times / (60 * 24)
  graphics::plot(dt ~ times, type="l", xlab="Time, in days",
       ylab=paste0(n, "-hour difference in glucose values, in ",
                   ifelse(unit=="mg", "mg/dL", "mmol/L")))
  graphics::abline(h=mean(dt, na.rm=TRUE), col="red")
  graphics::legend("topleft", legend=c(paste0("Mean ", n, "-hour difference")),
         col=c("red"), lty=1, cex=0.5)
}

#' Plot the symmetrized glucose values
#'
#' @param x vector of glucose readings
#' @param times vector of corresponding times, in minutes
#' @param unit "mg" if the units are mg/dL or "mmol" if the units are mmol/L. Null value is "mg".
#' @return A plot of symmetrized glucose values over time. These symmetrized values are used in the calculation of BGI and ADRR.
#' @export
#'
#' @examples
#' symm_plot(x=c(rep(100, 10),rep(120, 10), 105, 85), times=seq(0, 1260, 60), unit='mg')
symm_plot <- function(x, times, unit="mg") {
  if (unit == "mg") f_vec <- 1.509 * ( log(x)^1.084 - 5.381) # calculate f(x)
  if (unit == "mmol") f_vec <- 1.509 * ( log(x*18)^1.084 - 5.381) # calculate f(x)

  times <- times / (60 * 24)
  graphics::plot(f_vec ~ times, type="l", xlab="Time, in days",
       ylab=paste0("Symmetrized Glucose Values"), ylim=c(-max(min(f_vec), max(f_vec)),
                                                         max(min(f_vec), max(f_vec))))
  graphics::abline(h=0, col="red")
}
