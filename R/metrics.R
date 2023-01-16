#' Calculate continuous overall net glycemic action (CONGA)
#'
#' @param x vector of glucose readings
#' @param times vector of corresponding times, in minutes
#' @param n number of hours between "partner" observations. Null value is 1.
#' @param s number of minutes of slack used when searching for partners. Null value is 1.
#' @param method "manuscript" or "easy". Null value is "manuscript".
#'
#' @return The numeric CONGA value for a given dataset of glucose measurements and times.
#'
#' @export
#'
#' @examples
#' conga(x=c(rep(100, 10), rep(120, 10), 105, 85), times=seq(0, 1260, 60),
#' n=1, s=1, method="manuscript")
conga <- function(x, times, n=1, s=1, method="manuscript") {
  dt <- vector() # vector to store differences
  for (ti in times) {
    if (is.na(ti)) next
    index <- which(times == ti)
    partner <- which(times < ti - n * 60 + s & times > ti - n * 60 - s) # identify partner(s) for time t
    if (length(partner) > 1) dt[index] <- x[index] - mean(x[partner])
    if (length(partner) == 1) dt[index] <-  x[index] - x[partner] # calculate differences
    if (length(partner) == 0) dt[index] <-  NA # calculate differences
  }

  if (method == "manuscript") value <- stats::sd(dt, na.rm=TRUE)
  if (method == "easy") {
    keep <- which(!is.na(dt))
    value <- sqrt( sum((x[keep] - mean(abs(dt), na.rm=TRUE))^2)/ (length(keep) - 1) )
  }

  value
}

#' Calculate the lability index (LI)
#'
#' @param x vector of glucose readings
#' @param times vector of corresponding times, in minutes
#' @param k length of time (in minutes) used to find partners. Null value is 60.
#' @param s number of minutes of slack used when searching for partners. Null value is 1.
#'
#' @return The numeric value of the lability index for a given dataset of glucose measurements and times.
#' @export
#'
#' @examples
#' li(x=c(rep(100, 10), rep(120, 10), 105, 85), times=seq(0, 1260, 60), k=60, s=1)
li <- function(x, times, k=60, s=1) {
  num <- vector()
  denom <- vector()
  for (ti in times) {
    if (is.na(ti)) next
    index <- which(times == ti)
    partner <- which(times < ti - k + s & times > ti - k - s) # identify partner(s) for time t
    if (length(partner) > 1) {
      num[index] <-  (x[index] - mean(x[partner]))^2
      denom[index] <- (times[index] - mean(times[partner])) / 60 # calculate time difference, in hours
    }
    if (length(partner) == 0) num[index] <- NA; denom[index] <- NA
    if (length(partner) == 1) {
      num[index] <-  (x[index] - x[partner])^2
      denom[index] <- (times[index] - times[partner]) / 60 # calculate time difference, in hours
    }
  }
  value <- sum(num, na.rm=TRUE) / sum(denom, na.rm=TRUE)

  value
}

#' Calculate J-index
#'
#' @param x vector of glucose readings
#' @param unit "mg" if the units are mg/dL or "mmol" if the units are mmol/L. Null value is "mg".
#'
#' @return The numeric J-index value for a given dataset of glucose measurements.
#' @export
#'
#' @examples
#' j_index(x=c(rep(100, 10), rep(120, 10), 105, 85), unit='mg')
j_index <- function(x, unit="mg") {
  if (unit == "mg") value <- 0.001 * (mean(x, na.rm=TRUE) + stats::sd(x, na.rm=TRUE))^2
  if (unit == "mmol") value <- 0.324 * (mean(x, na.rm=TRUE) + stats::sd(x, na.rm=TRUE))^2
  value
}

#' Calculate Low / High Blood Glucose Index (LBGI, HBGI)
#'
#' @param x vector of glucose readings
#' @param unit "mg" if the units are mg/dL or "mmol" if the units are mmol/L. Null value is "mg".
#' @param method "manuscript", "corrected", or "easy". Null value is "manuscript".
#' @return A list containing the LBGI and HBGI values for a given dataset of glucose measurements.
#' @export
#'
#' @examples
#' bgi(x=c(rep(100, 10), rep(120, 10), 105, 85), unit='mg', method='manuscript')
bgi <- function(x, unit="mg", method="manuscript") {
  if (unit == "mg") f_vec <- 1.509 * ( log(x)^1.084 - 5.381 ) # calculate f(x)
  if (unit == "mmol") f_vec <- 1.509 * ( log(x*18)^1.084 - 5.381 ) # calculate f(x)
  rl <- ifelse(f_vec < 0, 10 * f_vec^2, 0)
  rh <- ifelse(f_vec > 0, 10 * f_vec^2, 0)

  if (method == "manuscript") {
    lbgi <- mean(rl, na.rm=TRUE)
    hbgi <- mean(rh, na.rm=TRUE)
  }

  if (method == "easy") {
    lbgi <- mean(rl[which(rl > 0)], na.rm=TRUE)
    hbgi <- mean(rh[which(rh > 0)], na.rm=TRUE)
  }

  if (method == "corrected") {
    rl <- ifelse(f_vec < 0, 22.77 * f_vec^2, 0)
    rh <- ifelse(f_vec > 0, 22.77 * f_vec^2, 0)
    lbgi <- mean(rl, na.rm=TRUE)
    hbgi <- mean(rh, na.rm=TRUE)
  }

  list("lbgi"=lbgi, "hbgi"=hbgi)
}

#' Calculate Glycemic Risk Assessment Diabetes Equation (GRADE)
#'
#' @param x vector of glucose readings
#' @param unit "mg" if the units are mg/dL or "mmol" if the units are mmol/L. Null value is "mg".
#' @param method "manuscript" or "easy". Null value is "manuscript".
#' @param c1 glucose value below which readings are considered hypoglycemic. Default is 70.2 mg/dL.
#' @param c2 glucose value above which readings are considered hyperglycemic. Default is 140.4 mg/dL.
#' @return A list containing the GRADE value and the percentage of the GRADE value due to euglycemia, hypoglycemia, and hyperglycemia for a given dataset of glucose measurements
#' @export
#'
#' @examples
#' grade(x=c(rep(100, 10), rep(120, 10), 105, 85), unit='mg', method='manuscript', c1=70.2, c2=140.4)
grade <- function(x, unit="mg", method="manuscript", c1=ifelse(unit=="mg", 70.2, 3.9),
                  c2=ifelse(unit=="mg", 140.4, 7.8)) {
  if (c2 < c1) stop('Incorrect c1 and c2 bounds entered')
  # Calculate g(x)
  if (unit == "mg" & method == "manuscript") vec <- pmin(425 * (log(log(x / 18, 10), 10) + 0.16)^2, 50)
  if (unit == "mmol" & method == "manuscript") vec <- pmin(425 * (log(log(x, 10), 10) + 0.16)^2, 50)
  if (unit == "mg" & method == "easy") vec <- pmin(425 * (log(log(x / 18, 10), 10) + 0.15554147)^2, 50)
  if (unit == "mmol" & method == "easy") vec <- pmin(425 * (log(log(x, 10), 10) + 0.15554147)^2, 50)

  if (method == "manuscript") value <- mean(vec, na.rm=TRUE)
  if (method == "easy") value <- stats::median(vec, na.rm=TRUE)

  # Calculate percentage hypo, hyper and eu
  if (unit == "mg") hypo <- sum(vec[which(x < c1)], na.rm=TRUE) / sum(vec, na.rm=TRUE)
  if (unit == "mg") hyper <- sum(vec[which(x > c2)], na.rm=TRUE) / sum(vec, na.rm=TRUE)

  if (unit == "mmol") hypo <- sum(vec[which(x < c1)], na.rm=TRUE) / sum(vec, na.rm=TRUE)
  if (unit == "mmol") hyper <- sum(vec[which(x > c2)], na.rm=TRUE) / sum(vec, na.rm=TRUE)

  eu <- 1 - hypo - hyper
  list("grade"=value, "hypo"=hypo, "eu"=eu, "hyper"=hyper)
}

#' Calculate Mean of Daily Differences (MODD)
#'
#' @param x vector of glucose readings
#' @param times vector of corresponding times, in minutes
#' @param s number of minutes of slack used when searching for partners. Null value is 1.
#' @param method "manuscript" or "easy". Null value is "manuscript".
#' @return The numeric MODD value for a given dataset of glucose measurements and times.
#' @export
#'
#' @examples
#' modd(x=c(rep(100, 10), rep(120, 10), 105, 85), times=seq(0, 1260, 60), s=1, method='manuscript')
modd <- function(x, times, s=1, method="manuscript") {
  comp <- vector() # vector of differences
  for (i in 1:length(x)) {
    cand <- x[which(abs(times + 24 * 60 - times[i]) < s)] # finds partner(s) for x
    comp[i] <- ifelse(length(cand) > 0, mean(cand), NA)
  }

  if (method=="easy") value <- mean(utils::head(abs(x - comp), -1), na.rm=TRUE)
  if (method=="manuscript") value <- mean(abs(x - comp), na.rm=TRUE)
  value
}

#' Calculate Mean Amplitude of Glycemic Excursions (MAGE)
#'
#' @param x vector of glucose readings
#' @param times vector of corresponding times, in minutes
#' @return The numeric MAGE value for a given dataset of glucose measurements and times.
#' @export
#'
#' @examples
#' mage(x=c(rep(100, 10), rep(120, 10), 105, 85), times=seq(0, 1260, 60))
mage <- function(x, times) {
  x=c(rep(100,10),rep(120,10),105,85)
  times=seq(0,1260,60)

  smoothed_values <- vector()
  smoothed_values[1:4] <- mean(stats::na.omit(x[1:4]))
  for (i in 5:(length(x)-4)){
    smoothed_values[i] <- (1*x[i - 4] + 2*x[i - 3] + 4*x[i - 2] + 8*x[i - 1] + 16*x[i] +
                             8*x[i + 1] + 4*x[i + 2] + 2*x[i + 3] + 1*x[i + 4]) / 46
  }
  ti <- length(smoothed_values)
  smoothed_values[(ti + 1):(ti + 4)] <- mean(x[(length(x) - 3):length(x)])
  sd <- stats::sd(x)

  u <- unique(smoothed_values)
  n <- length(u)
  v1 <- c(0,u[2:n] - u[1:(n - 1)])
  v2 <- c(u[2:n] - u[1:(n - 1)],0)

  peaks <- which(v1 > 0 & v2 < 0)
  troughs <- which(v1 < 0 & v2 > 0)
  if (length(peaks) > length(troughs)) troughs <- append(troughs, 1)
  if (length(peaks) < length(troughs)) peaks <- append(peaks, 1)

  differences <- u[peaks] - u[troughs]
  mean(stats::na.omit(differences[which(differences > sd)]))
}

### Measure 10: Average Daily Risk Range (ADRR)
# This function takes a vector of glucose values and a vector of times and calculates the ADRR
#' Calculate Average Daily Risk Range (ADRR)
#'
#' @param x vector of glucose readings
#' @param times vector of corresponding times, in minutes
#' @param unit "mg" if the units are mg/dL or "mmol" if the units are mmol/L. Null value is "mg".
#' @param method "manuscript", "corrected", or "easy". Null value is "manuscript".
#' @return The numeric ADRR value for a given dataset of glucose measurements and times.
#' @export
#'
#' @examples
#' adrr(x=c(rep(100, 10), rep(120, 10), 105, 85), times=seq(0, 1260, 60),
#' unit="mg", method='manuscript')
adrr <- function(x, times, unit="mg", method="manuscript") {
  if (unit == "mg") f_vec <- 1.509 * ( (log(x))^1.084 - 5.381 ) # calculate f(x)
  if (unit == "mmol") f_vec <- 1.509 * ( (log(x*18))^1.084 - 5.381 ) # calculate f(x)
  rl <- ifelse(f_vec < 0, 10 * f_vec^2, 0)
  rh <- ifelse(f_vec > 0, 10 * f_vec^2, 0)

  day <- ceiling(times / (24 * 60) + 0.005) # finds the day when each drawing occurred
  combined <- as.data.frame(cbind(rl, rh, day))
  lrd <- stats::aggregate(rl~day, data=combined, function(x) max=max(x, na.rm=TRUE))$rl # finds max rl by day
  hrd <- stats::aggregate(rh~day, data=combined, function(x) max=max(x, na.rm=TRUE))$rh # finds max rh by day
  if (method == "manuscript") value <- mean(lrd) + mean(hrd)
  if (method == "easy") value <- c(mean(lrd), mean(hrd))

  value
}

#' Calculate M-value
#'
#' @param x vector of glucose readings
#' @param unit "mg" if the units are mg/dL or "mmol" if the units are mmol/L. Null value is "mg".
#' @param index value to be considered a 'standard' blood glucose value, in mg/dL. Null value is 120.
#' @param method "manuscript", "corrected", or "easy". Null value is "manuscript".
#' @return The numeric M-value for a given dataset of glucose measurements.
#' @export
#'
#' @examples
#' m_value(x=c(rep(100, 10), rep(120, 10), 105, 85), unit='mg', index=120, method='manuscript')
m_value <- function(x, unit="mg", index=120, method="manuscript") {
  if (unit == "mg") mbs <- mean(abs(10 * log(x/index, 10))^3, na.rm=TRUE)
  if (unit == "mg") mw <- (max(x, na.rm=TRUE) - min(x, na.rm=TRUE)) / 20

  if (unit == "mmol") mbs <- mean(abs(10 * log(x * 18/index, 10))^3, na.rm=TRUE)
  if (unit == "mmol") mw <- (max(x * 18, na.rm=TRUE) - min(x * 18, na.rm=TRUE)) / 20

  if (method == "manuscript") value <- mbs + mw
  if (method == "easy") value <- mbs

  value
}

#' Calculate Mean Absolute Glucose (MAG)
#'
#' @param x vector of glucose readings
#' @param times vector of corresponding times, in minutes
#' @return The numeric MAG value for a given dataset of glucose measurements and times.
#' @export
#'
#' @examples
#' mag(x=c(rep(100, 10),rep(120, 10), 105, 85), times=seq(0, 1260, 60))
mag <- function(x, times) {
  value <- sum(abs(x[2:length(x)] - x[1:(length(x) - 1)]), na.rm=TRUE) /
    as.numeric((max(times, na.rm=TRUE) - min(times, na.rm=TRUE))/60)

  value
}

#' Calculate coefficient of variation (CV)
#'
#' @param x vector of glucose readings
#' @param times vector of corresponding times, in minutes
#' @param overall a logical, equal to TRUE you want the CV for the entire dataset, or equal to FALSE if you would prefer many CV values over a moving window
#' @param interval size (in hours) of the moving window to be used if overall is false. Null value is 1.
#' @return Either a numeric coefficient of variation over the entire dataset or a vector of CV values over windows of the data.
#' @export
#'
#' @examples
#' cv(x=c(rep(100, 10), rep(120, 10), 105, 85), times=seq(0, 1260, 60), overall=TRUE)
cv <- function(x, times, overall=TRUE, interval=1) {
  if (overall == TRUE) {
    value <- stats::sd(x, na.rm=TRUE) / mean(x, na.rm=TRUE)
    return(value)
  }

  if (overall == FALSE) {
    vec <- vector()
    for(ti in times[1:(length(times) - 2)]) {
      if (is.na(ti)) next
      index <- which(times == ti)
      time_int <- which(times >= ti & times <= interval * 60 + ti)
      if (length(time_int) < 2) stop("Not enough data within time interval")
      vec[index] <- stats::sd(x[time_int], na.rm=TRUE) / mean(x[time_int], na.rm=TRUE)
    }
    return(vec)
  }
}

#' Calculate standard deviation (SD)
#'
#' @param x vector of glucose readings
#' @param times vector of corresponding times, in minutes
#' @param overall a logical, equal to TRUE you want the CV for the entire dataset, or equal to FALSE if you would prefer many CV values over a moving window
#' @param interval size (in hours) of the moving window to be used if overall is false. Null value is 1.
#' @return Either a numeric standard deviation over the entire dataset or a vector of SD values over windows of the data.
#' @export
#'
#' @examples
#' st_dev(x=c(rep(100, 10), rep(120, 10), 105, 85), times=seq(0, 1260, 60), overall=TRUE)
st_dev <- function(x, times, overall=TRUE, interval=1) {
  if (overall == TRUE) {
    value <- stats::sd(x, na.rm=TRUE)
    return(value)
  }

  if (overall == FALSE) {
    vec <- vector()
    for(ti in times[1:(length(times) - 2)]) {
      if (is.na(ti)) next
      index <- which(times == ti)
      time_int <- which(times >= ti & times <= interval * 60 + ti)
      if (length(time_int) < 2) stop("Not enough data within time interval")
      vec[index] <- stats::sd(x[time_int], na.rm=TRUE)
    }
    return(vec)
  }
}

#' Calculate area under the curve (AUC)
#'
#' @param x vector of glucose readings
#' @param times vector of corresponding times, in minutes
#' @param thresh threshold above (or below) which you wish to calculate the AUC. Default is 100.
#' @param above logical indicating whether you wish to calculate area above the threshold value (TRUE) or below it (FALSE). Default is TRUE.
#' @return The numeric area under the curve value for a given dataset of glucose measurements and times.
#' @export
#'
#' @examples
#' cgm_auc(x=c(rep(100, 10), rep(120, 10), 105, 85), times=seq(0, 1260, 60), thresh=110, above=TRUE)
cgm_auc <- function(x, times, thresh=100, above=TRUE) {
  x_t <- x - thresh
  value <- 0
  for(index in 1:(length(x_t) - 1)) {
    x1 <- x_t[index]
    x2 <- x_t[index + 1]
    if (is.na(x1) | is.na(x2)) next
    ti <- times[index + 1] - times[index]
    if (is.na(ti)) next
    if (above == TRUE & x1 >= 0 & x2 >= 0) {
      value <- value + min(x1, x2) * ti + abs(x2 - x1) * ti / 2
    }
    if (above == FALSE & x1 <= 0 & x2 <= 0) {
      value <- value + abs(min(x1, x2) * ti + abs(x2 - x1) * ti / 2)
    }
  }
  value <- as.numeric(value) / ((max(times) - min(times)) / (60 * 24))

  value
}

#' Calculate time in range (TIR)
#'
#' @param x vector of glucose readings
#' @param low lower bound of the range. Default is 70
#' @param high upper bound of the range. Default is 180
#' @return The numeric time in range value for a given dataset of glucose measurements and times.
#' @export
#'
#' @examples
#' tir(x=c(rep(100, 10), rep(120, 10), 105, 85), low=70, high=80)
tir <- function(x, low=70, high=180) {
  value <- 100*length(x[which(x >= low & x <= high)]) / sum(!is.na(x))
  value
}

#' Calculate Glucose Management Indicator (GMI)
#'
#' @param x vector of glucose readings
#' @param unit "mg" if the units are mg/dL or "mmol" if the units are mmol/L. Null value is "mg".
#' @return The numeric GMI value for a given dataset of glucose measurements.
#' @export
#'
#' @examples
#' gmi(x=c(rep(100, 10), rep(120, 10), 105, 85), unit='mg')
gmi <- function(x, unit="mg") {
  if (unit == "mg") value <- 3.31 + 0.02392 * mean(x, na.rm=TRUE)
  if (unit == "mmol") value <- 12.71 + 4.70587 * mean(x, na.rm=TRUE)
  value
}

#' Find number of episodes below a given glucose value for a given amount of time
#'
#' @param x vector of glucose readings
#' @param times vector of corresponding times, in minutes
#' @param thresh glucoses below this threshold are considered as part of an episode. Default is 55
#' @param len minimum length of an episode. Default is 15
#' @param gap typical gap between CGM measurements, in minutes. Default is 5
#' @return The integer number of events for a given dataset of glucose measurements and times.
#' @export
#'
#' @examples
#' num_events(x=c(rep(100, 10), rep(120, 10), 105, 85), times=seq(0, 1260, 60),
#' thresh=55, len=15, gap=5)
num_events <- function(x, times, thresh=55, len=15, gap=5) {
  switch <- 0
  time1 <- NULL
  time2 <- NULL
  for (i in 1:length(x)) {
    if (switch == 0 & x[i] < thresh) {
      switch <- 1
      time1 <- c(time1, times[i])
    }
    if (switch == 1 & x[i] > thresh) {
      switch <- 0
      time2 <- c(time2, times[i-1])
    }
  }
  if (length(time2) < length(time1)) {
    time2 <- append(time2, max(times))
  }
  total_time <- time_on(times, gap=gap) / (60 * 24)
  length(which(time2 - time1 > len)) / total_time
}

#' Calculate Glycemic Variability Percentage (GVP)
#'
#' @param x vector of glucose readings
#' @param times vector of corresponding times, in minutes
#' @return The numeric GVP value for a given dataset of glucose measurements and times.
#' @export
#'
#' @examples
#' gvp(x=c(rep(100, 10), rep(120, 10), 105, 85), times=seq(0, 1260, 60))
gvp <- function(x, times) {
  n <- length(x)
  delta_gluc <- x[2:n] - x[1:(n - 1)]
  delta_time <- times[2:n] - times[1:(n - 1)]

  L <- sum(sqrt(delta_gluc^2 + delta_time^2))
  L0 <- sum(delta_time)
  value <- (L / L0 - 1) * 100
  value
}

#' Calculate distance travelled
#'
#' @param x vector of glucose readings
#' @return The numeric distance travelled value for a given dataset of glucose measurements.
#' @export
#'
#' @examples
#' dist_travelled(x=c(rep(100, 10), rep(120, 10), 105, 85))
dist_travelled <- function(x) {
  n <- length(x)
  value <- sum(abs(x[2:n] - x[1:(n - 1)]))
  value
}

#' Calculate amount of time that the CGM was active
#'
#' @param times vector of corresponding times, in minutes
#' @param gap typical gap between CGM measurements, in minutes. Default is 5
#' @return The numeric amount of time that the CGM device was active in a given dataset.
#' @export
#'
#' @examples
#' time_on(times=seq(0, 1260, 60), gap=5)
time_on <- function(times, gap=5) {
  total_time <- max(times) - min(times)
  gaps <- times[2:length(times)] - times[1:(length(times) - 1)]
  total_time <- total_time - sum(gaps[which(gaps > gap + 2)]) + gap*length(which(gaps > gap + 2))
  total_time
}
