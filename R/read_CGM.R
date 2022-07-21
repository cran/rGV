#' Read in continuous glucose monitor data
#'
#' @param file name of the CSV file to be read in
#' @param timezero set to "first" if the first glucose reading should be considered time zero and set to "midnight" if midnight of the day of the first reading should be considered time zero. Default is "first".
#' @param na.rm a logical that is TRUE if you wish to exclude all readings that are missing glucose values or time stamps and FALSE if not. Default is TRUE.
#' @param skip the number of lines in the data file to skip before beginning to read in data
#' @param calib_col the number or name of the column containing information regarding calibration status of each glucose entry
#' @param calib_tag the character value used to denote calibration rows in calib_col
#' @param mult_sensors a logical that is TRUE if you wish to split the data set into parts corresponding to different CGM sensors and FALSE if not. Default is FALSE.
#' @param sensor_times a vector of times (in the same format as the time data) that correspond to the beginning of a new CGM sensor. These times are used to split the data between multiple sensors if mult_sensors is TRUE. If sensor_times is NA, the data is split automatically at every gap of sensor_gap or more minutes.
#' @param sensor_gap a number specifying the minimum gap (in minutes) for which we should split the data into two pieces. Default is 120.
#' @param time_col the number or name of the column containing time data
#' @param gluc_col the number or name of the column containing glucose data
#' @param time_sep character that separates date from time in your time data
#' @param time_format specify date and time formats according to the specification used in the chron package. Default is c(dates = "m/d/y", times = "h:m:s").
#' @param high_ind character value that identifies high glucose values in the data. Default is "High".
#' @param high_value numeric value by which to replace glucose values equal to "high_ind". Default is 400.
#' @param low_ind character value that identifies low glucose values in the data. Default is "Low".
#' @param low_value numeric value by which to replace glucose values equal to "low_ind". Default is 40.
#'
#' @import chron
#'
#' @return A data frame with two columns: glucose values and time. This data frame can then be used with other rGV functions to calculate CGM metrics.
#' @export
read_cgm <- function(file, timezero="first", na.rm=TRUE, skip=0,
                     calib_col=NA, calib_tag,
                     mult_sensors=FALSE, sensor_times=NA, sensor_gap=120,
                     time_col, gluc_col, time_sep=" ",
                     time_format = c(dates = "m/d/y", times = "h:m:s"),
                     high_ind="High", high_value=400, low_ind="Low", low_value=40) {
  # file: specify path to the CSV file

  data <- utils::read.csv(file)

  # Skip a specified number of rows at the beginning of the data file
  if (skip > nrow(data)) stop("'skip' must be less than the number of rows in the data file")
  if (is.numeric(skip)) data <- data[(skip+1):nrow(data), ] # eliminate header rows
  if (!is.numeric(skip)) stop("Enter a numeric value for 'skip'")

  # Delete rows which are specified as calibration rows
  if (is.numeric(calib_col) & calib_col > ncol(data)) stop("'calib_col' must specify a column in the data file")
  if (is.character(calib_col) & !(calib_col %in% names(data))) stop("'calib_col' must specify a column in the data file")
  if (!is.na(calib_col)) data <- data[!data[, calib_col] == calib_tag, ] # eliminate rows where calibration occurred

  # Keep only the time and glucose columns specified
  if (is.numeric(time_col) & time_col > ncol(data)) stop("'time_col' must specify a column in the data file")
  if (is.numeric(gluc_col) & gluc_col > ncol(data)) stop("'gluc_col' must specify a column in the data file")
  if (is.character(time_col) & !(time_col %in% names(data))) stop("'time_col' must specify a column in the data file")
  if (is.character(gluc_col) & !(gluc_col %in% names(data))) stop("'gluc_col' must specify a column in the data file")
  data <- data[, c(time_col, gluc_col)] # keep only timestamp and glucose variables
  names(data) <- c("timestamp", "glucose")

  # Replace high and low values
  data$glucose <- as.character(data$glucose)
  data$glucose <- as.numeric(ifelse(data$glucose == high_ind, high_value,
                                    ifelse(data$glucose == low_ind, low_value, data$glucose)))

  # Split timestamp into date and time
  if (all(grepl(time_sep, data$timestamp))) {
    data$date <- gsub(paste0(time_sep, ".*$"), "", data$timestamp)
    data$time <- gsub(paste0(".*", time_sep), "", data$timestamp)
  } else stop("'time_sep' must appear in every time value in your data")

  # Remove rows without a timestamp
  data <- data[which(!is.na(data$timestamp)), ]

  # Add seconds to the data if it did not include seconds
  if (grepl("s", time_format[2]) == FALSE) {
    data$time <- paste0(data$time, ":00")
    time_format[2] <- paste0(time_format[2], ":s")
    switch <- 1
  }

  data$datetime <- chron::chron(dates=data$date, times=data$time, format=time_format) # interpretable date and time

  # Measure time, in minutes, since the first glucose measurement
  if (timezero == "first") data$elapsed <- as.numeric(difftime(data$datetime, data$datetime[1], units='mins'))
  # Measure time, in minutes, since midnight of first day of measurement
  if (timezero == "midnight") data$elapsed <- as.numeric(difftime(data$datetime, chron::chron(data$date[1], 00:00:00, format=time_format), units='mins'))
  # Measure time, in minutes, since user-specified time zero
  if (timezero != "first" & timezero != "midnight") {
    timezero.date <- gsub(paste0(time_sep, ".*$"), "", timezero)
    timezero.time <- gsub(paste0(".*", time_sep), "", timezero)
    if (switch <- 1) {
      timezero.time <- paste0(timezero.time, ":00")
    }
    timezero <- chron::chron(dates=timezero.date, times=timezero.time, format=time_format)
    data$elapsed <- as.numeric(difftime(data$datetime, timezero, units='mins'))
    # stop("Please specify a timezero value of 'first', 'midnight', or a specific time stamp")
  }

  data <- data[order(data$elapsed), ]
  data <- data[!duplicated(data$elapsed), ]

  if (mult_sensors == TRUE && !is.na(sensor_times)) {
    # split based on user-specified times
    # change sensor_times vector to chron format
    date_of_sensor <- gsub(paste0(time_sep, ".*$"), "", sensor_times)
    time_of_sensor <- gsub(paste0(".*", time_sep), "", sensor_times)
    if (switch <- 1) {
      time_of_sensor <- paste0(time_of_sensor, ":00")
    }
    sensor_chron <- c(data$datetime[1],
                      chron::chron(dates=date_of_sensor, times=time_of_sensor, format=time_format),
                      data$datetime[length(data$datetime)])
    data_list <- list()
    for (i in 1:(length(sensor_times)+1)) {
      cutpoint1 <- min(which(data$datetime >= sensor_chron[i]))
      cutpoint2 <- min(which(data$datetime >= sensor_chron[i+1])) - 1
      pertinent <- data[cutpoint1:cutpoint2, c(2, 6)]
      if (na.rm == TRUE) pertinent <- pertinent[which(!is.na(pertinent$glucose) & !is.na(pertinent$elapsed)), ]
      data_list[[i]] <- pertinent
    }
    return(data_list)
  }

  data <- data[, c(2, 6)] # keep only glucose value and elapsed time

  if (na.rm == TRUE) data <- data[which(!is.na(data$glucose) & !is.na(data$elapsed)), ]

  # Multiple sensors
  if (mult_sensors == FALSE) return(data)
  if (mult_sensors == TRUE & is.na(sensor_times)) {
    # look for two hour gaps
    gaps <- data$elapsed[2:length(data$elapsed)] - data$elapsed[1:(length(data$elapsed)-1)]
    cutpoints <- c(0, which(gaps > sensor_gap), length(data$elapsed))
    data.list <- list()
    for (i in 1:(length(cutpoints)-1)) {
      data.list[[i]] <- data[(cutpoints[i]+1):cutpoints[i+1], ]
    }
    return(data.list)
  }
}
