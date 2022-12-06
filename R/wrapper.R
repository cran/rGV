#' Calculate all glycemic variability metrics
#'
#' @param x vector of glucose readings
#' @param times vector of corresponding times, in minutes
#' @param unit "mg" if the units are mg/dL or "mmol" if the units are mmol/L. Null value is "mg".
#' @param m_index a value to be considered a 'standard' blood glucose value for calculation of M-value, in mg/dL. Null value is 120.
#' @param k length of time (in minutes) used to find partners. Null value is 60.
#' @param s number of minutes of slack used when searching for partners. Null value is 1.
#' @param conga_n number of hours between "partner" observations. Null value is 1.
#' @param interval size (in hours) of the moving window to be used if overall is false. Null value is 1.
#' @param thresh threshold above (or below) which you wish to calculate the AUC. Default is 100.
#' @param event_thresh glucoses below this threshold are considered as part of an episode. Default is 55
#'
#' @return A data frame containing the entire suite of rGV metrics for the given dataset.
#' @export
#'
#' @examples
#' GV(x=c(rep(100, 10), rep(120, 10), 105, 85), times=seq(0, 1260, 60), unit='mg',
#' m_index=120, k=60, s=1, conga_n=1, interval=1, thresh=100, event_thresh=55)
GV <- function(x, times, unit="mg", m_index=120, k=60, s=1, conga_n=1, interval=1, thresh=100,
               event_thresh=55) {
  m1 <- mean(x)
  m2 <- st_dev(x, times, overall=TRUE, interval=interval)
  m3_1 <- conga(x, times, n=conga_n, s=s, method="manuscript")
  m3_2 <- conga(x, times, n=conga_n, s=s, method="easy")
  m4 <- li(x, times, k=k, s=s)
  m5 <- j_index(x, unit=unit)
  m6_1 <- bgi(x,  unit=unit, method="manuscript")
  m6_2 <- bgi(x,  unit=unit, method="easy")
  m7_1 <- grade(x,  unit=unit, method="manuscript")
  m7_2 <- grade(x,  unit=unit, method="easy")
  m8_1 <- modd(x, times, method="manuscript")
  m8_2 <- modd(x, times, method="easy")
  m9 <- mage(x, times)
  m10_1 <- adrr(x, times, unit=unit, method="manuscript")
  m10_2 <- adrr(x, times, unit=unit, method="easy")
  m11_1 <- m_value(x, unit=unit, index=m_index, method="manuscript")
  m11_2 <- m_value(x, unit=unit, index=m_index, method="easy")
  m12 <- mag(x, times)
  m13 <- cv(x, times, overall=TRUE, interval=interval)
  m14 <- cgm_auc(x, times, thresh=thresh, above=TRUE)
  m15 <- tir(x, low=251, high=Inf)
  m16 <- tir(x, low=181, high=250)
  m17 <- tir(x, low=70, high=180)
  m18 <- tir(x, low=54, high=69)
  m19 <- tir(x, low=-Inf, high=53)
  m20 <- gmi(x, unit)
  m21 <- num_events(x, times, thresh=event_thresh, len=15)
  m22 <- gvp(x, times)
  m23 <- dist_travelled(x)

  vec1 <- c("Mean", "Standard Deviation", "CONGA", "Lability Index", "J-index",
            "LBGI", "HBGI", "GRADE", "GRADE Hypo %",
            "GRADE Eu %", "GRADE Hyper %",
            "MODD", "MAGE",
            "ADRR High", "ADRR Low", "M-value", "MAG", "CV", "AUC", "% readings above 250",
            "% readings between 181-250", "% readings between 70-180", "% readings between 54-69",
            "% readings below 54", "GMI", paste0("Number of Episodes below ", event_thresh),
            "GVP", "Distance Travelled")
  vec2 <- unlist(c(m1, m2, m3_1, m4, m5, m6_1[1], m6_1[2], m7_1[1], m7_1[2],
                   m7_1[3], m7_1[4], m8_1, m9, m10_1, NA, m11_1, m12, m13, m14, m15, m16, m17, m18, m19, m20, m21, m22, m23), use.names=FALSE)
  vec3 <- unlist(c(m1, m2, m3_2, m4, m5, m6_2[1], m6_2[2], m7_2[1], m7_2[2],
                   m7_2[3], m7_2[4], m8_2, m9, m10_2[1], m10_2[2], m11_2, m12, m13, m14, m15, m16, m17, m18, m19, m20, m21, m22, m23), use.names=FALSE)
  results <- data.frame(vec1, vec2, vec3)
  names(results) <- c("Measure", "Manuscript", "EasyGV")
  results$Manuscript <- round(results$Manuscript, 2)
  results$EasyGV <- round(results$EasyGV, 2)
  results
}
