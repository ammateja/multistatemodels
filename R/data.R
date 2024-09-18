#' Illness Death Model
#'
#' Simulated data from an illness death model
#' N=250 subjects
#' Subjects can be in state 1 (healthy), state 2 (illness), or state 3 (death, absorbing)
#' Subjects can transition from state 1 to 2, 1 to 3, or 2 to 3
#' Transition intensities are based on a Weibull distribution
#'
#'
#' @format ## `illness_death_dat`
#' A data frame with 1992 rows and 6 columns:
#' \describe{
#'   \item{id}{Subject ID, numeric, 1, 2, 3, etc.}
#'   \item{tstart}{Interval start time}
#'   \item{tstop}{Interval stop time}
#'   \item{statefrom}{State at interval start time}
#'   \item{stateto}{State at interval stop time}
#'   \item{obstype}{1=exactly observed data, 2=panel data}
#' }
#' @source /data-raw/DATASET.R
"illness_death_dat"
