% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/e1.ssr.R
\name{e1.ssr}
\alias{e1.ssr}
\title{Expected Number of Events for Stage 1 Subjects at Time of An Analysis With Planned Target Number of Events}
\usage{
e1.ssr(
  n1 = rep(50, 4),
  n2 = rep(200, 4),
  m = c(9, 9, 9, 9),
  Lambda1 = function(t) {
     (t/12) * as.numeric(t <= 12) + as.numeric(t > 12)
 },
  A1 = 12,
  Lambda2 = function(t) {
     (t/12) * as.numeric(t <= 12) + as.numeric(t > 12)
 },
  enrollment.hold = 4,
  targetEvents = 380
)
}
\arguments{
\item{n1}{Stage 1 sample size of each dose and control arm. length(n1) must be number of arms.}

\item{n2}{Stage 2 Sample size of the selected dose and control arm. length(n2) must be 2.}

\item{m}{Median survival time for each arm (dose 1, dose 2, ..., control). length(m) must be equal to length(n1)}

\item{Lambda1}{Enrollment distribution function (CDF) for stage 1.}

\item{A1}{Enrollment period for Stage 1}

\item{Lambda2}{Enrollment distribution function (CDF) for stage 2.}

\item{enrollment.hold}{Holding period in months after DCO1 of Stage 1 prior to enrollment of Stage 2 patients. 0 means seamless enrollment.}

\item{targetEvents}{Planned target number of events for Stage 2.}
}
\value{
Expected number of events for each dose arm + control for stage 1 subjects at time of analysis when the target events achieved.
}
\description{
This functions calculates the expected number of events for stage 1 subjects at time of analysis with planned target number of events.
This function is useful for determining a pre-fixed weight using Disjoint Subjects approach in Sample size re-estimation.
}
\examples{


e1.ssr(n1 = rep(50, 4), n2 = rep(200, 4), m = c(9,9, 9, 9), 
     Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, 
     A1 = 12, Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)},
     enrollment.hold=4, targetEvents = 380)

}
