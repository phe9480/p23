% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simu.single.arm.R
\name{simu.single.arm}
\alias{simu.single.arm}
\title{Simulate a single arm survival data with non-uniform enrollment distribution
and drop-off}
\usage{
simu.single.arm(
  n = 100,
  m = 10,
  Lambda = function(t) {
     (t/12) * as.numeric(t <= 12) + as.numeric(t > 12)
 },
  A = 12,
  drop = 0,
  DCO = 16,
  targetEvents = NULL
)
}
\arguments{
\item{n}{Sample size}

\item{m}{Median survival time}

\item{Lambda}{Enrollment distribution function (CDF)}

\item{A}{Enrollment period}

\item{drop}{Dropout rate per month}

\item{DCO}{Data cutoff date for Stage 1}

\item{targetEvents}{Planned target number of events. Either targetEvents or
DCO must be provided.}
}
\value{
Datasets including data1: Stage 1 data with multiple dose arms;
data2: Stage 2 data of the selected dose arm and control arm with data2((j))
for jth analysis. Each dataset includes variables
\describe{
\item{treatment}{treatment group with values of "control", "dose 1",
"dose 2", ...}
\item{enterTime}{Time of randomization in calendar time}
\item{calendarTime}{the time when event/censoring occurred in calendar time}
\item{survTime}{Survival time for analysis, = calendarTime - enterTime}
\item{cnsr}{censor status (0=event; 1=censor) before administrative
censoring due to data cut}
\item{calendarCutOff}{Data CutOff Time (DCO);}
\item{survTimeCut}{Survival time after cut}
\item{cnsrCut}{Censor status after cut}
\item{analysis}{Sequence of analysis}
}
}
\description{
This functions simulates trials that have multiple dose arms at Stage 1 and
the best dose is always selected after Stage 1 and perform multiple analyses
at Stage 2. The function returns multiple datasets including at end of Stage
1 (all dose arms) and at each of analysis at Stage 2
}
\examples{
 
Lambda = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}  
dat=simu.single.arm(n = 100, m = 10, Lambda = Lambda, A = 12, 
drop=0, DCO = 16, targetEvents = NULL)

dat=simu.single.arm(n = 100, m = 10, Lambda = Lambda, A = 12, 
drop=0, DCO = c(16, 28), targetEvents = NULL)

dat=simu.single.arm(n = 100, m = 10, Lambda = Lambda, A = 12, drop=0, 
DCO = NULL, targetEvents = 40)

dat=simu.single.arm(n = 100, m = 10, Lambda = Lambda, A = 12, drop=0, 
DCO = NULL, targetEvents = c(40, 60))

}
