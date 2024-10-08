% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simu.ph23data.R
\name{simu.ph23data}
\alias{simu.ph23data}
\title{Simulate randomized and controlled phase 2/3 dose optimization trial for survival endpoint}
\usage{
simu.ph23data(
  nSim = 1000,
  n1 = c(50, 50, 50, 50),
  n2 = c(200, 200),
  m = c(9, 11, 13, 8),
  Lambda1 = function(t) {
     (t/12) * as.numeric(t <= 12) + as.numeric(t > 12)
 },
  A1 = 12,
  DCO1 = 16,
  Lambda2 = function(t) {
     (t/12) * as.numeric(t <= 12) + as.numeric(t > 12)
 },
  A2 = 12,
  enrollment.hold = 4,
  targetEvents2 = c(300, 380),
  method = "Independent Incremental"
)
}
\arguments{
\item{nSim}{Number of simulated trials}

\item{n1}{Sample size of (dose arm, control arm) at Stage 1. length(n1) must be 2.}

\item{n2}{Sample size of the (selected dose arm, control arm). length(n2) must be 2.}

\item{m}{Median survival time for each arm (dose 1, dose 2, ..., control). length(m) must be equal to length(n1)}

\item{Lambda1}{Enrollment distribution function (CDF) for stage 1.}

\item{A1}{Enrollment period for Stage 1}

\item{DCO1}{Data cutoff date for Stage 1}

\item{Lambda2}{Enrollment distribution function (CDF) for stage 2.}

\item{A2}{Enrollment period for Stage 2}

\item{enrollment.hold}{Holding period in months after DCO1 of Stage 1 prior to enrollment of Stage 2 patients. 0 means seamless enrollment.}

\item{targetEvents2}{Planned target number of events for Stage 2. Either targetEvents2 must be provided.}

\item{method}{"Independent Incremental", "Disjoint Subjects". Currently, only "Independent Incremental" method is implemented.
The multiplicity adjustment is performed according to the extended followup data up to 1st analysis. This option is applicable to "Disjoined Subjects" method. Currently, not implemented yet.}
}
\value{
\describe{
\item{z1}{Stage 1 z values for each dose}
\item{z.c}{z value for each analysis at Stage 2 for the selected dose}
\item{z2}{Incremental z statistic at each anaysis at Stage 2, calculated from dose selection to the analysis}
\item{w}{Weight, calculated as sqrt (total events for selected dose + control at stage 1 / targetEvents at each analysis at Stage 2 for the selected dose + control combining stage 1 and stage 2 patients.)}
\item{selected.dose}{Selected dose}
\item{example.data}{A simulated dataset example}
}
}
\description{
This functions simulates trials that have multiple dose arms at Stage 1 and the best dose is always selected after Stage 1 and perform multiple analyses at Stage 2
The function returns multiple datasets including at end of Stage 1 (all dose arms) and at each of analysis at Stage 2.
Features include: (1) Stage 1 and Stage 2 have separate enrollment curves; (2) Allow enrollment gap between Stage 1 and Stage 2;
(3) Dose selection at stage 1 based on best performance in z statistic. (4) Stage 2 analyses are based on prespecified target events. (5) The function
also returns the weights for combination p value approach in the next step and the weights are sqrt of information fractions.
}
\examples{
#Example (1): Stage 1: 4 arms; 3 dose levels; each arm 50 patients.
#Stage 2: additional 200 patients per arm will be enrolled at stage 2
#medians for the 4 arms: 9, 11, 13 and control = 8 months
#Enrollment: 12 months uniform in stage 1; 12 months uniform in stage 2
#Holding period: 4 months between stage 1 and 2
#Dose selection will be based on data cut at 16 months
#Stage 2 has 2 planned analyses at 300 and 380 events respectively.

#Using O'Brien Fleming boundary, the rejection boundaries are: 
bd.z = actualBounds(planned.events=c(300, 380), act.events=c(300, 380), sf=gsDesign::sfLDOF, alpha=0.025)$actual.z
#2.268527 2.022098

o = simu.ph23data(nSim=100, n1 = c(50, 50, 50, 50), n2 = c(200, 200), m = c(9,9, 9, 9), 
Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A1 = 12,
DCO1 = 16, Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A2 = 12,
enrollment.hold=4, targetEvents2 = c(300, 380))

#Perform analysis at IA and FA using combination p method, with simes approach for multiplcity control at Stage 1 and closed testing procedure for controlling FWER
IA = ph23.comb.p(z1=o$z1,  z2 = o$z2[,1], bd.z=bd.z[1], w=o$w[,1])
FA = ph23.comb.p(z1=o$z1,  z2 = o$z2[,2], bd.z=bd.z[2], w=o$w[,2])

#power calculation using the standard group sequential boundaries
gsd.power(z = cbind(IA$comb.z, FA$comb.z), bd.z=bd.z)

}
