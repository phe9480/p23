% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/corr.DisjointSubjects.p23.ssr.R
\name{corr.DisjointSubjects.p23.ssr}
\alias{corr.DisjointSubjects.p23.ssr}
\title{Correlation between the test statistics at IA and FA Using Disjoint Subjects Method}
\usage{
corr.DisjointSubjects.p23.ssr(t11, t1, t21, t21star, t2star)
}
\arguments{
\item{t11}{Information fraction for Stage 1 subjects at IA.}

\item{t1}{Information fraction for both stage 1 and 2 subjects at IA.}

\item{t21}{Information fraction for Stage 1 subjects at FA prior to SSR.}

\item{t21star}{Information fraction for Stage 1 subjects at FA post SSR.}

\item{t2start}{Information fraction for both stage 1 and 2 subjects at FA post SSR. If there is increase of target events at FA after SSR, then t2star > 1.}
}
\value{
correlation between test statistics at IA and FA
}
\description{
This functions calculates the correlation between the test statistics (Z.tilde) at IA and FA Using disjoint subjects method.
The calculation is based on information fractions, relative to the information for both stage 1 and 2 subjects at FA (i.e. information = 1).
}
\examples{

#At IA: Stage 1 subjects 50 events (e11)
#       Stage 1 and 2 subjects 110 (e1)
#At FA (originally planned):
#       Stage 1 subjects are estimated having 70 events (e21)
#       Stage 1 and 2 subjects 180 events (e2)
#At FA (after SSR):
#       Stage 1 subjects 80 events (e21star)
#       Stage 1 and 2 subjects 196 (e2star)
#The information fractions: t11 = e11/e2 = 50/180; t1 = e1/e2 = 110/180
#t21 = e21/e2 = 70/180; t2 = 1; t21star = e21star/e2 = 80/180; t2star = e2star/e2 = 196/180

#The correlation between IA test statistic and FA test statistic
corr.DisjointSubjects.p23.ssr(t11=50/180, t1=110/180, t21=70/180, t21star=80/180, t2star=196/180)

#The correlation between IA and FA test statistics in the original design without SSR
sqrt(110/180) #sqrt(t1)

}
