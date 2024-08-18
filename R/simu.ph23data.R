#' Simulate randomized and controlled phase 2/3 dose optimization trial for survival endpoint[DO NOT USE]
#'
#' This functions simulates trials that have multiple dose arms at Stage 1 and the best dose is always selected after Stage 1 and perform multiple analyses at Stage 2
#' The function returns multiple datasets including at end of Stage 1 (all dose arms) and at each of analysis at Stage 2
#'
#' @param n1 Sample size of (dose arm, control arm) at Stage 1. length(n1) must be 2.
#' @param n2 Sample size of the (selected dose arm, control arm). length(n2) must be 2.
#' @param m Median survival time for each arm (dose 1, dose 2, ..., control). length(m) must be equal to length(n1)
#' @param Lambda2 Enrollment distribution function (CDF) for stage 2.
#' @param Lambda1 Enrollment distribution function (CDF) for stage 1.
#' @param DCO1 Data cutoff date for Stage 1
#' @param targetEvents Planned target number of events for Stage 2. Either targetEvents or DCO must be provided. 
#' 
#' @return Datasets including data1: Stage 1 data with multiple dose arms; data2: Stage 2 data of the selected dose arm and control arm with data2[[j]] for jth analysis. Each dataset includes variables
#' \describe{
#' \item{treatment}{treatment group with values of "control", "dose 1", "dose 2", ...}
#' \item{enterTime}{Time of randomization in calendar time}
#' \item{calendarTime}{the time when event/censoring occurred in calendar time}
#' \item{survTime}{Survival time for analysis, = calendarTime - enterTime}
#' \item{cnsr}{censor status (0=event; 1=censor) before administrative censoring due to data cut}
#' \item{calendarCutOff}{Data CutOff Time (DCO);}
#' \item{survTimeCut}{Survival time after cut}
#' \item{cnsrCut}{Censor status after cut}
#' }
#' The function also returns a dataframe decision with variables: 
#' \describe{
#' \item{z}{log rank z value for each dose compared to control using stage 1 data}
#' \item{p}{log rank p value for each dose compared to control using stage 1 data}
#' \item{adj.p}{adjusted log rank p value for selected dose arm}
#' \item{adj.z}{adjusted log rank p value for selected dose arm}
#' \item{selected}{Selected dose arm, applicable for data2 only}
#' }
#' 
#' @examples
#' #Example (1): 
#' 
#' p = c(0.01, 0.02, 0.03, 0.013)
#' simes(p)
#' 
#' @export 
#' 
 simu.ph23data = function(n1 = c(50, 50), n2 = c(200, 200), m = c(9, 11, 13, 8), 
                          Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A1 = 12,
                          Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A2 = 12,
                          DCO1 = 16, targetEvents = c(300, 380)){
 
   #number of analyses at Stage 2
   L = length(targetEvents)
   
   #(1) Generate exponential survival data
   n.doses = length(m) - 1
   m0 = m[length(m)] #control arm median
   
   #All dose arms have the same n
   t = matrix(NA, nrow=n1[1]+n2[1], ncol=n.doses) #dose arms
   t0 = genSurv(dist = "exponential", n = n1[2]+n2[2], lam=log(2)/m0)
   
   for (i in 1:n.doses){
     t[,i] = genSurv(dist = "exponential", n = n1[1]+n2[1], lam=log(2)/m[i])
   }
   
   #(2) Enrollment data
   r = (n1[1]+n2[1]) / (n1[2]+n2[2])
   
   nEachMonth1 = f.nEachMonth(N=sum(n1), A=12, w=NULL, r=r, Lambda=Lambda1)
   
   out = list(NULL)
   for (k in 2:(L+1)){out = c(out, list(NULL))}
   
   ############################
   #EnterTime, CalendarTime
   ############################
   enterTime = function(nPerMonth, A){
     enterT = rep(NA, n1[2])
     enterT[1:nEachMonth$n0[1]] = runif(nEachMonth$n0[1], min=0, max=1)
     if (A1 > 1) {for (m in 2:A1){
       LL = sum(nEachMonth$n0[1:(m-1)])+1
       UU = sum(nEachMonth$n0[1:m])
       enterTime0[LL:UU] = runif(nEachMonth$n0[m], min=m-1, max=m)
     }}
   }
   enterTime0 = rep(NA, n1[2])
   enterTime0[1:nEachMonth$n0[1]] = runif(nEachMonth$n0[1], min=0, max=1)
   if (A1 > 1) {for (m in 2:A1){
     LL = sum(nEachMonth$n0[1:(m-1)])+1
     UU = sum(nEachMonth$n0[1:m])
     enterTime0[LL:UU] = runif(nEachMonth$n0[m], min=m-1, max=m)
   }}
   
   enterTime1 = rep(NA, n1[1])
   enterTime1[1:nEachMonth$n1[1]] = runif(nEachMonth$n1[1], min=0, max=1)
   if (A1 > 1) {for (m in 2:A1){
     LL = sum(nEachMonth$n1[1:(m-1)])+1
     UU = sum(nEachMonth$n1[1:m])
     enterTime1[LL:UU] = runif(nEachMonth$n1[m], min=m-1, max=m)
   }}
   
   treatment = c(rep("control", n1[2]), rep("experimental", n1[1]))
   enterTime = c(enterTime0, enterTime1)
   survTime = as.numeric(c(t0, t[]))  
    
 }

