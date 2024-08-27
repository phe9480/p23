#' Expected Number of Events for Stage 1 Subjects at Time of An Analysis With Planned Target Number of Events
#'
#' This functions calculates the expected number of events for stage 1 subjects at time of analysis with planned target number of events.
#' This function is useful for determining a pre-fixed weight using Disjoint Subjects approach in Sample size re-estimation.
#' 
#' @param n1 Stage 1 sample size of each dose and control arm. length(n1) must be number of arms.
#' @param n2 Stage 2 Sample size of the selected dose and control arm. length(n2) must be 2.
#' @param m Median survival time for each arm (dose 1, dose 2, ..., control). length(m) must be equal to length(n1)
#' @param A1 Enrollment period for Stage 1
#' @param Lambda1 Enrollment distribution function (CDF) for stage 1.
#' @param enrollment.hold Holding period in months after DCO1 of Stage 1 prior to enrollment of Stage 2 patients. 0 means seamless enrollment.
#' @param Lambda2 Enrollment distribution function (CDF) for stage 2.
#' @param targetEvents Planned target number of events for Stage 2.
#' 
#' @return Expected number of events for each dose arm + control for stage 1 subjects at time of analysis when the target events achieved.
#' 
#' 
#' @examples
#' 
#' 
#' e1.ssr(n1 = rep(50, 4), n2 = rep(200, 4), m = c(9,9, 9, 9), 
#'      Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, 
#'      A1 = 12, Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)},
#'      enrollment.hold=4, targetEvents = 380)
#' 
#' @export 
#' 
e1.ssr = function(n1 = rep(50, 4), n2 = rep(200, 4), m = c(9,9, 9, 9), 
                  Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, 
                  A1 = 12,
                  Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)},
                  enrollment.hold=4, targetEvents = 380){
  n.arms = length(n1)
  E1 = rep(NA, n.arms-1)
  for (i in 1:(n.arms-1)) {
    #Find DCO to achieve original targetEvents events
    f.u = function(u){
      e1 = p23::fe(DCO = u, r = (n1[i] + n2[i])/(n1[n.arms] + n2[n.arms]), 
              h0 = function(t){log(2)/m[n.arms]}, S0 = function(t){exp(-log(2)/m[n.arms]*t)}, 
              h1 = function(t){log(2)/m[i]}, S1 = function(t){exp(-log(2)/m[i]*t)}, 
              Lambda = Lambda1, n = sum(n1[c(i, n.arms)]))$e
      e2 = p23::fe(DCO = max(u-enrollment.hold-A1, 0), r = (n1[i] + n2[i])/(n1[n.arms] + n2[n.arms]), 
              h0 = function(t){log(2)/m[n.arms]}, S0 = function(t){exp(-log(2)/m[n.arms]*t)}, 
              h1 = function(t){log(2)/m[i]}, S1 = function(t){exp(-log(2)/m[i]*t)}, 
              Lambda = Lambda2, n = sum(n2[c(i, n.arms)]))$e
      return(e1 + e2 - targetEvents)
    }
    DCO.FA0 = uniroot(f.u, c(16, 1000))$root
    
    #Expected events for Stage 1 at original FA: e1.FA0
    E1[i]=p23::fe(DCO = DCO.FA0, r = (n1[i] + n2[i])/(n1[n.arms] + n2[n.arms]), 
              h0 = function(t){log(2)/m[n.arms]}, S0 = function(t){exp(-log(2)/m[n.arms]*t)}, 
              h1 = function(t){log(2)/m[i]}, S1 = function(t){exp(-log(2)/m[i]*t)}, 
              Lambda = Lambda1, n = sum(n1[c(i, n.arms)]))$e
  }
  return(E1)
}


