#' Power Calculation by Simulations for Phase 2/3 Inferential Seamless Design with Sample Size Adjustment
#'
#' This functions calculates the cumulative power by Simulations for Phase 2/3 Inferential Seamless Design with Sample Size Adjustment.
#' The SSR is only performed at the first interim analysis at Stage 2. Regarding the pre-fixed weights for the FA combination Z,
#' (1) If the method is independent incremental, then weight maintains is based on the information fraction as in the original target events.
#' (2) If the method is disjoint subjects, then the weight is calculated based on the expected number of events for Stage 1 subjects 
#' at final analysis relative to the original target events at FA. (3) If the method is mixture, then
#' the weight for FA is same as the independent increment method.
#'
#' @param nSim Number of simulated trials
#' @param n1 Stage 1 sample size of each dose and control arm. length(n1) must be number of arms.
#' @param n2 Stage 2 Sample size of the selected dose and control arm. length(n2) must be 2.
#' @param m Median survival time for each arm (dose 1, dose 2, ..., control). length(m) must be equal to length(n1)
#' @param orr ORR for each arm. length(orr) = length(m). 
#' @param rho Correlation between ORR and time to event endpoint
#' @param dose_selection_endpoint  Dose selection end point: "ORR" or "not ORR"
#' @param ssr_HR_threshold The threshold in HR for sample size (target events) increase at first analysis in Stage 2.
#' @param events_increase Number of events increased to the original target number of events at final analysis when the HR crosses the threshold at IA
#' @param A1 Enrollment period for Stage 1
#' @param Lambda1 Enrollment distribution function (CDF) for stage 1.
#' @param DCO1 Data cutoff date for Stage 1
#' @param enrollment.hold Holding period in months after DCO1 of Stage 1 prior to enrollment of Stage 2 patients. 0 means seamless enrollment.
#' @param A2 Enrollment period for Stage 2
#' @param Lambda2 Enrollment distribution function (CDF) for stage 2.
#' @param targetEvents2 Originally planned target number of events for Stage 2. Only allow 1 IA. length(targetEvents2) must be 2.
#' @param alpha Type I error (one-sided) for testing the selected dose, usually 0.025.
#' @param sf Spending functions. acceptable options include all spending functions in gsDesign R package, for example, "gsDesign::sfLDOF"
#' @param multiplicity.method "simes" or "Dunnett"
#' @param method Options include "Independent Incremental": z1 at dose selection and z2 is from dose selection to kth analysis at stage 2; 
#' "Disjoint Subjects": z1 is at kth analysis for stage 1 subjects; z2 is at the kth analysis for stage 2 subjects. z1 will be adjusted by multiplicity and closed testing procedure at each analysis.
#' "Mixture": Only consider disjoint subjects at first analysis in stage 2. Starting from the 2nd analysis, consider independent incremental methods. Only z1 at 1st analysis will be adjusted by multiplicity and closed testing procedure.
#' 
#' @return An object with values:
#' \describe{
#' \item{bd.z}{z value rejection boundary at each analysis if method is independent incremental or mixture.}
#' \item{cum.pow}{Cumulative power}
#' \item{s}{Selected dose}
#' \item{selection}{Probability of selection among doses}
#' }
#' 
#' 
#' @examples
#' #Example (1): Stage 1: 4 arms; 3 dose levels; each arm 50 patients.
#' #Stage 2: additional 200 patients per arm will be enrolled at stage 2
#' #medians for the 4 arms: 9, 11, 13 and control = 8 months
#' #Enrollment: 12 months uniform in stage 1; 12 months uniform in stage 2
#' #Holding period: 4 months between stage 1 and 2
#' #Dose selection will be based on data cut at 16 months
#' #Stage 2 has 2 planned analyses at 300 and 380 events respectively.
#' #At 1st analysis in Stage 2, the target number of events will increase 30 
#' when the observed HR > 0.8.
#'
#' #Dose selection decision is NOT based on ORR.
#' simu.power.p23.ssr(nSim=1000, n1 = rep(50, 4), n2 = rep(200, 2), m = c(9, 9, 9, 9), 
#' orr = NULL, rho = NULL, dose_selection_endpoint = "not ORR",
#' ssr_HR_threshold = 0.8, events_increase = 30,
#' Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A1 = 12,
#' Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A2 = 12,
#' enrollment.hold=4, DCO1 = 16, targetEvents2=c(300, 380), sf=gsDesign::sfLDOF, 
#' alpha=0.025, method = "Independent Incremental")
#' 
#' #Example (2): #Dose selection decision based on ORR
#' simu.power.p23.ssr(nSim=10, n1 = rep(50, 4), n2 = rep(200, 4), m = c(9, 9, 9, 9), 
#' orr = c(0.25, 0.3, 0.4, 0.2), rho = 0.7, dose_selection_endpoint = "ORR",
#' ssr_HR_threshold = 0.8, events_increase = 30,
#' Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A1 = 12,
#' Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A2 = 12,
#' enrollment.hold=4, DCO1 = 16, targetEvents2=c(300, 380), sf=gsDesign::sfLDOF, 
#' alpha=0.025, method = "Independent Incremental")
#' 
#' simu.power.p23.ssr(nSim=10, n1 = rep(50, 4), n2 = rep(200, 4), m = c(9, 9, 9, 9), 
#' orr = c(0.25, 0.3, 0.3, 0.2), rho = 0.7, dose_selection_endpoint = "ORR",
#' ssr_HR_threshold = 0.8, events_increase = 30,
#' Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A1 = 12,
#' Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A2 = 12,
#' enrollment.hold=4, DCO1 = 16, targetEvents2=c(300, 380), sf=gsDesign::sfLDOF, 
#' alpha=0.025, method = "Disjoint Subjects")
#' 
#' simu.power.p23.ssr(nSim=10, n1 = rep(50, 4), n2 = rep(200, 4), m = c(9, 9, 9, 9), 
#' orr = c(0.25, 0.3, 0.3, 0.2), rho = 0.7, dose_selection_endpoint = "ORR",
#' ssr_HR_threshold = 0.8, events_increase = 30,
#' Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A1 = 12,
#' Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A2 = 12,
#' enrollment.hold=4, DCO1 = 16, targetEvents2=c(300, 380), sf=gsDesign::sfLDOF, 
#' alpha=0.025, method = "Mixture")
#' 
#' o=simu.power.p23.ssr(nSim=100, n1 = rep(50, 4), n2 = rep(200, 4), m = c(9, 9, 9, 9), 
#' orr = c(0.25, 0.3, 0.3, 0.2), rho = 0.7, dose_selection_endpoint = "ORR",
#' ssr_HR_threshold = 0.8, events_increase = 30,
#' Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, 
#' A1 = 12,Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, 
#' A2 = 12,enrollment.hold=4, DCO1 = 16, targetEvents2=c(300, 380), 
#' sf=gsDesign::sfLDOF, alpha=0.025, multiplicity.method = "dunnett", 
#' method = "Disjoint Subjects")
#' 
#' o=simu.power.p23.ssr(nSim=100, n1 = rep(50, 4), n2 = rep(200, 4), m = c(9, 12, 9, 9), 
#' orr = c(0.25, 0.3, 0.2, 0.2), rho = 0.7, dose_selection_endpoint = "ORR",
#' ssr_HR_threshold = 0.8, events_increase = 30,
#' Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, 
#' A1 = 12,Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, 
#' A2 = 12,enrollment.hold=4, DCO1 = 16, targetEvents2=c(300, 380), 
#' sf=gsDesign::sfLDOF, alpha=0.025, multiplicity.method = "dunnett", 
#' method = "Disjoint Subjects")
#' 
#' 
#' @importFrom gsDesign gsDesign
#' @importFrom stats uniroot
#' 
#' @export 
#' 
simu.power.p23.ssr = function(nSim=100, n1 = rep(50, 4), n2 = rep(200, 2), m = c(9,9, 9, 9), 
                          orr = c(0.25, 0.3, 0.4, 0.2), rho = 0.7, 
                          dose_selection_endpoint = "ORR",
                          ssr_HR_threshold = 0.8, events_increase = 30, 
                          Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A1 = 12,
                          Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A2 = 12,
                          enrollment.hold=4, DCO1 = 16, targetEvents2=c(300, 380), 
                          alpha=0.025, sf=gsDesign::sfLDOF, multiplicity.method="simes",
                          method = "Independent Incremental"){
  
  #Number of analyses in stage 2
  K = length(targetEvents2)
  if (K != 2) {stop; return("ERROR: Only 2 analyses in stage 2 are considered in SSR. The length of targetEvents2 must be 2.")}

  #Number of arms
  n.arms = length(n1)
  n2 = c(rep(n2[1], n.arms-1), n2[2])
  
  #Expected number of events for Stage 1 subjects at time of originally planned IA and FA
  e1.FA0 = e1.ssr(n1 = n1, n2 = n2, m = m, Lambda1 = Lambda1, A1 = A1,
                  Lambda2 = Lambda2, enrollment.hold=enrollment.hold, 
                  targetEvents = targetEvents2[2])
  
  #rejection boundary by GSD. The original GSD boundary is still valid in this SSR setting.
  #for independent incremental and mixture approaches
  if (method == "Independent Incremental" || method == "Mixture"){
    bd.z = gsDesign::gsDesign(k=K,alpha=alpha,timing=targetEvents2/targetEvents2[K],sfu=sf, test.type=1)$upper$bound
  }
  
  #Combination Z values
  comb.z = matrix(NA, nrow=nSim, ncol=K)
  s = rep(NA, nSim) #selected dose
  new.ss = rep(NA, nSim) #indicator whether sample size increased
  DS.cum.rej1 = DS.cum.rej2 = rep(NA, nSim) #cum rejection at IA or FA for DS method
  DS.bd.z = matrix(NA, nrow=nSim, ncol=K)
  
  for (i in 1:nSim){
    p23i = simu.p23trial(n1 = n1, n2 = n2, m = m, 
                         orr = orr, rho = rho, dose_selection_endpoint = dose_selection_endpoint,
                         Lambda1 = Lambda1, A1 = A1, 
                         Lambda2 = Lambda2, A2 = A2, enrollment.hold=enrollment.hold)
    #dose selection
    sel = select.dose.p23 (data=p23i, DCO1=DCO1, dose_selection_endpoint = dose_selection_endpoint)
    s[i] = sel$s
    
    #SSR
    ssr = ssr.p23(data=p23i, ssr_HR_threshold = ssr_HR_threshold, 
                  events_increase = events_increase, 
                  selected.dose=s[i], targetEvents2 = targetEvents2)
    
    #conduct p23
    o=conduct.p23.ssr(data=p23i, DCO1=DCO1, dose_selection_endpoint = dose_selection_endpoint, 
                      targetEvents.FA = ssr$targetEvents.FA, e1.FA = e1.FA0[s[i]],
                      targetEvents2 = targetEvents2, method = method, multiplicity.method=multiplicity.method)
    
    #Test statistics z.tilde
    comb.z[i, ] = o$z.tilde
    new.ss[i] = ssr$ssr
    
    #When method = DS, the traditional GSD boundary needs adjustment
    if (method == "Disjoint Subjects"){
      dat23 = p23i[p23i$group == 0 | p23i$group == s[i], ]
      dat23.1 = f.dataCut(data=dat23, targetEvents=targetEvents2[1])
      dat23.11 = dat23.1[dat23.1$stage == 1, ]
      
      #FA after SSR
      dat23.2star = f.dataCut(data=dat23, targetEvents=ssr$targetEvents.FA)
      dat23.21star = dat23.2star[dat23.2star$stage == 1, ] #Analysis 2 (FA) stage 1
      
      #Observed events for Stage 1 subjects at IA and FA
      N11 = sum(1-dat23.11$cnsrCut)
      N21star = sum(1-dat23.21star$cnsrCut)
      
      t11 = N11 / targetEvents2[2]; 
      t1 = targetEvents2[1] / targetEvents2[2]
      t21 = e1.FA0[s[i]] / targetEvents2[2]; #pre-fixed weight for FA at time of IA
      t21star =  N21star/ targetEvents2[2];
      t2star = ssr$targetEvents.FA / targetEvents2[2]
      
      correl = corr.DisjointSubjects.p23.ssr(t11=t11, t1=t1, t21=t21, t21star=t21star, t2star=t2star)
      
      bd.zi = gsDesign::gsDesign(k=K,alpha=alpha,timing=correl^2,sfu=sf, test.type=1)$upper$bound
      
      #IA
      DS.cum.rej1[i] = as.numeric(o$z.tilde[1] >=  bd.zi[1])
      #FA
      DS.cum.rej2[i] = as.numeric(o$z.tilde[1] >=  bd.zi[1] || o$z.tilde[2] >=  bd.zi[2])
      
      DS.bd.z[i, ] = bd.zi
    }
  }
  
  cum.pow = rep(NA, K)
  if (method == "Disjoint Subjects") {
    cum.pow[1] = sum(DS.cum.rej1)/nSim
    cum.pow[2] = sum(DS.cum.rej2)/nSim
  } else {
    cum.pow=gsd.power(z = comb.z, bd.z=bd.z)
  }
  
  o = list()
  o$cum.pow = cum.pow
  if (method == "Independent Incremental" || method == "Mixture") {
    o$bd.z = bd.z
  } else if (method == "Disjoint Subjects"){
    #only output 100 rows for rejection boundaries
    o$bd.z = DS.bd.z[1:min(nrow(DS.bd.z), 100),]
  }
  o$multiplicity.method = multiplicity.method
  o$method = method
  
  selection = rep(NA, n.arms-1)
  for (j in 1:(n.arms-1)) {
    selection[j] = sum(s == j) / nSim
  }
  
  o$selection = selection
  o$ssr = sum(new.ss) / nSim
  #o$s=s
  #Calculate the generalized power by simulation, defined as the correct selection of the best dose in OS and H0 rejected
  
  #Best dose by design
  doses.m = m[1:(n.arms-1)]
  max.m = max(doses.m)
  
  if (sum(m == max.m) == 1) {
    #There is a best dose in OS.
    best.dose = (1:(n.arms-1))[doses.m == max.m]
    if (sum(s == best.dose) > 0) {
      correct.selection = (1:nSim)[s == best.dose]
      generalized.pow = rep(NA, K)
    
      if (method == "Disjoint Subjects") {
        correct.DS.cum.rej1 = DS.cum.rej1[correct.selection]
        correct.DS.cum.rej2 = DS.cum.rej2[correct.selection]
      
        generalized.pow[1] = sum(correct.DS.cum.rej1)/nSim
        generalized.pow[2] = sum(correct.DS.cum.rej2)/nSim
      } else {
        correct.comb.z = comb.z[correct.selection, ]
        generalized.pow=gsd.power(z = correct.comb.z, bd.z=bd.z) * length(correct.selection) / nSim
      }
    } else {generalized.pow = 0}
    
    o$best.dose = best.dose 
    o$generalized.pow = generalized.pow
  }
  
  return(o)
}


