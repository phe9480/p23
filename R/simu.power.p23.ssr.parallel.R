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
#' @param nCore Number of cores distributed for simulation;
#' @param seed An integer, or nCore number of integers as random seed for reproducibility;
#' 
#' 
#' @return An object with values:
#' \describe{
#' \item{bd.z}{z value rejection boundary at each analysis}
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
#' simu.power.p23.ssr(nSim=10, n1 = rep(50, 4), n2 = rep(200, 2), m = c(9, 9, 9, 9), 
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
#' simu.power.p23.ssr(nSim=10, n1 = rep(50, 4), n2 = rep(200, 4), m = c(9, 9, 9, 9), 
#' orr = c(0.25, 0.3, 0.3, 0.2), rho = 0.7, dose_selection_endpoint = "ORR",
#' ssr_HR_threshold = 0.8, events_increase = 30,
#' Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, 
#' A1 = 12,Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, 
#' A2 = 12,enrollment.hold=4, DCO1 = 16, targetEvents2=c(300, 380), 
#' sf=gsDesign::sfLDOF, alpha=0.025, multiplicity.method = "dunnett", 
#' method = "Disjoint Subjects")
#' 
#' @importFrom gsDesign gsDesign
#' @importFrom stats uniroot
#' @export 
#' 
simu.power.p23.ssr.parallel <- function(nSim=10, n1 = rep(50, 4), n2 = rep(200, 2), m = c(9,9, 9, 9), 
                              orr = c(0.25, 0.3, 0.4, 0.2), rho = 0.7, dose_selection_endpoint = "ORR",
                              ssr_HR_threshold = 0.8, events_increase = 30, 
                              Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A1 = 12,
                              Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A2 = 12,
                              enrollment.hold=4, DCO1 = 16, targetEvents2=c(300, 380), 
                              alpha=0.025, sf=gsDesign::sfLDOF, multiplicity.method="simes",
                              method = "Independent Incremental", nCore=NULL, seed=123){
  
  simu.power.p23.ssr.onecore <- function(nSim=10, n1 = rep(50, 4), n2 = rep(200, 2), m = c(9,9, 9, 9), 
                                        orr = c(0.25, 0.3, 0.4, 0.2), rho = 0.7, dose_selection_endpoint = "ORR",
                                        ssr_HR_threshold = 0.8, events_increase = 30, 
                                        Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A1 = 12,
                                        Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A2 = 12,
                                        enrollment.hold=4, DCO1 = 16, targetEvents2=c(300, 380), 
                                        alpha=0.025, sf=gsDesign::sfLDOF, multiplicity.method="simes",
                                        method = "Independent Incremental", bd.z=bd.z){
    #Number of analyses in stage 2
    K = length(targetEvents2)
    if (K != 2) {stop; return("ERROR: Only 2 analyses in stage 2 are considered in SSR. The length of targetEvents2 must be 2.")}
    
    #Number of arms
    n.arms = length(n1)
    
    
    #Combination Z values
    comb.z = matrix(NA, nrow=nSim, ncol=K)
    s = rep(NA, nSim) #selected dose
    new.ss = rep(NA, nSim) #indicator whether sample size increased
    
    n2 = c(rep(n2[1], n.arms-1), n2[2])
    
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
      
      #Find DCO to achieve original targetEvents events
      f.u = function(u){
        e1 = fe(DCO = u, r = (n1[s[i]] + n2[s[i]])/(n1[n.arms] + n2[n.arms]), 
                h0 = function(t){log(2)/m[n.arms]}, S0 = function(t){exp(-log(2)/m[n.arms]*t)}, 
                h1 = function(t){log(2)/m[s[i]]}, S1 = function(t){exp(-log(2)/m[s[i]]*t)}, 
                Lambda = Lambda1, n = sum(n1[c(s[i], n.arms)]))$e
        e2 = fe(DCO = max(u-enrollment.hold-A1, 0), r = (n1[s[i]] + n2[s[i]])/(n1[n.arms] + n2[n.arms]), 
                h0 = function(t){log(2)/m[n.arms]}, S0 = function(t){exp(-log(2)/m[n.arms]*t)}, 
                h1 = function(t){log(2)/m[s[i]]}, S1 = function(t){exp(-log(2)/m[s[i]]*t)}, 
                Lambda = Lambda2, n = sum(n2[c(s[i], n.arms)]))$e
        return(e1+e2 - targetEvents2[2])
      }
      DCO.FA0 = uniroot(f.u, c(16, 1000))$root
      
      #Expected events for Stage 1 at original FA: e1.FA0
      e1.FA0=fe(DCO = DCO.FA0, r = (n1[s[i]] + n2[s[i]])/(n1[n.arms] + n2[n.arms]), 
                h0 = function(t){log(2)/m[n.arms]}, S0 = function(t){exp(-log(2)/m[n.arms]*t)}, 
                h1 = function(t){log(2)/m[s[i]]}, S1 = function(t){exp(-log(2)/m[s[i]]*t)}, 
                Lambda = Lambda1, n = sum(n1[c(s[i], n.arms)]))$e
      
      o=conduct.p23.ssr(data=p23i, DCO1=DCO1, dose_selection_endpoint = dose_selection_endpoint, 
                        targetEvents.FA = ssr$targetEvents.FA, e1.FA = e1.FA0,
                        targetEvents2 = targetEvents2, method = method, multiplicity.method=multiplicity.method)
      
      comb.z[i, ] = o$z.tilde
      new.ss[i] = ssr$ssr
    }

  
    # selection = rep(NA, n.arms-1)
    # for (j in 1:(n.arms-1)) {
    #   selection[j] = sum(s == j) 
    # }
    # o$selection = selection
    # o$s=s
    
    re = list(comb.z=comb.z, s=s, new.ss=new.ss)
    
    return(re)
  }
  
  
  ## Start simulation 
  if(is.null(nCore)){
    nCore = max(parallel::detectCores()-4, 8)
  }
  #nsim_per_cluster <- 1000
  cl <- parallel::makeCluster(nCore, type = "PSOCK")
  doParallel::registerDoParallel(cl)
  
  nsim_per_cluster = ceiling(nSim/nCore)
  
  if(length(seed==1)) seed=(1:nCore)*seed
  else if(length(seed)!=nCore) stop("The number of seeds should match the number of cores")
  
  
  #Number of analyses in stage 2
  K = length(targetEvents2)
  
  
  #rejection boundary by GSD. The original GSD boundary is still valid in this SSR setting.
  bd.z = gsDesign::gsDesign(k=K,alpha=alpha,timing=targetEvents2/targetEvents2[K],sfu=sf, test.type=1)$upper$bound
  
  
  # Use parLapply to run in parallel
  results <- parallel::parLapply(cl, seed, fun = simu.power.p23.ssr.onecore,
                                 nSim=nsim_per_cluster,
                                 n1 = n1, n2 = n2, m = m, 
                                 orr = orr, rho = rho, dose_selection_endpoint = dose_selection_endpoint,
                                 Lambda1 = Lambda1, A1 = A1,
                                 Lambda2 = Lambda2, A2 = A2,
                                 enrollment.hold=enrollment.hold, DCO1 = DCO1, targetEvents2=targetEvents2, 
                                 alpha=alpha, sf=sf, multiplicity.method=multiplicity.method,
                                 method = method, bd.z=bd.z
  )
  
  comb.zall <- c()
  s.all <- c()
  newss.all <- c()
  for( i in 1:length(results)){
    comb.zall = rbind(comb.zall, results[[i]]$comb.z)
    # select.all = rbind(select.all, results[[i]]$selection)
    s.all <- c(s.all, results[[i]]$s)
    newss.all <- c(newss.all, results[[i]]$new.ss)
  }
  
  cum.pow=gsd.power(z = comb.zall, bd.z=bd.z)
  
  o = list()
  o$cum.pow = cum.pow
  o$bd.z = bd.z
  
  #o$selection = colSums(select.all)/nSim
  o$multiplicity.method = multiplicity.method
  o$method = method
  #o$s=s.all
  
  n.arms = length(n1)
  selection = rep(NA, n.arms-1)
  for (j in 1:(n.arms-1)) {
    selection[j] = sum(s.all == j) / (nsim_per_cluster * nCore)
  }
  o$selection = selection
  
  o$ssr = sum(newss.all) / (nsim_per_cluster * nCore)
  
  # Stop the cluster
  parallel::stopCluster(cl)
  
  return(o)
  
  
 
}

