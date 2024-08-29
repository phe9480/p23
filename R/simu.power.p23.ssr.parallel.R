#' #' Power Calculation by Simulations for Phase 2/3 Inferential Seamless Design with Sample Size Adjustment
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
#' simu.power.p23.ssr.parallel(nSim=10, n1 = rep(50, 4), n2 = rep(200, 2), m = c(9, 9, 9, 9), 
#' orr = NULL, rho = NULL, dose_selection_endpoint = "not ORR",
#' ssr_HR_threshold = 0.8, events_increase = 30,
#' Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A1 = 12,
#' Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A2 = 12,
#' enrollment.hold=4, DCO1 = 16, targetEvents2=c(300, 380), sf=gsDesign::sfLDOF, 
#' alpha=0.025, method = "Independent Incremental")
#' 
#' #Example (2): #Dose selection decision based on ORR
#' simu.power.p23.ssr.parallel(nSim=10, n1 = rep(50, 4), n2 = rep(200, 4), m = c(9, 9, 9, 9), 
#' orr = c(0.25, 0.3, 0.4, 0.2), rho = 0.7, dose_selection_endpoint = "ORR",
#' ssr_HR_threshold = 0.8, events_increase = 30,
#' Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A1 = 12,
#' Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A2 = 12,
#' enrollment.hold=4, DCO1 = 16, targetEvents2=c(300, 380), sf=gsDesign::sfLDOF, 
#' alpha=0.025, method = "Independent Incremental")
#' 
#' simu.power.p23.ssr.parallel(nSim=10, n1 = rep(50, 4), n2 = rep(200, 4), m = c(9, 9, 9, 9), 
#' orr = c(0.25, 0.3, 0.3, 0.2), rho = 0.7, dose_selection_endpoint = "ORR",
#' ssr_HR_threshold = 0.8, events_increase = 30,
#' Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A1 = 12,
#' Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A2 = 12,
#' enrollment.hold=4, DCO1 = 16, targetEvents2=c(300, 380), sf=gsDesign::sfLDOF, 
#' alpha=0.025, method = "Disjoint Subjects")
#' 
#' simu.power.p23.ssr.parallel(nSim=10, n1 = rep(50, 4), n2 = rep(200, 4), m = c(9, 9, 9, 9), 
#' orr = c(0.25, 0.3, 0.3, 0.2), rho = 0.7, dose_selection_endpoint = "ORR",
#' ssr_HR_threshold = 0.8, events_increase = 30,
#' Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A1 = 12,
#' Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A2 = 12,
#' enrollment.hold=4, DCO1 = 16, targetEvents2=c(300, 380), sf=gsDesign::sfLDOF, 
#' alpha=0.025, method = "Mixture")
#' 
#' simu.power.p23.ssr.parallel(nSim=10, n1 = rep(50, 4), n2 = rep(200, 4), m = c(9, 9, 9, 9), 
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
                                        orr = c(0.25, 0.3, 0.4, 0.2), rho = 0.7, 
                                        dose_selection_endpoint = "ORR",
                                        ssr_HR_threshold = 0.8, events_increase = 30, 
                                        Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A1 = 12,
                                        Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A2 = 12,
                                        enrollment.hold=4, DCO1 = 16, targetEvents2=c(300, 380), 
                                        alpha=0.025, sf=gsDesign::sfLDOF, multiplicity.method="simes",
                                        method = "Independent Incremental", nCore=NULL, seeds=123){
  
  simu.power.p23.ssr.onecore <- function(seed=123, nSim=10, n1 = rep(50, 4), n2 = rep(200, 2), m = c(9,9, 9, 9), 
                                         orr = c(0.25, 0.3, 0.4, 0.2), rho = 0.7, 
                                         dose_selection_endpoint = "ORR",
                                         ssr_HR_threshold = 0.8, events_increase = 30, 
                                         Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A1 = 12,
                                         Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A2 = 12,
                                         enrollment.hold=4, DCO1 = 16, targetEvents2=c(300, 380), 
                                         alpha=0.025, sf=gsDesign::sfLDOF, multiplicity.method="simes",
                                         method = "Independent Incremental"){
    
    set.seed(seed)
    
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
    
    re = list(comb.z=comb.z, s=s, new.ss=new.ss,
              DS.cum.rej1=DS.cum.rej1, DS.cum.rej2=DS.cum.rej2,
              DS.bd.z=DS.bd.z)
    
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
  
  if(length(seeds==1)) seeds=(1:nCore)*seeds
  else if(length(seeds)!=nCore) stop("The number of seeds should match the number of cores")
  
  
  #Number of analyses in stage 2
  K = length(targetEvents2)
  
  
  #rejection boundary by GSD. The original GSD boundary is still valid in this SSR setting.
  bd.z = gsDesign::gsDesign(k=K,alpha=alpha,timing=targetEvents2/targetEvents2[K],sfu=sf, test.type=1)$upper$bound
  
  n.arms = length(n1)
  
  # Use parLapply to run in parallel
  results <- parallel::parLapply(cl, seeds, fun = simu.power.p23.ssr.onecore,
                                 nSim=nsim_per_cluster,
                                 n1 = n1, n2 = n2, m = m, 
                                 orr = orr, rho = rho, dose_selection_endpoint = dose_selection_endpoint,
                                 Lambda1 = Lambda1, A1 = A1,
                                 Lambda2 = Lambda2, A2 = A2,
                                 enrollment.hold=enrollment.hold, DCO1 = DCO1, targetEvents2=targetEvents2, 
                                 alpha=alpha, sf=sf, multiplicity.method=multiplicity.method,
                                 method = method
  )
  
  comb.zall <- c()
  s.all <- c()
  newss.all <- c()
  DS.cum.rej1.all <- DS.cum.rej2.all <- c()
  DS.bd.z.all <- c()
  for( i in 1:length(results)){
    comb.zall = rbind(comb.zall, results[[i]]$comb.z)
    # select.all = rbind(select.all, results[[i]]$selection)
    s.all <- c(s.all, results[[i]]$s)
    newss.all <- c(newss.all, results[[i]]$new.ss)
    DS.cum.rej1.all <- c(DS.cum.rej1.all, results[[i]]$DS.cum.rej1)
    DS.cum.rej2.all <- c(DS.cum.rej2.all, results[[i]]$DS.cum.rej2)
    DS.bd.z.all <- rbind(DS.bd.z.all, results[[i]]$DS.bd.z)
  }
  
  nSim_actual = (nsim_per_cluster * nCore)
  cum.pow = rep(NA, K)
  if (method == "Disjoint Subjects") {
    cum.pow[1] = sum(DS.cum.rej1.all)/nSim_actual
    cum.pow[2] = sum(DS.cum.rej2.all)/nSim_actual
  } else {
    cum.pow=gsd.power(z = comb.zall, bd.z=bd.z)
  }
  
  o = list()
  o$cum.pow = cum.pow
  if (method == "Independent Incremental" || method == "Mixture") {
    o$bd.z = bd.z
  } else if (method == "Disjoint Subjects"){
    #only output 100 rows for rejection boundaries
    o$bd.z = DS.bd.z.all[1:min(nrow(DS.bd.z.all), 100),]
  }
  o$multiplicity.method = multiplicity.method
  o$method = method
  
  selection = rep(NA, n.arms-1)
  for (j in 1:(n.arms-1)) {
    selection[j] = sum(s.all == j) / nSim_actual
  }
  
  o$ssr = sum(newss.all) / nSim_actual
  #o$s=s
  o$ssr.all = newss.all

  
  # Stop the cluster
  parallel::stopCluster(cl)
  
  return(o)
  
  
 
}

