#' Conduct a phase 2/3 dose optimization trial for survival endpoint
#'
#' This function calculates the stage 2 z statistics according to preplanned target number of events at each analysis.
#' The function returns (1) Z statistics at Stage 2 based on prespecified target events for pooled stage 1 (selected dose) and stage 2 patients. (2) The function
#' also returns the weights for combination p value approach in the next step and the weights are sqrt of information fractions. 
#'
#' @param data Dataset produced by simu.p23trial function
#' @param DCO1 Data cutoff for Stage 1
#' @param targetEvents2 Planned target number of events for Stage 2. Either targetEvents2 must be provided.
#' @param dose_selection_endpoint Endpoint for dose selection at Stage 1. "ORR" or "not ORR"
#' @param targetEvents.FA Target number of events after SSR
#' @param e1.FA Expected number of events at FA for stage 1 subjects
#' @param method Options include "Independent Incremental": z1 at dose selection and z2 is from dose selection to kth analysis at stage 2; 
#' "Disjoint Subjects": z1 is at kth analysis for stage 1 subjects; z2 is at the kth analysis for stage 2 subjects. z1 will be adjusted by multiplicity and closed testing procedure at each analysis.
#' "Mixture": Only consider disjoint subjects at first analysis in stage 2. Starting from the 2nd analysis, consider independent incremental methods. Only z1 at 1st analysis will be adjusted by multiplicity and closed testing procedure.
#' @param multiplicity.method "simes", "Dunnett". 
#' 
#' @return If method is "Independent Incremental", return an object with variables
#' \describe{
#' \item{s}{Selected dose}
#' \item{z1}{Stage 1 z values for each dose}
#' \item{z2}{Incremental z statistic at each anaysis at Stage 2, calculated from dose selection to the analysis}
#' \item{w}{Weight, calculated as sqrt (total events for selected dose + control at stage 1 / targetEvents at each analysis at Stage 2 for the selected dose + control combining stage 1 and stage 2 patients.)}
#' \item{z.c}{z value for each analysis at Stage 2 for the selected dose}
#' }
#' If method is "Disjoint Subjects", return an object with variables
#' #' \describe{
#' \item{s}{Selected dose}
#' \item{z1}{Stage 1 z values for unselected doses combined with the selected dose at kth analysis at Stage 2. z1(k, unselected doses Stage 1 followed by selected dose at k)}
#' \item{z2}{Incremental z statistic at each anaysis at Stage 2, calculated from dose selection to the analysis}
#' \item{w}{Weight, calculated as sqrt (total events for selected dose + control at stage 1 / targetEvents at each analysis at Stage 2 for the selected dose + control combining stage 1 and stage 2 patients.)}
#' }
#' If method is "Mixture", return an object with variables
#' \describe{
#' \item{s}{Selected dose}
#' \item{z1.unselected}{Stage 1 z values for unselected doses}
#' \item{z11}{Stage 1 subjects z value at 1st analysis}
#' \item{z21}{Stage 2 subjects z value at 1st analysis}
#' \item{omega1}{Weight for 1st analysis combining z11 and z21 to calcualte z1.tilde}
#' \item{z.tilde}{Test statistics for group sequential design at Stage 2}
#' \item{z.c}{z value for each analysis at Stage 2 for the selected dose}
#' \item{z2}{Incremental z statistic at each anaysis at Stage 2, calculated from 1st analysis to each following analysis}
#' \item{v}{Weight for 2nd and further analysis combining z1.tilde and z2(k) to calculate z.tilde(k)}
#' }
#' 
#' @examples
#' #Example (1): Stage 1: 4 arms; 3 dose levels; each arm 50 patients.
#' #Stage 2: additional 200 patients per arm will be enrolled at stage 2
#' #medians for the 4 arms: 9, 11, 13 and control = 8 months
#' #Enrollment: 12 months uniform in stage 1; 12 months uniform in stage 2
#' #Holding period: 4 months between stage 1 and 2
#' #Dose selection will be based on data cut at 16 months
#' #Stage 2 has 2 planned analyses at 300 and 380 events respectively.
#'  
#' p23trial = simu.p23trial(n1 = rep(50, 4), n2 = rep(200, 4), m = c(9,9, 9, 9), 
#' orr = NULL, rho = NULL, dose_selection_endpoint = "not ORR",
#' Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A1 = 12,
#' Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A2 = 12,
#' enrollment.hold=4)
#' 
#' select.dose.p23 (data=p23trial, DCO1=16, dose_selection_endpoint = "not ORR")
#' 
#' #Independent incremental; dose not selected by ORR
#' conduct.p23.ssr(data=p23trial, DCO1=16, dose_selection_endpoint = "not ORR",
#' targetEvents.FA = 410,e1.FA = 93,
#' targetEvents2 = c(300, 380), method = "Independent Incremental")
#' 
#' #Disjoint Subjects; dose not selected by ORR
#' conduct.p23.ssr(data=p23trial, DCO1=16, dose_selection_endpoint = "not ORR",
#' targetEvents.FA = 410,e1.FA = 93,
#' targetEvents2 = c(300, 380), method = "Disjoint Subjects")
#' 
#' #Mixture; dose not selected by ORR
#' conduct.p23.ssr(data=p23trial, DCO1=16, dose_selection_endpoint = "not ORR",
#' targetEvents.FA = 410,e1.FA = 93,
#' targetEvents2 = c(300, 380), method = "Mixture")
#' 
#' #Example (2): Stage 1: 4 arms; 3 dose levels; each arm 50 patients.
#' #Stage 2: additional 200 patients per arm will be enrolled at stage 2
#' #medians for the 4 arms: 9, 11, 13 and control = 8 months
#' #Enrollment: 12 months uniform in stage 1; 12 months uniform in stage 2
#' #Holding period: 4 months between stage 1 and 2
#' #Dose selection will be based on ORR with data cut at 16 months
#' #Stage 2 has 2 planned analyses at 300 and 380 events respectively.
#'
#' 
#' p23trial = simu.p23trial(n1 = rep(50, 4), n2 = rep(200, 4), m = c(9,9, 9, 9), 
#' orr = c(0.25, 0.3, 0.2, 0.2), rho = 0.7, dose_selection_endpoint = "ORR",
#' Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A1 = 12,
#' Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A2 = 12,
#' enrollment.hold=4)
#' 
#' sel=select.dose.p23 (data=p23trial, DCO1=16, dose_selection_endpoint = "ORR")
#' 
#' #Independent Incremental, dose selected by ORR
#' conduct.p23.ssr(data=p23trial, DCO1=16, dose_selection_endpoint = "ORR",
#' targetEvents.FA = 410,
#' targetEvents2 = c(300, 380), method = "Independent Incremental")
#' 
#' #Disjoint Subjects; dose not selected by ORR
#' d = conduct.p23.ssr(data=p23trial, DCO1=16, dose_selection_endpoint = "ORR",
#' targetEvents.FA = 410,e1.FA = 93,
#' targetEvents2 = c(300, 380), method = "Disjoint Subjects")
#' 
#' #Mixture; dose not selected by ORR
#' conduct.p23.ssr(data=p23trial, DCO1=16, dose_selection_endpoint = "ORR",
#' targetEvents.FA = 410,e1.FA = 93,
#' targetEvents2 = c(300, 380), method = "Mixture")
#' 
#' 
#' @importFrom survival coxph Surv
#' @export 
#' 
conduct.p23.ssr = function(data=NULL, DCO1=16, targetEvents2 = c(300, 380), 
                           dose_selection_endpoint = "ORR",
                           targetEvents.FA = 410, e1.FA = 93,
                           method = "Independent Incremental", 
                           multiplicity.method="simes"){

  #1. Dose selection  
  sel = select.dose.p23 (data=data, DCO1=DCO1, dose_selection_endpoint = dose_selection_endpoint)
  s=sel$s
  
  #2. Assemble the trial data combining stage 1 and stage 2 for selected dose + control
  dat23 = data[data$group == 0 | data$group == s, ]
  
  K = length(targetEvents2) #K analyses in Stage 2
  if (K != 2) {stop; return("ERROR: Only 2 analyses in stage 2 are considered in SSR. The length of targetEvents2 must be 2.")}
  
  n.arms = length(unique(data$group))
  
  if (method == "Independent Incremental"){
    #1st component in weighted z as of dose selection IAd.
    #2nd component in weighted z as of incremental from IAd to IA
    
    #At IAd
    z1 = sel$z1; e1=sel$e1; 
    
     #At IA
    
      #4. data cut for each analysis
      dat23.1 = f.dataCut(data=dat23, targetEvents=targetEvents2[1])
      
      #5. Calculate z statistics for logrank test combining Stage 1 and stage 2
      z.c1 = logrank.one.sided(time=dat23.1$survTimeCut, cnsr=dat23.1$cnsrCut, group=dat23.1$group)$z
      
      #2nd component z statistic and weight for combining z1 and z2 in next step
      frac.1 = e1[s]/targetEvents2[1]  
      z2.1 = (z.c1[1] - sqrt(frac.1)*z1[s])/sqrt(1-frac.1)
      w1 = sqrt(frac.1)
      
      #calculate Z.tilde

      #adjusted z1 by CTP and simes method
      z1.tilde = comb.pvalue.p23(z1=matrix(z1, nrow=1),  z2 = z2.1, selected.dose = s, bd.z=Inf, w=w1, method=multiplicity.method)$comb.z

    #At FA
      #weight is fixed according to the original target events at FA
      #the weight is for constructing FA test statistic z2.tilde = w2*z1.tilde + sqrt(1-w2^2)*z2
      w2 = sqrt(targetEvents2[1]/targetEvents2[2])
      
      dat23.FA = f.dataCut(data=dat23, targetEvents=targetEvents.FA)
      z.FA = logrank.one.sided(time=dat23.FA$survTimeCut, cnsr=dat23.FA$cnsrCut, group=dat23.FA$group)$z
      #z2: independent incremental
      frac.new = targetEvents2[1]/targetEvents.FA
      
      z2 = (z.FA - sqrt(frac.new)*z.c1)/sqrt(1-frac.new)
      z2.tilde = w2 * z1.tilde + sqrt(1-w2^2)*z2
    
  } else if (method == "Disjoint Subjects") {
    #1st component z1 in weighted z as of Stage 1 subjects at each analysis eg z11
    #For each analysis, z1 has unselected doses at IAd for multiplicity adjustment and z for stage 1 subjects
    #2nd component in weighted z as of Stage 2 subjects at each analysis eg z12

    #At IAd
    #multiplicity adjustment needs z from IAd for the unselected doses
    z1.IAd = sel$z1[-s]
    
    #At IA
      #4. data cut
      dat23.1 = f.dataCut(data=dat23, targetEvents=targetEvents2[1])
      
      #stage 1 subjects
      dat23.11 = dat23.1[dat23.1$stage == 1, ] 
      
      #z statistic for stage 1 subjects at analysis 1
      z11 = logrank.one.sided(time=dat23.11$survTimeCut, cnsr=dat23.11$cnsrCut, group=dat23.11$group)$z
      
      #stage 2 subjects at analysis 1
      dat23.12 = dat23.1[dat23.1$stage == 2, ] 
      
      #z statistic for stage 2 subjects at analysis 1
      z12 = logrank.one.sided(time=dat23.12$survTimeCut, cnsr=dat23.12$cnsrCut, group=dat23.12$group)$z
      
      #Number of events for stage 1 subjects
      e11 = sum(1-dat23.11$cnsrCut)
      frac.11 = e11/targetEvents2[1]
      
      #weight w for combining z1 and z2 in next step
      w11 = sqrt(frac.11)
      
      #adjusted z11 by CTP and simes method (trick the program with selected dose = last one, because z11 is the last)
      z1.tilde = comb.pvalue.p23(z1=matrix(c(z1.IAd, z11), nrow=1),  z2 = z12, 
                                 selected.dose = length(c(z1.IAd, z11)), bd.z=Inf, 
                                 w=w11, method=multiplicity.method)$comb.z

    #At FA
      #weight is fixed according to the original target events at FA and 
      #originally expected events of stage 1 subjects at FA
      #the weight is for constructing FA test statistic z2.tilde = w2*z1.tilde + sqrt(1-w2^2)*z2
      #e1.FA: expected number of events for stage 1 subjects at analysis 2 for the original targetEvents2 at FA
      
      ####The following provides an estimate based on a simulated dataset for programming purpose,
      #but e1.FA should be provided as a parameter entry, calculated as the expected events
      
      #dat23.FA.old = f.dataCut(data=dat23, targetEvents=targetEvents2[2])
      #dat23.FA.old.1 = dat23.FA.old[dat23.FA.old$stage == 1, ]
      #e1.FA = sum(1-dat23.FA.old.1$cnsrCut)
      w2 = sqrt(e1.FA/targetEvents2[2])
      
      dat23.FA = f.dataCut(data=dat23, targetEvents=targetEvents.FA)
      dat23.1FA = dat23.FA[dat23.FA$stage==1, ]
      dat23.2FA = dat23.FA[dat23.FA$stage==2, ]
      
      z.1FA = logrank.one.sided(time=dat23.1FA$survTimeCut, cnsr=dat23.1FA$cnsrCut, group=dat23.1FA$group)$z
      z.2FA = logrank.one.sided(time=dat23.2FA$survTimeCut, cnsr=dat23.2FA$cnsrCut, group=dat23.2FA$group)$z
      
      #multiplicity adjustment for z.1FA by IAd
      #adjusted z11 by CTP and simes method (trick the program with selected dose = last one, because z11 is the last)
      z2.tilde = comb.pvalue.p23(z1=matrix(c(z1.IAd, z.1FA), nrow=1),  z2 = z.2FA, 
                                 selected.dose = length(c(z1.IAd, z11)), bd.z=Inf, 
                                 w=w2, method=multiplicity.method)$comb.z
      
  } else if (method == "Mixture") {
    #Only adjust multiplicity for Stage 1 patients at 1st analysis, 
    #then construct following z as independent incrementals
    
    #Multiplicity adjustment needs z from IAd for the unselected doses
    z1.IAd = sel$z1[-s]
    
    #At IA
    ##################################
    
    #data cut
    dat23.1 = f.dataCut(data=dat23, targetEvents=targetEvents2[1])
    
    #stage 1 subjects
    dat23.11 = dat23.1[dat23.1$stage == 1, ]
    
    #z statistic for stage 1 subjects
    z11 = logrank.one.sided(time=dat23.11$survTimeCut, cnsr=dat23.11$cnsrCut, group=dat23.11$group)$z
    
    #stage 2 subjects
    dat23.12 = dat23.1[dat23.1$stage == 2, ]
    
    #z statistic for stage 2 subjects
    z12 = logrank.one.sided(time=dat23.12$survTimeCut, cnsr=dat23.12$cnsrCut, group=dat23.12$group)$z
    
    #Number of events for stage 1 subjects
    e11 = sum(1-dat23.11$cnsrCut)
    frac1 = e11/targetEvents2[1]
    
    #weight w for combining z1 and z2 in next step
    omega1 = sqrt(frac1)
    
    #adjusted z1 by CTP and simes method
    z1.tilde = comb.pvalue.p23(z1=matrix(c(z1.IAd, z11), nrow=1),  z2 = z12, selected.dose = length(c(z1.IAd, z11)), bd.z=Inf, w=omega1, method=multiplicity.method)$comb.z

    #At FA (Independent Incremental)
    ##################################
    
    #z statistic at IA
    z.c1 = logrank.one.sided(time=dat23.1$survTimeCut, cnsr=dat23.1$cnsrCut, group=dat23.1$group)$z
    
    #4. data cut 
    dat23.FA = f.dataCut(data=dat23, targetEvents=targetEvents.FA)
        
    #5. Calculate z statistics for logrank test combining Stage 1 and stage 2
    z.c2 = logrank.one.sided(time=dat23.FA$survTimeCut, cnsr=dat23.FA$cnsrCut, group=dat23.FA$group)$z
      
    #2nd component z statistic and weight for combining z1 and z2 in next step
    #weight must be prefixed in SSR, here according to the original targetEvents2
    frac = targetEvents2[1]/targetEvents2[2] #prefixed
    
    frac.new = targetEvents2[1]/targetEvents.FA
    #incremental z from z.1c to z.2c
    z2 = (z.c2 - sqrt(frac.new)*z.c1)/sqrt(1-frac.new)
      
    z2.tilde = sqrt(frac) * z1.tilde + sqrt(1-frac) * z2

  }
  
  o = list()
  o$s = s
  
  z.tilde = c(z1.tilde, z2.tilde)
  o$z.tilde = z.tilde

  if (dose_selection_endpoint == "ORR") {
    o$orr.diff = matrix(sel$orr.diff, nrow=1)
  }
  o$dose.selection.endpoint=dose_selection_endpoint
  o$method = method
  o$multiplicity.method = multiplicity.method
  
  return(o)
}

  

