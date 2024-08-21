#' Conduct a phase 2/3 dose optimization trial for survival endpoint
#'
#' This function calculates the stage 2 z statistics according to preplanned target number of events at each analysis.
#' The function returns (1) Z statistics at Stage 2 based on prespecified target events for pooled stage 1 (selected dose) and stage 2 patients. (2) The function
#' also returns the weights for combination p value approach in the next step and the weights are sqrt of information fractions. 
#'
#' @param data Dataset produced by simu.p23trial() function
#' @param targetEvents2 Planned target number of events for Stage 2. Either targetEvents2 must be provided. 
#' @param method "Independent Incremental", "Disjoint Subjects". Currently, only "Independent Incremental" method is implemented.
#' 
#' @return 
#' \describe{
#' \item{z1}{Stage 1 z values for each dose}
#' \item{z.c}{z value for each analysis at Stage 2 for the selected dose}
#' \item{z2}{Incremental z statistic at each anaysis at Stage 2, calculated from dose selection to the analysis}
#' \item{w}{Weight, calculated as sqrt (total events for selected dose + control at stage 1 / targetEvents at each analysis at Stage 2 for the selected dose + control combining stage 1 and stage 2 patients.)}
#' \item{selected.dose}{Selected dose}
#' \item{example.data}{A simulated dataset example}
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
#' Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A1 = 12,
#' Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A2 = 12,
#' enrollment.hold=4)
#' 
#' sel = select.dose (data=p23trial, DCO1=16)
#' 
#' o=conduct.p23(data=p23trial, DCO1=16, targetEvents = c(300, 380), method = "Independent Incremental")
#' 
#' 
#' 
#' @export 
#' 
conduct.p23 = function(data=p23trial, DCO1=16, targetEvents = c(300, 380), method = "Independent Incremental"){

  #1. Dose selection  
  sel = select.dose.p23 (data=p23trial, DCO1=16)
  z1 = sel$z1; e1=sel$e1; s=sel$s
  
  #2. Assemble the trial data combining stage 1 and stage 2 for selected dose + control
  dat23 = data[data$group == 0 | data$group == s, ]
  
  #3. Determine z2
  K = length(targetEvents2) #K analyses in Stage 2
  z.c = rep(NA, ncol=K) #z value at each of Stage 2 analyses using combined data
  z2  = rep(NA, ncol=K) # incremental z value from stage 2 compared to Stage 1
  w = rep(NA, ncol=K) #weights as sqrt(information fraction) = sqrt(events at dose selection / events at each stage 2 analysis) at each analysis.

  if (method == "Independent Incremental"){
    for (k in 1:K){
      #4. data cut for each analysis
      dat23k = f.dataCut(data=dat23, targetEvents=targetEvents2[k])
      
      #5. Calculate z statistics for logrank test combining Stage 1 and stage 2
      z.c[k] = logrank.one.sided(time=dat23k$survTimeCut, cnsr=dat23k$cnsrCut, group=dat23k$group)$z
      frac.k = e1[s]/targetEvents2[k]  
      z2[k] = (z.c[k] - sqrt(frac.k)*z1[s])/sqrt(1-frac.k)
      w[k] = sqrt(frac.k)
    }
    
  }
  
  o = list()
  o$z1 = matrix(z1, nrow=1)
  o$e1 = matrix(e1, nrow=1)
  o$s = s
  o$z2 = matrix(z2, nrow=1)
  o$z.c = matrix(z.c, nrow=1)
  o$w = matrix(w, nrow=1)
  
  return(o)
}

  

