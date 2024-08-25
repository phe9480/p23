#' Perform one-sided logrank test
#'
#' This functions performs the one-sided logrank test. The standard logrank test in survival package only produces two-sided test. This function can facilitate one-sided logrank test.
#'
#' @param time Survival time
#' @param cnsr Censoring status (0 = event, 1 = censor)
#' @param group Group indicator (0 = control, 1 = experimental arm)
#' 
#' @return
#' \describe{
#' \item{z}{Test statistics z value}
#' \item{chisq}{Test statistics chisq value}
#' \item{two.sided.p}{two.sided p value}
#' \item{one.sided.p}{one sided p value}
#' }
#' 
#' @examples
#' 
#' n = 100
#' time = c(rexp(n, rate=log(2)/12), rexp(n, rate=log(2)/12*0.7))
#' cnsr = c(as.numeric(runif(n) > 0.7), as.numeric(runif(n) > 0.7))
#' group = c(rep(0, n), rep(1, n))
#' 
#' data = as.data.frame(cbind(time, cnsr, group))
#' 
#' logrank.one.sided(time, cnsr, group)
#' 
#' @importFrom survival survdiff
#' @export 
#' 
logrank.one.sided = function(time, cnsr, group){

    lr.test = survival::survdiff(survival::Surv(time, 1-cnsr) ~ group)
    
    #convert to z value in correct direction: z>0 means better experimental arm.
    better = as.numeric(lr.test$obs[2] < lr.test$exp[2])
    sign = 2*better - 1
    z = sqrt(lr.test$chisq) * sign
    
    o = list()
    o$z = z
    o$chisq = lr.test$chisq
    o$two.sided.p = lr.test$pvalue
    o$one.sided.p = 1-pnorm(z)
    
    return(o)
}

