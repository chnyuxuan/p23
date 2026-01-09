#' Phase 2/3 Inferential Seamless Design Power calculation by simulations
#'
#' This functions calculates the cumulative power and overall power for a group sequential design by simulations.
#'
#' @param nSim Number of simulated trials
#' @param n1 Stage 1 sample size of each dose and control arm. length(n1) must be number of arms. The last component is control.
#' @param n2 Stage 2 Sample size of the selected dose and control arm. length(n2) must be 2. The second component is control.
#' @param m Median survival time for each arm (dose 1, dose 2, ..., control). length(m) must be equal to length(n1). The last component in m is control.
#' @param orr ORR for each arm. length(orr) = length(m). 
#' @param rho Correlation between ORR and time to event endpoint
#' @param dose_selection_endpoint  Dose selection end point: "ORR" or "not ORR"
#' @param A1 Enrollment period for Stage 1
#' @param Lambda1 Enrollment distribution function (CDF) for stage 1.
#' @param DCO1 Data cutoff date for Stage 1
#' @param enrollment.hold Holding period in months after DCO1 of Stage 1 prior to enrollment of Stage 2 patients. 0 means seamless enrollment.
#' @param A2 Enrollment period for Stage 2
#' @param Lambda2 Enrollment distribution function (CDF) for stage 2.
#' @param targetEvents2 Planned target number of events for Stage 2. Either targetEvents2 must be provided. 
#' @param alpha Type I error (one-sided) for testing the selected dose, usually 0.025.
#' @param sf Spending functions. acceptable options include all spending functions in gsDesign R package, for example, "gsDesign::sfLDOF"
#' @param multiplicity.method "simes" or "Dunnett"
#' @param method Options include "Independent Incremental": z1 at dose selection and z2 is from dose selection to kth analysis at stage 2; 
#' "Disjoint Subjects": z1 is at kth analysis for stage 1 subjects; z2 is at the kth analysis for stage 2 subjects. z1 will be adjusted by multiplicity and closed testing procedure at each analysis.
#' "Mixture": Only consider disjoint subjects at first analysis in stage 2. Starting from the 2nd analysis, consider independent incremental methods. Only z1 at 1st analysis will be adjusted by multiplicity and closed testing procedure.
#' @param ORRdiff The tolerated difference between ORR of highest dose and lower doses, with default value 0. 
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
#'
#' #Dose selection decision is NOT based on ORR.
#' simu.power.p23(nSim=10, n1 = rep(50, 4), n2 = rep(200, 2), m = c(9, 9, 9, 9), 
#' orr = NULL, rho = NULL, dose_selection_endpoint = "not ORR",
#' Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A1 = 12,
#' Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A2 = 12,
#' enrollment.hold=4, DCO1 = 16, targetEvents2=c(300, 380), sf=gsDesign::sfLDOF, 
#' alpha=0.025, method = "Independent Incremental")
#' 
#' #Example (2): #Dose selection decision based on ORR
#' simu.power.p23(nSim=10, n1 = rep(50, 4), n2 = rep(200, 4), m = c(9, 9, 9, 9), 
#' orr = c(0.25, 0.3, 0.4, 0.2), rho = 0.7, dose_selection_endpoint = "ORR",
#' Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A1 = 12,
#' Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A2 = 12,
#' enrollment.hold=4, DCO1 = 16, targetEvents2=c(300, 380), sf=gsDesign::sfLDOF, 
#' alpha=0.025, method = "Independent Incremental")
#' 
#' simu.power.p23(nSim=10, n1 = rep(50, 4), n2 = rep(200, 4), m = c(9, 9, 9, 9), 
#' orr = c(0.25, 0.3, 0.3, 0.2), rho = 0.7, dose_selection_endpoint = "ORR",
#' Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A1 = 12,
#' Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A2 = 12,
#' enrollment.hold=4, DCO1 = 16, targetEvents2=c(300, 380), sf=gsDesign::sfLDOF, 
#' alpha=0.025, method = "Disjoint Subjects")
#' 
#' simu.power.p23(nSim=10, n1 = rep(50, 4), n2 = rep(200, 4), m = c(9, 9, 9, 9), 
#' orr = c(0.25, 0.3, 0.3, 0.2), rho = 0.7, dose_selection_endpoint = "ORR",
#' Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A1 = 12,
#' Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A2 = 12,
#' enrollment.hold=4, DCO1 = 16, targetEvents2=c(300, 380), sf=gsDesign::sfLDOF, 
#' alpha=0.025, method = "Mixture")
#' 
#' simu.power.p23(nSim=10, n1 = rep(50, 4), n2 = rep(200, 4), m = c(9, 9, 9, 9), 
#' orr = c(0.25, 0.3, 0.3, 0.2), rho = 0.7, dose_selection_endpoint = "ORR",
#' Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, 
#' A1 = 12,Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, 
#' A2 = 12,enrollment.hold=4, DCO1 = 16, targetEvents2=c(300, 380), sf=gsDesign::sfLDOF, 
#' alpha=0.025, multiplicity.method = "dunnett", method = "Disjoint Subjects")
#' 
#' 
#' simu.power.p23(nSim=1000, n1 = rep(50, 4), n2 = rep(200, 4), m = c(9, 9, 12, 9), 
#' orr = c(0.25, 0.3, 0.3, 0.2), rho = 0.7, dose_selection_endpoint = "ORR",
#' Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, 
#' A1 = 12,Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, 
#' A2 = 12,enrollment.hold=4, DCO1 = 16, targetEvents2=c(300, 380), sf=gsDesign::sfLDOF, 
#' alpha=0.025, multiplicity.method = "dunnett", method = "Disjoint Subjects")
#' 
#' @importFrom gsDesign gsDesign
#' @export 
#' 
simu.power.p23 = function(nSim=10, n1 = rep(50, 4), n2 = rep(200, 2), m = c(9,9, 9, 9), 
                          orr = c(0.25, 0.3, 0.4, 0.2), rho = 0.7, dose_selection_endpoint = "ORR",
                          Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A1 = 12,
                          Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A2 = 12,
                          enrollment.hold=4, DCO1 = 16, targetEvents2=c(300, 380), 
                          e1 = NULL,
                          alpha=0.025, sf=gsDesign::sfLDOF, multiplicity.method="simes",
                          method = "Independent Incremental", ORRdiff = 0
                          ){
  
  
  #Number of analyses in stage 2
  K = length(targetEvents2)
  
  #Number of arms
  n.arms = length(n1)
  
  # #rejection boundary by traditional GSD
  # if (K == 1) {bd.z = qnorm(1-alpha)} else {
  #   bd.z = gsDesign::gsDesign(k=K,alpha=alpha,timing=targetEvents2/targetEvents2[K],sfu=sf, test.type=1)$upper$bound
  # }
  
  #Combination Z values
  comb.z = bd.z = matrix(NA, nrow=nSim, ncol=K)
  s = rep(NA, nSim) #selected dose
  
  n2 = c(rep(n2[1], n.arms-1), n2[2])
  
  # calculate pre-specified weights YC ============================
  if(is.null(e1)){
    e1 = matrix(rep(0,length(targetEvents2)*(n.arms-1)), nrow=length(targetEvents2))
    for(k in 1:length(targetEvents2)){
      e1[k,] = e1.ssr(n1 = n1, n2 = n2, m = m, 
                      A1 = 12, Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, 
                      Lambda2 = function(t){(t/10)*as.numeric(t<= 10) + as.numeric(t > 10)}, 
                      enrollment.hold=enrollment.hold, targetEvents = targetEvents2[k])
    }
  }
  
  for (i in 1:nSim){
    if(i%% 1000 == 0){print(i)}
    p23i = simu.p23trial(n1 = n1, n2 = n2, m = m, 
                             orr = orr, rho = rho, dose_selection_endpoint = dose_selection_endpoint,
                             Lambda1 = Lambda1, A1 = A1, 
                             Lambda2 = Lambda2, A2 = A2, enrollment.hold=enrollment.hold)
    o=conduct.p23(data=p23i, DCO1=DCO1, 
                  dose_selection_endpoint = dose_selection_endpoint, 
                  targetEvents2 = targetEvents2, method = method, 
                  multiplicity.method=multiplicity.method,
                  e1=e1, ORRdiff=ORRdiff)
    s[i] = o$s
    
    if(o$method=="NA"){ # deal with IA exceeds FA YC =============================
      comb.z[i,]=c(NA, o$z) # must be before gsDesign function in calculating the boundary bd.z!
      next
    }
    
    #rejection boundary by traditional GSD
    if (K == 1) {bd.z[i] = qnorm(1-alpha)} else {
      if(o$method == "Disjoint Subjects"){
        corr.z = o$w[1]*o$w[2]*sqrt(o$actualEventsS1[1]/o$actualEventsS1[K]) + 
          sqrt(1-o$w[1]^2)*sqrt(1-o$w[2]^2)*sqrt((o$actualEvents[1]-o$actualEventsS1[1])/(o$actualEvents[K]-o$actualEventsS1[K]))
        
        bd.z[i,] = gsDesign::gsDesign(k=K,alpha=alpha,
                                      sfu=sf, test.type=1,
                                      n.I = c(o$actualEvents[1], o$actualEvents[1]/corr.z^2),
                                      maxn.IPlan = targetEvents2[K])$upper$bound
      }else{
        bd.z[i,] = gsDesign::gsDesign(k=K,alpha=alpha,timing=o$actualEvents/o$actualEvents[K],sfu=sf, test.type=1)$upper$bound
      }
      }
    
    if (o$method == "Independent Incremental") {
      for (j in 1:K){
        oj = comb.pvalue.p23(z1=o$z1,  z2 = o$z2[,j], bd.z=bd.z[i,j], w=o$w[,j], selected.dose = s[i], method=multiplicity.method)
        comb.z[i, j] = oj$comb.z; 
      }
    } else if (o$method == "Disjoint Subjects") {
      for (j in 1:K){
        cr.j = build_cor_from_events(o$n.events_trt[j,], o$n.events_ctrl[j,])
        oj = comb.pvalue.p23(z1=matrix(o$z1[j, ], nrow=1),  z2 = o$z2[,j], bd.z=bd.z[i,j], w=o$w[,j], selected.dose = n.arms-1, method=multiplicity.method, 
                             cr=array(cr.j, dim = c(nrow(cr.j), ncol(cr.j), 1)))
        comb.z[i, j] = oj$comb.z; 
      }
    } else if (o$method == "Mixture") {
      comb.z[i, ] = o$z.tilde
    }
  }
  
    
  selection = rep(NA, n.arms-1)
  for (j in 1:(n.arms-1)) {
    selection[j] = sum(s == j) / nSim
  }
  
  cum.pow=gsd.power(z = comb.z, bd.z=bd.z)
  
  o = list()
  o$cum.pow = cum.pow
  o$bd.z = bd.z[1:5,]
  o$multiplicity.method = multiplicity.method
  o$method = method
  
  o$selection = selection
  
  #Calculate the generalized power by simulation, defined as the correct selection of the best dose in OS and H0 rejected
  
  #Best dose by design
  doses.m = m[1:(n.arms-1)]
  max.m = max(doses.m)
  
  if (max.m >= m[length(m)]) {
    #There is a best dose in OS.
    best.dose = (1:(n.arms-1))[doses.m == max.m]
    if (sum(s %in% best.dose) > 0) {
      correct.selection = (1:nSim)[s %in% best.dose]
      correct.comb.z = matrix(comb.z[correct.selection, ], ncol = K)
      correct.bd.z = matrix(bd.z[correct.selection, ], ncol = K)
      generalized.pow=gsd.power(z = correct.comb.z, bd.z=correct.bd.z) * length(correct.selection) / nSim
    } else {generalized.pow = 0}
   
    o$best.dose = best.dose 
    o$generalized.pow = generalized.pow
  }
  #o$s=s
  
  return(o)
}

