#' p value combination approach using closed testing procedure (only up to 3 dose levels)
#'
#' This functions implements p-value combination approach using closed testing procedure for phase 2/3 inferential design with multiple dose levels at Stage 1.
#' Multiple testing adjustment method is Simes method for each family wise hypothesis H_J at stage 1. Then closed testing procedure is applied for weighted z combining stage 1 and stage 2.
#' 
#' 
#' @param z1 A matrix of z values at stage 1. z1 has dimensions of (ntrials, ndoses). ndoses must be 2 or 3.
#' @param p1 A matrix of p values at stage 1. p1 = 1-pnorm(z1). either z1 or p1 must be provided.
#' @param z2 Incremental z test statistic from z1. The two statistics z1 and z2 must be independent. length(z2) = ntrials
#' @param p2 stage 2 p value. p2 = 1-pnorm(z2). Either z2 or p2 must be provided.
#' @param w weight for stage 1 z statistic. If w is a single number, all trials use the same weight in p value combination.
#' @param bd.z rejection boundary in z scale for the combination test
#' @param selected.dose Selected dose. If NULL, the dose will be determined based on max(z1(i, )) for each trial i
#' @param method "simes", "Dunnett".
#' @param cr correlation matrix for Dunnett
#' 
#' @return Returned values include:
#' \describe{
#' \item{selected.dose}{Selected dose}
#' \item{unadj.ps}{unadjusted p value for selected dose s}
#' \item{adj.ps}{adjusted p value for selected dose s in Stage 1 after closed testing procedure using Simes method}
#' \item{comb.p}{combination p value combining stage 1 and stage 2}
#' \item{unadj.zs}{unadjusted z value for selected dose s}
#' \item{adj.zs}{adjusted z value for selected dose s in Stage 1 after closed testing procedure using Simes method}
#' \item{comb.z}{Combination z value combining stage 1 and stage 2}
#' \item{RejPct}{Percentage of rejection, applicable for simulations when z1 has a large number of rows (trials)}
#' \item{Selection}{Probability of selection of each dose, applicable for simulations when z1 has a large number of rows (trials)}
#' }
#' 
#' @examples
#' 
#' #Example (1): Stage 1 has 3 z statistics: z1 = c(0.784, 1.118, 1.941)
#' Select 3rd dose to proceed to Stage 2. At Stage 2, the incremental z test statistic is z2 = 2.1 and the rejection boundary is 1.96.
#' 
#' 
#' comb.pvalue.p23(z1=matrix(c(0.784, 1.118, 1.941), nrow=1),  z2 = 2.1, bd.z=1.96, w=0.2, method="simes")
#' 
#' comb.pvalue.p23(z1=matrix(c(0.784, 1.118, 1.941), nrow=1),  z2 = 2.1, bd.z=1.96, w=0.2, method="dunnett")
#' 
#' comb.pvalue.p23(z1=matrix(c(0.784, 1.118), nrow=1),  z2 = 2.1, bd.z=1.96, w=0.2, method="simes")
#' 
#' comb.pvalue.p23(z1=matrix(c(0.784, 1.118), nrow=1),  z2 = 2.1, bd.z=1.96, w=0.2, method="dunnett")
#' 
#' comb.pvalue.p23(z1=matrix(c(0.784), nrow=1),  z2 = 2.1, bd.z=1.96, w=0.2, method="simes")
#' 
#' comb.pvalue.p23(z1=matrix(c(0.784), nrow=1),  z2 = 2.1, bd.z=1.96, w=0.2, method="dunnett")
#' 
#' #Example (2): check type I error by simulations for a weighted z based on two independent z1 and z2.
#' 
#' nSim = 100000
#' z1 = matrix(NA, nrow=nSim, ncol=3)
#' 
#' for (j in 1:3){z1[ ,j] = rnorm(nSim)}
#' z2 = rnorm(nSim)
#' 
#' o=comb.pvalue.p23(z1=z1,  z2 = z2, bd.z=1.96, w=0.2, method="simes")
#' 
#' #Example (3). w is random between 0 and 1
#' o=comb.pvalue.p23(z1=z1,  z2 = z2, bd.z=1.96, w=runif(nSim), method="simes")
#' 
#' #Example (4): Stage 1: 4 arms; 3 dose levels; each arm 50 patients.
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
#' sel = select.dose.p23 (data=p23trial, DCO1=16)
#' 
#' o=conduct.p23(data=p23trial, DCO1=16, targetEvents = c(300, 380), method = "Independent Incremental")
#' 
#' #Using O'Brien Fleming boundary, the rejection boundaries are: 
#' bd.z = actualBounds(planned.events=c(300, 380), act.events=c(300, 380), sf=gsDesign::sfLDOF, alpha=0.025)$actual.z
#' #2.268527 2.022098
#' 
#' IA=comb.pvalue.p23(z1=matrix(o$z1, nrow=1),  z2 = o$z2[1], bd.z=bd.z[1], w=o$w[1], method="simes")
#' FA=comb.pvalue.p23(z1=matrix(o$z1, nrow=1),  z2 = o$z2[2], bd.z=bd.z[2], w=o$w[2], method="simes")
#' 
#' IA=comb.pvalue.p23(z1=matrix(o$z1, nrow=1),  z2 = o$z2[1], bd.z=bd.z[1], w=o$w[1], method="dunnett")
#' FA=comb.pvalue.p23(z1=matrix(o$z1, nrow=1),  z2 = o$z2[2], bd.z=bd.z[2], w=o$w[2], method="dunnett")
#' 
#' gsd.power(z = cbind(IA$comb.z, FA$comb.z), bd.z=bd.z)
#' 
#' #Force the selected dose = 2, which may be based on ORR selection at Stage 1
#' comb.pvalue.p23(z1=matrix(o$z1, nrow=1),  z2 = o$z2[2], selected.dose = 2, bd.z=bd.z[2], w=o$w[2], method="simes")
#' 
#' comb.pvalue.p23(z1=matrix(o$z1, nrow=1),  z2 = o$z2[2], selected.dose = 2, bd.z=bd.z[2], w=o$w[2], method="dunnett")
#' 
#' @importFrom stats pnorm qnorm
#' @export 
#' 
comb.pvalue.p23 = function(z1, z2, p1=NULL, p2=NULL, bd.z=1.96, w=0.2, selected.dose = NULL, method="simes", cr=NULL){
  n.doses = ncol(z1) #only works for up to 3 dose levels in this program.
  N = nrow(z1) #number of trials
  
  if (!is.null(p1)) {z1 = qnorm(1-p1)}
  if (!is.null(p2)) {z2 = qnorm(1-p2)}
  
  #familywise P(reject H_J | H_J) assuming independent incremental
  
  ##################
  #Stage 1
  ##################
  zs = rep(NA, N) #raw z for Stage 1 selected dose s 
  s = rep(NA, N) #Stage 1 selected dose s
  if (is.null(selected.dose)) {
    for (i in 1:N){
      tmp = sort(z1[i,], index.return = TRUE)
      zs[i] = tmp$x[n.doses]
      s[i] = tmp$ix[n.doses]
    }
  } else {
    s = selected.dose
    for (i in 1:N){
      zs[i] = z1[i, s[i]]
    }  
  }
  
  #cbind(z, zs, ps, s)
  ps = 1-pnorm(zs) #raw p value for selected dose
  
  #Simes and dunnett methods for adjustment
  p = 1 - pnorm(z1)
  p1s = rep(NA, N) # final adjusted p for stage 1 on the selected dose s with Closed testing procedure
  
  if (n.doses > 1) {
    if (method == "simes") {
      simes.p3 = rep(NA, N) #adjusted p by all doses for stage 1 on the selected dose s
      simes.p2 = matrix(NA, nrow=N, ncol=n.doses-1)
    
      for (i in 1:N){
        #adjustment for all p values
        simes.p3[i] = simes(p[i,])
    
        #adjustment for any 2 p values that includes s
        unselected = (1:n.doses)[-s[i]]
        for (j in 1:length(unselected)){simes.p2[i,j] = simes(c(p[i, s[i]], p[i,unselected[j]]))}
        #final adjusted p by CTP
        p1s[i] = max(simes.p3[i], simes.p2[i,], p[i, s[i]])
      }
    } else if (method == "dunnett"){
      dunnett.p3 = rep(NA, N) #adjusted p by all doses for stage 1 on the selected dose s
      dunnett.p2 = matrix(NA, nrow=N, ncol=n.doses-1)
      
      for (i in 1:N){
        #adjustment for all p values
        cr_i = cr[,,i] 
        dunnett.p3[i] = dunnett(p[i,], cr=cr_i)
        
        #adjustment for any 2 p values that includes s
        unselected = (1:n.doses)[-s[i]]
        for (j in 1:length(unselected)){dunnett.p2[i,j] = dunnett(c(p[i, s[i]], p[i,unselected[j]]), cr=cr_i[c(s[i],unselected[j]),c(s[i],unselected[j])])}
        #final adjusted p by CTP
        p1s[i] = max(dunnett.p3[i], dunnett.p2[i,], p[i, s[i]])
      }
    }
  } else {
    p1s = ps
  }
  
  z1s = qnorm(1-p1s) #adjusted z for stage 1 on the selected dose s
  #cbind(s, p, simes.p3, simes.p2, ps, p1s, z1s)
  
  ##################
  #Stage 2, independent incremental z for the selected dose
  ##################

  #Final test statistic
  comb.z = w*z1s + sqrt(1-w^2)*z2
  comb.p = 1-pnorm(comb.z)
  
  o = list()
  o$selected.dose = s
  o$unadj.ps = ps
  o$adj.ps = p1s
  o$comb.p = comb.p
  
  o$unadj.zs = zs
  o$adj.zs = z1s
  o$comb.z = comb.z
  o$method = method
  
  selection = rep(NA, n.doses)
  if (N > 1){ 
    RejPct = sum(comb.z > bd.z) / N
    o$RejPct = RejPct
    for (j in 1:n.doses){
      selection[j] = sum(s == j) / N
    }
    o$selection = selection
  }
  
  return(o)
}

