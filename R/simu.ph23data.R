#' Simulate randomized and controlled phase 2/3 dose optimization trial for survival endpoint
#'
#' This functions simulates trials that have multiple dose arms at Stage 1 and the best dose is always selected after Stage 1 and perform multiple analyses at Stage 2
#' The function returns multiple datasets including at end of Stage 1 (all dose arms) and at each of analysis at Stage 2.
#' Features include: (1) Stage 1 and Stage 2 have separate enrollment curves; (2) Allow enrollment gap between Stage 1 and Stage 2; 
#' (3) Dose selection at stage 1 based on best performance in z statistic. (4) Stage 2 analyses are based on prespecified target events. (5) The function
#' also returns the weights for combination p value approach in the next step and the weights are sqrt of information fractions.
#'
#' @param nSim Number of simulated trials
#' @param n1 Sample size of (dose arm, control arm) at Stage 1. length(n1) must be 2.
#' @param n2 Sample size of the (selected dose arm, control arm). length(n2) must be 2.
#' @param m Median survival time for each arm (dose 1, dose 2, ..., control). length(m) must be equal to length(n1)
#' @param A1 Enrollment period for Stage 1
#' @param Lambda1 Enrollment distribution function (CDF) for stage 1.
#' @param DCO1 Data cutoff date for Stage 1
#' @param drop Drop off rate per month
#' @param enrollment.hold Holding period in months after DCO1 of Stage 1 prior to enrollment of Stage 2 patients. 0 means seamless enrollment.
#' @param A2 Enrollment period for Stage 2
#' @param Lambda2 Enrollment distribution function (CDF) for stage 2.
#' @param targetEvents2 Planned target number of events for Stage 2. Either targetEvents2 must be provided. 
#' @param method "Independent Incremental", "Disjoint Subjects". Currently, only "Independent Incremental" method is implemented.
#' The multiplicity adjustment is performed according to the extended followup data up to 1st analysis. This option is applicable to "Disjoined Subjects" method. Currently, not implemented yet.
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
#' #Using O'Brien Fleming boundary, the rejection boundaries are: 
#' bd.z = actualBounds(planned.events=c(300, 380), act.events=c(300, 380), sf=gsDesign::sfLDOF, alpha=0.025)$actual.z
#' #2.268527 2.022098
#' 
#' o = simu.ph23data(nSim=100, n1 = c(50, 50, 50, 50), n2 = c(200, 200), m = c(9,9, 9, 9), 
#' Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A1 = 12,
#' DCO1 = 16, Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A2 = 12,
#' enrollment.hold=4, targetEvents2 = c(300, 380))
#' 
#' #Perform analysis at IA and FA using combination p method, with simes approach for multiplcity control at Stage 1 and closed testing procedure for controlling FWER
#' IA = ph23.comb.p(z1=o$z1,  z2 = o$z2[,1], bd.z=bd.z[1], w=o$w[,1])
#' FA = ph23.comb.p(z1=o$z1,  z2 = o$z2[,2], bd.z=bd.z[2], w=o$w[,2])
#' 
#' #power calculation using the standard group sequential boundaries
#' gsd.power(z = cbind(IA$comb.z, FA$comb.z), bd.z=bd.z)
#' 
#' @importFrom dplyr rename
#' 
#' @export 
#' 
simu.ph23data = function(nSim=1000, n1 = c(50, 50, 50, 50), n2 = c(200, 200), m = c(9, 11, 13, 8), 
                          Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A1 = 12,
                          DCO1 = 16, Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A2 = 12,
                         drop = 0,
                          enrollment.hold=4, targetEvents2 = c(300, 380), method = "Independent Incremental"){
   
  if (method == "Independent Incremental"){
   ######################
   #Stage 1
   ######################
   
   n.arms = length(n1)
   z1 = e1 = matrix(NA, nrow=nSim, ncol=n.arms-1)
   
   #e1: total number of events for by DCO1 for each dose arm + control
   s = rep(NA, nSim) #selected dose level
   z.c = matrix(NA, nrow=nSim, ncol=length(targetEvents2)) #z value at each of Stage 2 analyses using combined data
   z2  = matrix(NA, nrow=nSim, ncol=length(targetEvents2)) # incremental z value from stage 2 compared to Stage 1
   w = matrix(NA, nrow=nSim, ncol=length(targetEvents2)) #weights as sqrt(information fraction) = sqrt(events at dose selection / events at each stage 2 analysis) at each analysis.

   for (i in 1:nSim){
     dat1 = list(NULL) #control arm is the last one
     
     #1. simulate Stage 1 survival data
     for (j in 1:n.arms){
       dat1[[j]] = simu.single.arm(n=n1[j], m=m[j], Lambda=Lambda1, A=A1, drop=drop, DCO=DCO1)[[1]]
       dat1[[j]]$group = j
       if (j == n.arms){dat1[[j]]$group = 0}
     }
     
     #2. Calculate z statistics and e
     for (j in 1:(n.arms-1)){
       datj = rbind(dat1[[j]], dat1[[n.arms]])
       z1[i, j] = logrank.one.sided(time=datj$survTimeCut, cnsr=datj$cnsrCut, group=datj$group)$z
       e1[i, j] = sum(1-datj$cnsrCut)
     }
     
     #3. Dose selection
     tmp = sort(z1[i,], index.return = TRUE)
     s[i] = tmp$ix[n.arms-1]
   
     #cbind(z1, s)
   
     ######################
     #Stage 2
     ######################
     #After decision is made in Stage 1 for dose selection, then enroll Stage 2 patients.
     #If there is enrollment hold
      
     #4. simulate Stage 2 data (control and selected dose s)
     dat20 = simu.single.arm(n=n2[2], m=m[n.arms], Lambda=Lambda2, A=A2, drop=drop, DCO=Inf)[[1]]
     dat2s = simu.single.arm(n=n2[1], m=m[s[i]], Lambda=Lambda2, A=A2, drop=drop, DCO=Inf)[[1]]
     dat20$group=0; dat2s$group = s[i]

     #5. Enrollment gap
     dat2a = rbind(dat20, dat2s)
     dat2a$enrollment.hold = enrollment.hold
   
     #Rename variables
     dat2b <- subset(dat2a, select = -c(calendarCutoff, survTimeCut, cnsrCut))
     dat2c = dplyr::rename(dat2b,
                 enterTime.original ="enterTime",
                 calendarTime.original ="calendarTime")
     dat2c$enterTime = dat2c$enterTime.original + enrollment.hold + A1
     dat2c$calendarTime = dat2c$enterTime + dat2c$survTime
   
     #6. Combine Stage 1 and Stage 2 data
     dat.s = rbind(dat1[[s[i]]], dat1[[n.arms]])
     dat.s = subset(dat.s, select = -c(calendarCutoff, survTimeCut, cnsrCut))
     dat.s$stage = 1 #flag this is stage 2 data
     
     dat2d = as.data.frame(cbind(dat2c$enterTime, dat2c$calendarTime, dat2c$survTime, dat2c$cnsr, dat2c$group))
     dat2d = dplyr::rename(dat2d,
                         enterTime ="V1",
                         calendarTime ="V2",
                         survTime="V3",
                         cnsr="V4",
                         group="V5")
     dat2d$stage = 2 #flag this is stage 2 data
   
     dat.comb0 = rbind(dat.s, dat2d)
     
     #7. Data cut according to K analyses in stage 2
   
     K = length(targetEvents2) #K analyses in Stage 2
   
     for (ii in 1:K){
       dat.comb = f.dataCut(data=dat.comb0, targetEvents=targetEvents2[ii])
       
       #Calculate z statistics for logrank test combining Stage 1 and stage 2
       z.c[i, ii] = logrank.one.sided(time=dat.comb$survTimeCut, cnsr=dat.comb$cnsrCut, group=dat.comb$group)$z
       frac = e1[i, s[i]]/targetEvents2[ii]     
       z2[i, ii] = (z.c[i, ii] - sqrt(frac)*z1[i, s[i]])/sqrt(1-frac)
       w[i, ii] = sqrt(frac)
     }
   }
   #################################
   #################################
  } else if (method == "Disjoint Subjects"){
    
    #Extended followup; multiplicity adjustment is based on Stage 1 data up to 1st analysis at Stage 2
    #In this approach, the data cutoff for multiplicity adjustment is IA1 (1st analysis), but the dose selection
    #is at IAd. There could be misalignment between max(z1 at IAd) and max(z1 at IA1). However, the decision is based on
    #max(z1 at IAd).
    
    ######################
    #Stage 1
    ######################
    
    n.arms = length(n1)
    z1 = e1 = matrix(NA, nrow=nSim, ncol=n.arms-1)
    
    #e1: total number of events for by DCO1 for each dose arm + control
    s = rep(NA, nSim) #selected dose level
    z.c = matrix(NA, nrow=nSim, ncol=length(targetEvents2)) #z value at each of Stage 2 analyses using combined data
    z2  = matrix(NA, nrow=nSim, ncol=length(targetEvents2)) # incremental z value from stage 2 compared to Stage 1
    w = matrix(NA, nrow=nSim, ncol=length(targetEvents2)) #weights as sqrt(information fraction) = sqrt(events at dose selection / events at each stage 2 analysis) at each analysis.
    
    for (i in 1:nSim){
      dat1 = list(NULL) #control arm is the last one
      
      #1. simulate Stage 1 survival data
      for (j in 1:n.arms){
        dat1[[j]] = simu.single.arm(n=n1[j], m=m[j], Lambda=Lambda1, A=A1, drop=drop, DCO=DCO1)[[1]]
        dat1[[j]]$group = j
        if (j == n.arms){dat1[[j]]$group = 0}
      }
      
      #2. Calculate z statistics and e
      for (j in 1:(n.arms-1)){
        datj = rbind(dat1[[j]], dat1[[n.arms]])
        z1[i, j] = logrank.one.sided(time=datj$survTimeCut, cnsr=datj$cnsrCut, group=datj$group)$z
        e1[i, j] = sum(1-datj$cnsrCut)
      }
      
      #3. Dose selection
      tmp = sort(z1[i,], index.return = TRUE)
      s[i] = tmp$ix[n.arms-1]
      
      #cbind(z1, s)
      
      ######################
      #Stage 2
      ######################
      #After decision is made in Stage 1 for dose selection, then enroll Stage 2 patients.
      #If there is enrollment hold
      
      #4. simulate Stage 2 data (control and selected dose s)
      dat20 = simu.single.arm(n=n2[2], m=m[n.arms], Lambda=Lambda2, A=A2, drop=drop, DCO=Inf)[[1]]
      dat2s = simu.single.arm(n=n2[1], m=m[s[i]], Lambda=Lambda2, A=A2, drop=drop, DCO=Inf)[[1]]
      dat20$group=0; dat2s$group = s[i]
      
      #5. Enrollment gap
      dat2a = rbind(dat20, dat2s)
      dat2a$enrollment.hold = enrollment.hold
      
      #Rename variables
      dat2b <- subset(dat2a, select = -c(calendarCutoff, survTimeCut, cnsrCut))
      dat2c = dplyr::rename(dat2b,
                            enterTime.original ="enterTime",
                            calendarTime.original ="calendarTime")
      dat2c$enterTime = dat2c$enterTime.original + enrollment.hold + A1
      dat2c$calendarTime = dat2c$enterTime + dat2c$survTime
      
      #6. Combine Stage 1 and Stage 2 data
      dat.s = rbind(dat1[[s[i]]], dat1[[n.arms]])
      dat.s = subset(dat.s, select = -c(calendarCutoff, survTimeCut, cnsrCut))
      dat.s$stage = 1 #flag this is stage 2 data
      
      dat2d = as.data.frame(cbind(dat2c$enterTime, dat2c$calendarTime, dat2c$survTime, dat2c$cnsr, dat2c$group))
      dat2d = dplyr::rename(dat2d,
                            enterTime ="V1",
                            calendarTime ="V2",
                            survTime="V3",
                            cnsr="V4",
                            group="V5")
      dat2d$stage = 2 #flag this is stage 2 data
      
      dat.comb0 = rbind(dat.s, dat2d)
      
      #7. Data cut according to K analyses in stage 2
      
      K = length(targetEvents2) #K analyses in Stage 2
      
      for (ii in 1:K){
        dat.comb = f.dataCut(data=dat.comb0, targetEvents=targetEvents2[ii])
        
        #Calculate z statistics for logrank test combining Stage 1 and stage 2
        z.c[i, ii] = logrank.one.sided(time=dat.comb$survTimeCut, cnsr=dat.comb$cnsrCut, group=dat.comb$group)$z
        frac = e1[i, s[i]]/targetEvents2[ii]     
        z2[i, ii] = (z.c[i, ii] - sqrt(frac)*z1[i, s[i]])/sqrt(1-frac)
        w[i, ii] = sqrt(frac)
      }
    }
    
  }
  
   o=list()
   o$z1 = z1
   o$z2 = z2
   o$z.c = z.c
   o$w = w
   o$selected.dose = s
   
   return(o)
 }
