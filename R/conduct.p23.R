#' Conduct a phase 2/3 dose optimization trial for survival endpoint
#'
#' This function calculates the stage 2 z statistics according to preplanned target number of events at each analysis.
#' The function returns (1) Z statistics at Stage 2 based on prespecified target events for pooled stage 1 (selected dose) and stage 2 patients. (2) The function
#' also returns the weights for combination p value approach in the next step and the weights are sqrt of information fractions. 
#'
#' @param data Dataset produced by simu.p23trial() function
#' @param DCO1 Data cutoff for Stage 1
#' @param targetEvents2 Planned target number of events for Stage 2. Either targetEvents2 must be provided.
#' @param dose_selection_endpoint Endpoint for dose selection at Stage 1. "ORR" or "not ORR"
#' @param method Options include "Independent Incremental": z1 at dose selection and z2 is from dose selection to kth analysis at stage 2; 
#' "Disjoint Subjects": z1 is at kth analysis for stage 1 subjects; z2 is at the kth analysis for stage 2 subjects. z1 will be adjusted by multiplicity and closed testing procedure at each analysis.
#' "Mixture": Only consider disjoint subjects at first analysis in stage 2. Starting from the 2nd analysis, consider independent incremental methods. Only z1 at 1st analysis will be adjusted by multiplicity and closed testing procedure.
#' @param method "simes", "Dunnett". 
#' @param ORRdiff The tolerated difference between ORR of highest dose and lower doses, with default value 0. 
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
#' \item{v}{Weight for 2nd and further analysis combining z1.tilde and z2[k] to calculate z.tilde[k]}
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
#' conduct.p23(data=p23trial, DCO1=16, dose_selection_endpoint = "not ORR",
#' targetEvents2 = c(300, 380), method = "Independent Incremental")
#' 
#' #Disjoint Subjects; dose not selected by ORR
#' conduct.p23(data=p23trial, DCO1=16, dose_selection_endpoint = "not ORR",
#' targetEvents2 = c(300, 380), method = "Disjoint Subjects")
#' 
#' #Mixture; dose not selected by ORR
#' conduct.p23(data=p23trial, DCO1=16, dose_selection_endpoint = "not ORR",
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
#' conduct.p23(data=p23trial, DCO1=16, dose_selection_endpoint = "ORR",
#' targetEvents2 = c(300, 380), method = "Independent Incremental")
#' 
#' #Disjoint Subjects; dose not selected by ORR
#' d = conduct.p23(data=p23trial, DCO1=16, dose_selection_endpoint = "ORR",
#' targetEvents2 = c(300, 380), method = "Disjoint Subjects")
#' 
#' #Mixture; dose not selected by ORR
#' conduct.p23(data=p23trial, DCO1=16, dose_selection_endpoint = "ORR",
#' targetEvents2 = c(300, 380), method = "Mixture")
#' 
#' #Rejection boundary
#' bd.z = gsDesign::gsDesign(k=2,alpha=0.025,timing=targetEvents2/targetEvents2[K],sfu=gsDesign::sfLDOF, test.type=1)$upper$bound
#' 
#' #Combination p value
#' comb.pvalue.p23(z1=matrix(d$z1[1, ], nrow=1),  z2 = d$z2[,1], bd.z=bd.z[1], w=d$w[,1], selected.dose = sel$s, method="simes")
#' @import dplyr 
#' @export 
#' 
conduct.p23 = function(data=NULL, DCO1=16, targetEvents2 = c(300, 380), dose_selection_endpoint = "ORR",
                       method = "Independent Incremental", multiplicity.method="simes", e1=NULL, ORRdiff=0){

  #1. Dose selection  
  sel = select.dose.p23 (data=data, DCO1=DCO1, dose_selection_endpoint = dose_selection_endpoint, ORRdiff=ORRdiff)
  s=sel$s

  

  if(method=="Mixture"){
    # Determine IA timing YC============
    n1 = data %>% filter(group==0, stage==1) %>% nrow()
    DCO_IAstage1<- data %>%
      arrange(calendarTime) %>%               # Sort by time
      dplyr::filter(cnsr == 0, stage==1, group%in%c(0, s)) %>%          # Keep only rows with events in selected arm and control arms, `group` is 0,1,2,... but `s` starts with 1
      slice(round(n1*2*0.637)) %>%                   # Select the DCO event for IA
      pull(calendarTime)                      # Extract the time of the DCO event for IA
    
    targetEvents.IA.candidate = data %>%
      arrange(calendarTime) %>%               # Sort by time
      dplyr::filter(cnsr==0, group%in%c(0, s), calendarTime<=DCO_IAstage1) %>% 
      nrow() 
    
    DCO_FA <- data %>%
      arrange(calendarTime) %>%               # Sort by time
      dplyr::filter(cnsr == 0, group%in%c(0, s)) %>%          # Keep only rows with events in selected arm and control arms
      slice(targetEvents2[length(targetEvents2)]) %>%                   # Select the FA event
      pull(calendarTime)                      # Extract the time of the DCO event for FA
    
    if(DCO_IAstage1 >= DCO_FA){
      warning("Warning: The timing of IA exceeds timing of FA.")
      o = list()
      o$s = s
      o$dose.selection.endpoint=dose_selection_endpoint
      o$method = "NA"
      dat23 = data[data$group == 0 | data$group == s, ]
      dat23k = f.dataCut(data=dat23, targetEvents=targetEvents2[length(targetEvents2)])
      o$z = logrank.one.sided(time=dat23k$survTimeCut, cnsr=dat23k$cnsrCut, group=dat23k$group)$z
      o$actualEvents = rep(targetEvents2[length(targetEvents2)],2)
      return(o)
    }
    targetEvents2[1] = max(targetEvents2[1], targetEvents.IA.candidate)
  }
  
  
  
  #2. Assemble the trial data combining stage 1 and stage 2 for selected dose + control
  dat23 = data[data$group == 0 | data$group == s, ]
  
  K = length(targetEvents2) #K analyses in Stage 2
  n.arms = length(unique(data$group))
  
  #3. Determine z1 and z2
  if (method == "Independent Incremental"){
    #1st component in weighted z as of dose selection IAd.
    z1 = sel$z1; e1=sel$e1; 
    
    #2nd component in weighted z as of incremental from IAd to IA
    z.c = rep(NA, K) #z value at each of Stage 2 analyses using combined data
    z2  = rep(NA, K) # incremental z value from stage 2 compared to Stage 1 for "independent incremental" approach.
    w = rep(NA, K) #weights as sqrt(information fraction) = sqrt(events at dose selection / events at each stage 2 analysis) at each analysis.
    
    for (k in 1:K){
      #4. data cut for each analysis
      dat23k = f.dataCut(data=dat23, targetEvents=targetEvents2[k])
      
      #5. Calculate z statistics for logrank test combining Stage 1 and stage 2
      z.c[k] = logrank.one.sided(time=dat23k$survTimeCut, cnsr=dat23k$cnsrCut, group=dat23k$group)$z
      
      #2nd component z statistic and weight for combining z1 and z2 in next step
      frac.k = e1[s]/targetEvents2[k]  
      z2[k] = (z.c[k] - sqrt(frac.k)*z1[s])/sqrt(1-frac.k)
      
      w[k] = sqrt(frac.k)
    }
    
  } else if (method == "Disjoint Subjects") {
    
    # identify Events from stage 1 patients at IA and FA with K==2 YC====================
    actualEvents_S1_IA <- dat23 %>%
      arrange(calendarTime) %>%               # Sort by time
      dplyr::filter(cnsr == 0) %>%     
      mutate(nevents=cumsum(cnsr==0)) %>% 
      filter(nevents<= targetEvents2[1], stage==1) %>% nrow()    
    
    dat.FA <- dat23 %>%
      arrange(calendarTime) %>%               # Sort by time
      dplyr::filter(cnsr == 0) %>%
      mutate(nevents=cumsum(cnsr==0)) %>% 
      filter(nevents<= targetEvents2[length(targetEvents2)])
    
    actualEvents_S1_FA <- dat.FA %>%
      filter(stage==1) %>% nrow()    
    
    actualTime_FA <- max(dat.FA$calendarTime)
    
    #1st component in weighted z as of Stage 1 subjects at each analysis
    z11 = rep(NA, K)    
    
    #For each analysis k, z1 has unselected doses at IAd for multiplicity adjustment and z for stage 1 subjects
    z1 = matrix(NA, nrow=K, ncol=n.arms-1)
    
    #2nd component in weighted z as of Stage 2 subjects at each analysis
    z2 = z21 = w = rep(NA, K)
    
    # number of events at each analysis YC ==============================
    n.events = matrix(NA, nrow=K, ncol=2)

    #multiplicity adjustment needs z from IAd for the unselected doses
    z1.IAd = sel$z1[-s]
    # # modified by YC ==========================================================
    # z1s = rep(NA, n.arms-1)
    # dat1cut.IAD = f.dataCut(data=data[data$stage == 1 & data$group!=0,], DCO=DCO1)
    # 
    # dat1cutIA = f.dataCut(data=data[data$group%in%c(0,s)], targetEvents = targetEvents2[1]) %>% 
    #   dplyr::filter(stage==1, group==0)
    # 
    # dat1cut = rbind(dat1cut.IAD, dat1cutIA)
    # 
    # for (j in 1:(n.arms-1)){
    #   datj = dat1cut[(dat1cut$group == 0 | dat1cut$group == j), ]
    #   z1s[j] = logrank.one.sided(time=datj$survTimeCut, cnsr=datj$cnsrCut, group=datj$group)$z
    # }
    # z1.IAd = z1s[-s]
    # # END of modification =====================================================
    
    
    for (k in 1:K){
      #4. data cut for kth analysis
      dat23k = f.dataCut(data=dat23, targetEvents=targetEvents2[k])
      
      n.events[k,] = t(table(dat23k$group, dat23k$cnsrCut))[1,]
      
      #stage 1 subjects at analysis k
      dat1k = dat23k[dat23k$stage == 1, ] 
      
      #z11 z statistic for stage 1 subjects at analysis k
      z11[k] = logrank.one.sided(time=dat1k$survTimeCut, cnsr=dat1k$cnsrCut, group=dat1k$group)$z
      
      #1st component z statistic
      z1[k,] = c(z1.IAd, z11[k]) # ??? last element is selected arm? =================
      
      #stage 2 subjects at analysis k
      dat2k = dat23k[dat23k$stage == 2, ] 
      
      #z21 z statistic for stage 1 subjects at analysis k
      z21[k] = logrank.one.sided(time=dat2k$survTimeCut, cnsr=dat2k$cnsrCut, group=dat2k$group)$z
      
      #Number of events for stage 1 subjects
      # e1k = sum(1-dat1k$cnsrCut)
      frac.k = e1[k,s]/targetEvents2[k] # pre-specified weight YC ===================
      
      #weight w for combining z1 and z2 in next step
      w[k] = sqrt(frac.k)
    }
    z2 = z21
  } else if (method == "Mixture") {
    #Only adjust multiplicity for Stage 1 patients at 1st analysis, 
    #then construct following z as independent incrementals
    
    #Multiplicity adjustment needs z from IAd for the unselected doses
    #multiplicity adjustment needs z from IAd for the unselected doses
    # ??? this also needs to be changed?? YC ==========================================
    # z1.IAd = sel$z1[-s] 
    # modified by YC ==========================================================
    z1s = rep(NA, n.arms-1)
    dat1cut.IAD = f.dataCut(data=data[data$stage == 1 & data$group!=0,], DCO=DCO1)
    
    dat1cutIA = f.dataCut(data=data[data$group%in%c(0,s)], targetEvents = targetEvents2[1]) %>% 
      dplyr::filter(stage==1, group==0)
    
    dat1cut = rbind(dat1cut.IAD, dat1cutIA)
    
    for (j in 1:(n.arms-1)){
      datj = dat1cut[(dat1cut$group == 0 | dat1cut$group == j), ]
      z1s[j] = logrank.one.sided(time=datj$survTimeCut, cnsr=datj$cnsrCut, group=datj$group)$z
    }
    z1.IAd = z1s[-s]
    # END of modification =====================================================
    
    ##############
    #1st analysis
    ##############
    
    #data cut
    dat23.1 = f.dataCut(data=dat23, targetEvents=targetEvents2[1])
    
    #stage 1 subjects
    dat23.11 = dat23.1[dat23.1$stage == 1, ]
    
    #z statistic for stage 1 subjects
    z11 = logrank.one.sided(time=dat23.11$survTimeCut, cnsr=dat23.11$cnsrCut, group=dat23.11$group)$z
    
    #stage 2 subjects
    dat23.21 = dat23.1[dat23.1$stage == 2, ]
    
    #z statistic for stage 2 subjects
    z21 = logrank.one.sided(time=dat23.21$survTimeCut, cnsr=dat23.21$cnsrCut, group=dat23.21$group)$z
    
    #Number of events for stage 1 subjects
    # e11 = sum(1-dat23.11$cnsrCut)
    # frac1 = e11/targetEvents2[1]
    frac1 = e1[1,s]/targetEvents2[1] # pre-specified weight YC ===================
    
    #weight w for combining z1 and z2 in next step
    omega1 = sqrt(frac1)
    
    z.tilde = rep(NA, K) #test statistics
    
    #adjusted z1 by CTP and simes method
    z.tilde[1] = comb.pvalue.p23(z1=matrix(c(z1.IAd, z11), nrow=1),  z2 = z21, selected.dose = n.arms-1, bd.z=Inf, w=omega1, method=multiplicity.method)$comb.z

    ##################
    #>= 2 analysis
    #Independent incrementals
    ##################
    
    if (K > 1) {
      z.c = rep(NA, K) #z value at each of Stage 2 analyses using combined data
      z.c[1] = logrank.one.sided(time=dat23.1$survTimeCut, cnsr=dat23.1$cnsrCut, group=dat23.1$group)$z
      z2 = rep(NA, K-1) #independent incrementals
      v =  rep(NA, K-1) #weights for independent incrementals
      
      for (k in 2:K){
        #4. data cut for kth analysis
        dat23k = f.dataCut(data=dat23, targetEvents=targetEvents2[k])
        
        #5. Calculate z statistics for logrank test combining Stage 1 and stage 2
        z.c[k] = logrank.one.sided(time=dat23k$survTimeCut, cnsr=dat23k$cnsrCut, group=dat23k$group)$z
        
        #2nd component z statistic and weight for combining z1 and z2 in next step
        frac.k = targetEvents2[1]/targetEvents2[k]
        z2[k-1] = (z.c[k] - sqrt(frac.k)*z.c[1])/sqrt(1-frac.k)
        
        v[k-1] = sqrt(frac.k)
        
        z.tilde[k] = v[k-1] * z.tilde[1] + sqrt(1-v[k-1]^2) * z2[k-1]
      }
    }
    
  }
  
  o = list()
  o$s = s
  
  if (method == "Independent Incremental"){
    o$z1 = matrix(z1, nrow=1)
    o$z2 = matrix(z2, nrow=1)
    o$z.c = matrix(z.c, nrow=1)
    o$w = matrix(w, nrow=1)
  } else if (method == "Disjoint Subjects") {
    o$actualEventsS1 = c(actualEvents_S1_IA, actualEvents_S1_FA)
    o$z1 = z1; o$z2 = matrix(z2, nrow=1); o$w = matrix(w, nrow=1)
    o$actualTime_FA = actualTime_FA
    o$n.events_trt = matrix(rep(sel$n.events_IAD[-1], K), byrow=TRUE, nrow = K)
    o$n.events_trt[,s] = n.events[,2]
    o$n.events_ctrl = matrix(rep(sel$n.events_IAD[1], (n.arms-1)*K), byrow=TRUE, nrow = K)
    o$n.events_ctrl[,s] = n.events[,1]
  } else if (method == "Mixture") {
    o$z1.unselected = z1.IAd
    o$z11 = z11; o$z21 = z21; o$omega1 = omega1
    o$z.tilde=z.tilde; o$z.c = z.c; o$z2 = z2; o$v = v
  }
  if (dose_selection_endpoint == "ORR") {
    o$orr.diff = matrix(sel$orr.diff, nrow=1)
  }
  o$dose.selection.endpoint=dose_selection_endpoint
  o$method = method
  
  o$actualEvents = targetEvents2
  
  return(o)
}

  

