#' Simulate randomized and controlled phase 2/3 dose optimization trial for survival endpoint without interim dose selection
#'
#' This function simulates a trial that has multiple dose arms at Stage 1 and Stage 2. No adaptation yet.
#' The function returns a data frame including all dose arms in Stage 1 and Stage 2.
#' Features include: (1) Stage 1 and Stage 2 have separate enrollment curves; 
#' (2) Allow enrollment gap between Stage 1 and Stage 2.
#'
#' @param n1 Sample size of each arm at Stage 1.
#' @param n2 Sample size of each arm at Stage 2.
#' @param m Median survival time for each arm (dose 1, dose 2, ..., control). length(m) must be equal to length(n1)
#' @param orr Objective response (binary: 1 = response, 0 = non-response). Must be available for all doses. length(orr) = number of arms. The last one is for control arm.
#' @param rho Correlation between ORR and time to event endpoint
#' @param dose_selection_endpoint Endpoint for dose selection. "ORR" or "not ORR"
#' @param A1 Enrollment period for Stage 1
#' @param Lambda1 Enrollment distribution function (CDF) for stage 1.
#' @param enrollment.hold Holding period in months after DCO1 of Stage 1 prior to enrollment of Stage 2 patients. 0 means seamless enrollment.
#' @param A2 Enrollment period for Stage 2
#' @param Lambda2 Enrollment distribution function (CDF) for stage 2.
#' @param drop Drop out rate per month 
#' 
#' @return A phase 2/3 trial dataframe with variables
#' \describe{
#' \item{enterTime}{Time of randomization to trial}
#' \item{calendarTime}{calendarTime = enterTime + survTime}
#' \item{survTime}{Survival time for analysis calculated from randomization to event/censoring prior to cut off (DCO)}
#' \item{cnsr}{Censoring status prior to cutoff (DCO)}
#' \item{group}{Treatment group}
#' \item{stage}{Stage 1 or 2 for each subject}
#' \item{enrollment.period}{Enrollment period of Stage 1 or Stage 2}
#' \item{enrollment.hold}{Enrollment hold between Stage 1 and Stage 2}
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
#' #Dose selection decision is based on ORR. 
#' p23trial = simu.p23trial(n1 = rep(50, 4), n2 = rep(200, 4), m = c(9,9, 9, 9), 
#' orr = c(0.25, 0.3, 0.4, 0.2), rho = 0.7, dose_selection_endpoint = "ORR",
#' Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A1 = 12,
#' Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A2 = 12,
#' enrollment.hold=4)
#' 
#' #Example (2): 
#' #Dose selection decision is NOT based on ORR.
#' p23trial = simu.p23trial(n1 = rep(50, 4), n2 = rep(200, 4), m = c(9,9, 9, 9), 
#' orr = NULL, rho = NULL, dose_selection_endpoint = "not ORR",
#' Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A1 = 12,
#' Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A2 = 12,
#' enrollment.hold=4)
#' 
#' head(p23trial)
#' 
#' table(p23trial$stage)
#' 
#' @importFrom data.table rbindlist
#' 
#' @export 
#' 
simu.p23trial = function(n1 = c(50, 50, 50, 50), n2 = c(200, 200, 200, 200), m = c(9, 11, 13, 8), 
                         orr = c(0.25, 0.3, 0.4, 0.2), rho = 0.7, dose_selection_endpoint = "ORR",
                         Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A1 = 12,
                         Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A2 = 12,
                         drop = 0,
                         enrollment.hold=4){
   
   ######################
   #Stage 1
   ######################
   
   n.arms = length(n1)
   
   dat1 = list(NULL) #control arm is the last one
      
   #1. simulate Stage 1 survival data
   for (j in 1:n.arms){
       if (dose_selection_endpoint == "ORR"){
         dat1[[j]] = simu.single.arm.hz(n=n1[j], m=m[j], orr = orr[j], rho = rho, Lambda=Lambda1, A=A1, drop=drop, DCO=Inf)[[1]]
       } else {
         dat1[[j]] = simu.single.arm(n=n1[j], m=m[j], Lambda=Lambda1, A=A1, drop=drop, DCO=Inf)[[1]]
       }
       dat1[[j]] <- subset(dat1[[j]], select = -c(calendarCutoff, survTimeCut, cnsrCut))
       dat1[[j]]$group = j; dat1[[j]]$stage = 1
       dat1[[j]]$enrollment.period = A1
       if (j == n.arms){dat1[[j]]$group = 0}
   }
    #combine all treatment group data together 
  dat1 =  data.table::rbindlist(dat1, fill = TRUE)
  
  ######################
  #Stage 2:  
  ######################
  #After decision is made in Stage 1 for dose selection, then enroll Stage 2 patients.
  #If there is enrollment hold
      
  #4. simulate Stage 2 data (control and all dose levels)
  dat2 = list(NULL) #control arm is the last one 
  for (j in 1:n.arms){
    if (dose_selection_endpoint == "ORR"){
      dat2[[j]] = simu.single.arm.hz(n=n2[j], m=m[j], orr = orr[j], rho = rho, Lambda=Lambda2, A=A2, drop=drop, DCO=Inf)[[1]]
    } else {
      dat2[[j]] = simu.single.arm(n=n2[j], m=m[j], Lambda=Lambda2, A=A2, drop=drop, DCO=Inf)[[1]]
    } 
     dat2[[j]]$group = j; 
     if (j == n.arms){dat2[[j]]$group = 0}; 
     dat2[[j]]$stage = 2
     dat2[[j]]$enrollment.period = A2
     
     dat2[[j]] <- subset(dat2[[j]], select = -c(calendarCutoff, survTimeCut, cnsrCut))
     
     #5. Enrollment gap
     dat2[[j]]$enrollment.hold = enrollment.hold
     dat2[[j]] = dplyr::rename(dat2[[j]],
                           enterTime.original ="enterTime",
                           calendarTime.original ="calendarTime")
     dat2[[j]]$enterTime = dat2[[j]]$enterTime.original + enrollment.hold + A1
     dat2[[j]]$calendarTime = dat2[[j]]$enterTime + dat2[[j]]$survTime
     
     dat2[[j]] <- subset(dat2[[j]], select = -c(enterTime.original, calendarTime.original))
  }
  #combine all treatment group data together 
  dat2 =  data.table::rbindlist(dat2, fill = TRUE)
  
  #6. Combine Stage 1 and Stage 2 data
  dat =  data.table::rbindlist(list(dat1, dat2), fill = TRUE)

  return(dat)
}

