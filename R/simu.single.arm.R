#' Simulate a single arm survival data with non-uniform enrollment distribution 
#' and drop-off
#'
#' This functions simulates trials that have multiple dose arms at Stage 1 and 
#' the best dose is always selected after Stage 1 and perform multiple analyses 
#' at Stage 2. The function returns multiple datasets including at end of Stage 
#' 1 (all dose arms) and at each of analysis at Stage 2
#'
#' @param n Sample size
#' @param m Median survival time
#' @param Lambda Enrollment distribution function (CDF) 
#' @param A Enrollment period
#' @param drop Dropout rate per year
#' @param DCO Data cutoff date for Stage 1
#' @param targetEvents Planned target number of events. Either targetEvents or 
#' DCO must be provided. 
#' 
#' @return Datasets including data1: Stage 1 data with multiple dose arms; 
#' data2: Stage 2 data of the selected dose arm and control arm with data2((j)) 
#' for jth analysis. Each dataset includes variables
#' \describe{
#' \item{treatment}{treatment group with values of "control", "dose 1", 
#' "dose 2", ...}
#' \item{enterTime}{Time of randomization in calendar time}
#' \item{calendarTime}{the time when event/censoring occurred in calendar time}
#' \item{survTime}{Survival time for analysis, = calendarTime - enterTime}
#' \item{cnsr}{censor status (0=event; 1=censor) before administrative 
#' censoring due to data cut}
#' \item{calendarCutOff}{Data CutOff Time (DCO);}
#' \item{survTimeCut}{Survival time after cut}
#' \item{cnsrCut}{Censor status after cut}
#' \item{analysis}{Sequence of analysis}
#' }
#' 
#' @examples
#'  
#' Lambda = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}  
#' dat=simu.single.arm(n = 100, m = 10, Lambda = Lambda, A = 12, 
#' drop=0, DCO = 16, targetEvents = NULL)
#' 
#' dat=simu.single.arm(n = 100, m = 10, Lambda = Lambda, A = 12, 
#' drop=0, DCO = c(16, 28), targetEvents = NULL)
#' 
#' dat=simu.single.arm(n = 100, m = 10, Lambda = Lambda, A = 12, drop=0, 
#' DCO = NULL, targetEvents = 40)
#' 
#' dat=simu.single.arm(n = 100, m = 10, Lambda = Lambda, A = 12, drop=0, 
#' DCO = NULL, targetEvents = c(40, 60))
#' 
#' @importFrom stats rexp runif
#' @export 
#' 
simu.single.arm = function(n = 100, m = 10, 
                           Lambda = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)},
                           A = 12, drop=0, DCO = 16, targetEvents = NULL){
  
  #(1) Survival data
  t = genSurv(dist = "exponential", n = n, lam=log(2)/m)

  #(2) Drop Off data
  ############################
  if (drop > 0) {W = rexp(n, rate=-log(1-drop)/12)} else {W = rep(Inf, n)}

  ############################
  #Censor data from Drop-off
  ############################
  Y = apply(cbind(t, W), 1, min)
  
  event = as.numeric(t < Inf)
  event[W < t] = 0
  
  #(3) Enrollment data
  #Trick the program to generate a single arm data
  nEachMonth = p23::f.nEachMonth(N=2*n, A=A, w=NULL, r=1, Lambda=Lambda)$n0
  
  ############################
  #EnterTime, CalendarTime
  ############################
  
  enterTime = rep(NA, n)
  enterTime[1:nEachMonth[1]] = runif(nEachMonth[1], min=0, max=1)
  if (A > 1){for (j in 2:A){
    LL = sum(nEachMonth[1:(j-1)])+1
    UU = sum(nEachMonth[1:j])
    enterTime[LL:UU] = runif(nEachMonth[j]) + j - 1
  }}
  enterTime = enterTime * A / floor(A)
  survTime = as.numeric(Y)
  
  #trick infinity
  survTime[survTime > 1e6] = 1e6
  calendarTime = as.numeric(enterTime) + as.numeric(survTime)
  cnsr = 1-event
  
  dati = data.frame(cbind(enterTime, calendarTime, survTime, cnsr))
  
  ############################
  #(4) Cut data
  ############################
  L = ifelse(!is.null(DCO), length(DCO), length(targetEvents))
  
  dat.cut = NULL
  for (ii in 1:L){
    dat.cut[[ii]] = f.dataCut(data=dati, targetEvents=targetEvents[ii], DCO = DCO[ii])
    dat.cut[[ii]]$analysis = ii
  }
  return(dat.cut)
}


