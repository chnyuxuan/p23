#' Select dose at Stage 1 in a phase 2/3 dose optimization trial for survival endpoint
#'
#' This function calculates the z statistics (z1) of all arms in stage 1, and select the dose based on maximum z in stage 1.
#'
#' @param data Dataset produced by simu.p23trial() function.
#' @param DCO1 Data cutoff date for stage 1 dose selection
#' @param dose_selection_endpoint  Dose selection end point: "ORR" or "not ORR"
#' @param ORRdiff The tolerated difference between ORR of highest dose and lower doses, with default value 0. 
#' 
#' @return 
#' \describe{
#' \item{z1}{Stage 1 z values for each dose}
#' \item{e1}{Stage 1 number of events for each dose}
#' \item{s}{Selected dose}
#' }
#' 
#' @examples
#'
#' #Example (1): Stage 1: 4 arms; 3 dose levels; each arm 50 patients.
#' #Stage 2: additional 200 patients per arm will be enrolled at stage 2
#' #medians for the 4 arms: 9, 11, 13 and control = 8 months
#' #Enrollment: 12 months uniform in stage 1; 12 months uniform in stage 2
#' #Holding period: 4 months between stage 1 and 2
#' #Dose selection will be based on data cut at 16 months
#' #Stage 2 has 2 planned analyses at 300 and 380 events respectively.
#'
#' 
#' p23trial = simu.p23trial(n1 = rep(50, 4), n2 = rep(200, 4), m = c(9,9, 9, 9), 
#' orr = NULL, rho = NULL, dose_selection_endpoint = "not ORR",
#' Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A1 = 12,
#' Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A2 = 12,
#' enrollment.hold=4)
#' 
#' select.dose.p23 (data=p23trial, DCO1=16, dose_selection_endpoint = "not ORR")
#' 
#' #' #Example (2): Stage 1: 4 arms; 3 dose levels; each arm 50 patients.
#' #Stage 2: additional 200 patients per arm will be enrolled at stage 2
#' #medians for the 4 arms: 9, 11, 13 and control = 8 months
#' #Enrollment: 12 months uniform in stage 1; 12 months uniform in stage 2
#' #Holding period: 4 months between stage 1 and 2
#' #Dose selection will be based on ORR with data cut at 16 months
#' #Stage 2 has 2 planned analyses at 300 and 380 events respectively.
#'
#' 
#' p23trial = simu.p23trial(n1 = rep(50, 4), n2 = rep(200, 4), m = c(9,9, 9, 9), 
#' orr = c(0.25, 0.3, 0.4, 0.2), rho = 0.7, dose_selection_endpoint = "ORR",
#' Lambda1 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A1 = 12,
#' Lambda2 = function(t){(t/12)*as.numeric(t<= 12) + as.numeric(t > 12)}, A2 = 12,
#' enrollment.hold=4)
#' 
#' 
#' select.dose.p23 (data=p23trial, DCO1=16, dose_selection_endpoint = "ORR")
#' 
#' @export 
#' 
select.dose.p23 = function(data=NULL, DCO1=16, dose_selection_endpoint = "ORR", ORRdiff=0){
   
   ######################
   #Stage 1
   ######################
   
   n.arms = length(unique(data$group))
   
   #1. Stage 1 data cut
   dat1cut = f.dataCut(data=data[data$stage == 1,], DCO=DCO1)
   
   n.events_IAD = t(table(dat1cut$group, dat1cut$cnsrCut))[1,] # add actual events at IAD YC================

   #2. Calculate z statistics and e

   #e1: total number of events for by DCO1 for each dose arm + control
   #z1: z statistics by DCO1 for each dose arm + control
   z1 = e1 = rep(NA, ncol=n.arms-1)
   if (dose_selection_endpoint == "ORR"){
     orr.control <- mean(dat1cut[dat1cut$group == 0, ]$response)
     orr.diff = rep(NA, n.arms-1)
    }
   
   for (j in 1:(n.arms-1)){
     datj = dat1cut[(dat1cut$group == 0 | dat1cut$group == j), ]
     z1[j] = logrank.one.sided(time=datj$survTimeCut, cnsr=datj$cnsrCut, group=datj$group)$z
     e1[j] = sum(1-datj$cnsrCut)
     
     if (dose_selection_endpoint == "ORR"){
       orr.diff[j] = mean(datj$response) - orr.control
     }   
   }
    
   #3. Dose selection
   if (dose_selection_endpoint == "ORR"){
     tmp = sort(orr.diff, index.return = TRUE)
     if(ORRdiff > 0 && orr.diff[1]==max(orr.diff) && ((orr.diff[1]-max(orr.diff[-1]))<=ORRdiff) ){
       # tolerate ORRdiff% difference between ORR of highest dose and lower doses: default first dose is the highest dose
         s=which.max(orr.diff[-1])+1
     }else{
       s = tmp$ix[n.arms-1]
     }
   } else {
     tmp = sort(z1, index.return = TRUE)
     s = tmp$ix[n.arms-1]
   }
   
   o = list()
   o$z1 = z1
   o$e1 = e1
   o$s = s
   o$n.events_IAD = n.events_IAD
   if (dose_selection_endpoint == "ORR") {o$orr.diff = orr.diff}
   
   return(o)
}
