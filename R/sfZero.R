#' Zero Interim Alpha Spending Function
#'
#' This function defines a custom alpha spending function for a group sequential design 
#' where exactly 0 alpha is spent at all interim analyses ($t < 1$), and the full 
#' overall alpha is spent at the final analysis ($t = 1$). 
#'
#' @param alpha Type I error (one-sided) to be spent overall, usually 0.025.
#' @param t A vector of points with values between 0 and 1 inclusive, indicating the information fraction or timing at each analysis.
#' @param param A placeholder parameter kept for compatibility with the `gsDesign` spending function template. It is not used in the calculation. Default is NULL.
#'
#' @return An object of class `spendfn` with values:
#' \describe{
#' \item{name}{The character string "Zero Interim Spending"}
#' \item{param}{The parameter value (NULL for this function)}
#' \item{parname}{The character string "none"}
#' \item{sf}{The function `sfZero` itself}
#' \item{spend}{A vector of cumulative alpha spent at each analysis}
#' \item{bound}{Boundary values (returns NULL as it is calculated by `gsDesign`)}
#' \item{prob}{Probabilities associated with bounds (returns NULL)}
#' }
#' 
#' @examples
#' # Example: 3-stage design using sfZero for the upper bound
#' library(gsDesign)
#' 
#' design <- gsDesign(
#'   k = 3, 
#'   test.type = 4, 
#'   alpha = 0.025, 
#'   beta = 0.1, 
#'   sfu = sfZero, 
#'   sfupar = NULL, 
#'   sfl = sfHSD, 
#'   sflpar = -2
#' )
#' 
#' print(design)
#' 
#' @export

sfZero <- function(alpha, t, param = NULL) {
  # 1. Cap information fraction at 1
  t[t > 1] <- 1
  
  # 2. Define the spending logic:
  # Get the total number of analyses passed to the function
  k <- length(t)
  
  # Create a vector of 0s for all analyses
  spend <- rep(0, k)
  
  # Force the VERY LAST analysis to spend the full alpha
  spend[k] <- alpha
  
  # 3. Construct the output list matching gsDesign's expectation
  x <- list(
    name = "Zero Interim Spending", 
    param = param, 
    parname = "none", 
    sf = sfZero, 
    spend = spend, 
    bound = NULL, 
    prob = NULL
  )
  
  class(x) <- "spendfn"
  return(x)
}
