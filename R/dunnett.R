#' Adjusted p value for testing H_J using Dunnett method
#'
#'This functions provides the adjusted p value for testing a family-wise hypothesis H_J based on the p values of individual raw p values in J.
#'
#' @param p A vector of individual raw p values in J. 
#'
#' @return Adjusted p values using Dunnett procedure
#' 
#' @importFrom mvtnorm pmvnorm
#' @importFrom stats qnorm
#' @export
#'
#' @examples
#' p = c(0.01, 0.02, 0.03, 0.013)
#' dunnett(p)
dunnett <- function(p, cr=NULL){
  num_trt <- length(p)
  # correlation
  if(is.null(cr)) {
    cr <- matrix(0.5, num_trt, num_trt)
    diag(cr) <- 1
  }
  # equal weight
  w <- rep(1, num_trt)/num_trt
  padjusted <- p.dunnet.ph23(p = p, cr = cr,w = w, upscale = FALSE)
  ans <- min(min(padjusted), 1)
  return(ans)
}

# Build a Dunnett correlation matrix for logrank Z's from event counts.

# nE_trt[i] = #events in treatment arm i at its DCO (stage 1 subjects here)

# nE_con[i] = #events in control arm at the same DCO used for comparison i

build_cor_from_events <- function(nE_trt, nE_con, make_PD = TRUE, tol = 1e-10) {
  
  stopifnot(length(nE_trt) == length(nE_con))
  
  k <- length(nE_trt)
  
  # Denominator term (treatment + control events) per comparison
  
  denom <- nE_trt + nE_con
  
  # Fraction of shared control information between comparisons i and j:
  
  # nested-cutoff approximation => shared control events = min(n0_i, n0_j)
  
  frac_shared <- outer(nE_con, nE_con, pmin) / outer(nE_con, nE_con, pmax)
  
  # Core Dunnett-type correlation term due to sharing the same control:
  
  #   (n_i n_j)/((n_i+n0_i)(n_j+n0_j))
  
  core <- outer(nE_trt, nE_trt) / outer(denom, denom)
  
  # Combine the two components and take square root to get correlations
  
  cr <- sqrt(frac_shared * core)
  
  # Guard against 0/0 or Inf from zero event counts
  
  cr[!is.finite(cr)] <- 0
  
  # Enforce correlation matrix structure
  
  diag(cr) <- 1
  
  cr <- (cr + t(cr)) / 2  # exact symmetry
  
  # Optional: enforce positive-definiteness (PD) for MVN calculations
  
  # (often needed for mvtnorm / multivariate Dunnett computations)
  
  if (make_PD) {
    
    ev_min <- min(eigen(cr, symmetric = TRUE, only.values = TRUE)$values)
    
    if (ev_min < -tol) {
      
      if (requireNamespace("Matrix", quietly = TRUE)) {
        
        cr <- as.matrix(Matrix::nearPD(cr, corr = TRUE)$mat)
        
        diag(cr) <- 1
        
        cr <- (cr + t(cr)) / 2
        
      } else {
        
        warning("Matrix not installed; returning a possibly non-PD correlation matrix.")
        
      }
      
    }
    
  }
  
  cr
  
}


## direct copy of gMCP:::p.dunnet, but added seed = 20240716 for pmvnorm so that the results are reproducible. 
p.dunnet.ph23 <- function(p,cr,w,upscale, alternatives="less"){
  if(length(cr)>1){
    conn <- conn.comp.ph23(cr)
  } else {
    conn <- 1
  }
  twosided <- alternatives==rep("two.sided", length(w))
  lconn <- sapply(conn,length)
  conn <- lapply(conn,as.numeric)
  e <- sapply(1:length(p),function(i){
    sum(sapply(conn,function(edx){
      if(length(edx)>1){
        if (upscale=="o3") {
          return((1-mvtnorm::pmvnorm(
            lower=ifelse(twosided[edx],qnorm(pmin(1,(w[edx]*p[i]/(w[i]*sum(w))))/2),-Inf),
            upper=ifelse(twosided[edx],qnorm(1-pmin(1,(w[edx]*p[i]/(w[i]*sum(w))))/2),qnorm(1-pmin(1,(w[edx]*p[i]/(w[i]*sum(w)))))),
            corr=cr[edx,edx], seed = 20240716, abseps=10^-5)))
        } else {
          return((1-mvtnorm::pmvnorm(
            lower=ifelse(twosided[edx],qnorm(pmin(1,(w[edx]*p[i]/(w[i])))/2),-Inf),
            upper=ifelse(twosided[edx],qnorm(1-pmin(1,(w[edx]*p[i]/(w[i])))/2),qnorm(1-pmin(1,(w[edx]*p[i]/(w[i]))))),
            corr=cr[edx,edx], seed = 20240716,abseps=10^-5))/ifelse(upscale,1,sum(w)))
        }
      } else {
        if(upscale=="o3" || !upscale){
          return((w[edx]*p[i]/(w[i]*sum(w))))
        } else {
          return((w[edx]*p[i]/(w[i])))
        }
      }
    }))})
  
  e <- pmin(e,1)
  e
}

## direct copy of gMCP:::conn.comp
conn.comp.ph23 <- function(m){
  N <- 1:ncol(m)
  M <- numeric(0)
  out <- list()
  while(length(N)>0){
    Q <- setdiff(N,M)[1]
    while(length(Q)>0){
      w <- Q[1]
      M <- c(M,w)
      Q <- setdiff(unique(c(Q,which(!is.na(m[w,])))),M)
    }
    out <- c(out,list(M))
    N <- setdiff(N,M)
    M <- numeric(0)
  }
  return(out)
}