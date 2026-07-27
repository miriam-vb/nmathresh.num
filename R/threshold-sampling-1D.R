#' Interval Bisection Sampling of Bias Thresholds (1D)
#' 
#' This function uses interval bisection to approximate decision-invariant
#' bias adjustment thresholds for network meta-analysis, with options to use 
#' either preset or user-defined decision functions. Thresholds represent 
#' the amount of adjustment needed in an individual data point before the 
#' treatment decision changes.
#'
#' @param data  Object (data frame, list, etc.) containing the NMA data that is 
#'    passed as an argument to decision_function
#' @param decision_function  Function accepting NMA data and bias adjustment
#'    used to implement the decision rule at each step of the boundary finding 
#'    method
#' @param indices  Numerical vector indicating the indices of the sequential list
#'    of data points for which bias thresholds will be estimated
#' @param admin  Administrative cutoff value for bias adjustment beyond which 
#'    decision invariance will not be assessed
#' @param tol  Tolerance for the absolute difference between converging boundary
#'    estimates
#' @param preset  Numeric value determining whether a specific preset 
#'    decision_function should be implemented rather than a user-supplied function
#' @param future_plan  Function call of future::plan that specifies how futures 
#'    are to be resolved. Defaults to sequential evaluation with one worker
#'
#' @return  List containing thresh.df, a data frame of thresholds and new 
#'    recommended treatments with columns \code{- Bias Thresh}, \code{- New Rec}, 
#'    \code{+ Bias Thresh}, and \code{+ New Rec}, best, the vector of optimal 
#'    treatments recommended by the decision function before bias adjustment,
#'    and args, a list of the arguments defined in the original function call
#'    
#' @import doFuture
#' @importFrom future plan
#' @importFrom foreach foreach
#' 
#' @export 


bias_thresh_1D <- function(data, decision_function = NULL, indices, admin = 10, 
                        tol = 10**(-3), preset = 1, future_plan = plan("sequential")) {
  
  # set decision_function to frequentist threshold analysis using the 
  # projection matrix with max efficacy as default
  if (preset == 1 || preset == 2) {
    n <- length(data$vec)
    decision_function <- function(data, bias = rep(0,n)) {
      if ("C" %in% names(data)) {
        # allow for arm-level bias assessment using the matrix mapping arms to
        # contrasts
        H <- solve(t(data$X)%*%(data$W)%*%(data$X))%*%t(data$X)%*%(data$W)%*%
          (data$C)
      } else {
        H <- solve(t(data$X)%*%(data$W)%*%(data$X))%*%t(data$X)%*%(data$W)
      }
      
      trtm_estimates <- H%*%(data$vec+bias)
      if (preset == 1){
        # maximal treatment effect is optimal
        if (all(trtm_estimates < 0)) {
          best <- c(1)
        } else {
          best <- c(which.max(trtm_estimates) + 1)
        }
      } else {
        # minimal treatment effect is optimal
        if (all(trtm_estimates > 0)) {
          best <- c(1)
        } else {
          best <- c(which.min(trtm_estimates) + 1)
        }
      }
      
      return(best)
    }
  }
  
  # output warning if an invalid preset was selected and no decision function
  # was supplied
  if (is.null(decision_function)) {
    stop("Invalid preset selected with no alternative decision function provided")
  }
  
  best <- decision_function(data)
  
  # implement random bias adjustment at each selected data point and use
  # IVT and interval bisection method to compute thresholds
  thresh.df <- data.frame(matrix(ncol=4,nrow=0))
  
  # define function for implementing threshold convergence
  thresh_conv <- function(ind) {
    
    b0 <- 0
    b1 <- admin/2
    b2 <- admin
    best0 <- best
    bvec <- rep(0,n)
    bvec[ind] <- b1
    best1 <- decision_function(data, bias = bvec)
    bvec[ind] <- b2
    best2 <- decision_function(data, bias = bvec)
    
    if (setequal(best2,best)) {
      # record administrative threshold if recommendation doesn't shift
      u <- admin
      trtU <- "Admin"
    } else {
      # iterate until biases are within tolerance
      while(abs(b2 - b0) > tol) {
        # select interval [b0,b1] or [b1,b2], then obtain midpoint of 
        # chosen interval and update variables
        if (!setequal(best0,best1)) {
          min <- b0
          max <- b1
          mid <- min + (b1 - b0)/2
          b0 <- min
          b1 <- mid
          b2 <- max
        } else if (!setequal(best1,best2)) {
          min <- b1
          max <- b2
          mid <- min + (b2 - b1)/2
          b0 <- min
          b1 <- mid
          b2 <- max
        }
        # update treatment recommendations for each point of bias
        bvec[ind] <- b0
        best0 <- decision_function(data, bias = bvec)
        bvec[ind] <- b1
        best1 <- decision_function(data, bias = bvec)
        bvec[ind] <- b2
        best2 <- decision_function(data, bias = bvec)
      }
      u <- b0
      trtU <- paste0(best2, collapse = ", ")
    }
    
    ## repeat for negative bias threshold
    
    b0 <- 0
    b1 <- -admin/2
    b2 <- -admin
    best0 <- best
    bvec <- rep(0,n)
    bvec[ind] <- b1
    best1 <- decision_function(data, bias = bvec)
    bvec[ind] <- b2
    best2 <- decision_function(data, bias = bvec)
    
    if (setequal(best2,best)) {
      # record negative administrative threshold if recommendation doesn't shift
      l <- -admin
      trtL <- "Admin"
    } else {
      # iterate until biases are within tolerance
      while(abs(b2 - b0) > tol) {
        # select interval [b0,b1] or [b1,b2], then obtain midpoint of 
        # chosen interval and update variables
        if (!setequal(best0,best1)){
          min <- b0
          max <- b1
          mid <- min + (b1 - b0)/2
          b0 <- min
          b1 <- mid
          b2 <- max
        } else if (!setequal(best1,best2)) {
          min <- b1
          max <- b2
          mid <- min + (b2 - b1)/2
          b0 <- min
          b1 <- mid
          b2 <- max
        }
        # update treatment recs for each point of bias
        bvec[ind] <- b0
        best0 <- decision_function(data, bias = bvec)
        bvec[ind] <- b1
        best1 <- decision_function(data, bias = bvec)
        bvec[ind] <- b2
        best2 <- decision_function(data, bias = bvec)
      }
      l <- b0
      trtL <- paste0(best2, collapse = ", ")
    }
    
    # store bias threshold and admin indicator/new superior treatment
    # report the point just inside the invariant region
    row <- c(l,trtL,u,trtU)
    return(row)
  }
  
  # allow for parallelization of boundary convergence method
  if (is.null(future_plan)) {
    with(plan("sequential"), local = TRUE)
  } else {
    with(future_plan, local = TRUE)
  }
  thresh <- foreach(ind = indices, .options.future = 
                    list(seed = TRUE, package = structure(TRUE, 
                    add = c("nmathresh.num")))) %dofuture% {
    thresh_conv(ind)
  }
  # reform the data.frame using futures
  thresh.df <- t(as.data.frame(thresh))
  rownames(thresh.df) <- NULL
  thresh.df <- as.data.frame(thresh.df)
  

  # convert bias columns to numeric
  thresh.df[, c(1,3)] <- apply(thresh.df[, c(1,3)], 2, 
                               function(x) as.numeric(as.character(x)))
  
  # return bias values and index of treatment that became superior at each switch 
  # (or that it was admin cutoff, indicating a potential invariant region)
  colnames(thresh.df) <- c("- Bias Thresh", "- New Rec", 
                           "+ Bias Thresh", "+ New Rec")
  return(list(thresh.df = thresh.df, best = best, 
              args = list(data = data, decision_function = decision_function, 
                          indices = indices, admin = admin, tol = tol, 
                          preset = preset, future_plan = future_plan)))
}


