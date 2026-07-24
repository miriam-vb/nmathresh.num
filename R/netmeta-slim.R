#' ----------------------------------------------------------------------------
#' Network meta-analysis using a graph-theoretical estimation method
#' 
#' The netmeta function from the netmeta package (authors Gerta Rücker 
#' \email{gerta.ruecker@@uniklinik-freiburg.de}, Guido Schwarzer 
#' \email{guido.schwarzer@@uniklinik-freiburg.de}) employs a graph-theoretical
#' approach to network meta-analysis, which is equivalent to the Frequentist 
#' approach based on weighted least squares regression. This slimmed version
#' of the function was copied from the netmeta package and modified to focus 
#' on efficiently estimating the relative treatment effect estimates and their 
#' standard errors for both a fixed effects and random effects model.
#' 
#' Source: https://github.com/guido-s/netmeta/blob/develop/R/netmeta.R
#' License: GPL (>=2)
#'
#' @param TE  Vector containing study-level estimated treatment effects
#' @param seTE  Vector containing standard errors of study-level treatment effects
#' @param treat1  Vector indicating the label/number corresponding to the first
#'   treatment of the pairwise comparison
#' @param treat2  Vector indicating the label/number corresponding to the second
#'   treatment of the pairwise comparison
#' @param studlab  Vector containing study labels
#' @param data  Optional data frame containing the study information
#' @param subset  Optional vector specifying a subset of studies to be used
#' @param correlated  Optional logical vector specifying whether treatment arms 
#'    of a multi-arm study are correlated (see Details of the netmeta function)
#' @param sm  Character string indicating underlying summary measure, e.g., "RD",
#'    "RR", "OR", "ASD", "HR", "MD", "SMD", or "ROM"
#' @param reference.group  Reference treatment (first treatment used as default)
#' @param baseline.reference  Logical indicating whether results should be 
#'    expressed as comparisons of other treatments versus the reference treatment
#'    (default) or vice versa. Only considered if reference.group is specified
#' @param all.treatments  Logical or NULL value. If TRUE, matrices with all 
#'    treatment effects will be printed
#' @param seq  Character or numerical vector specifying the sequence of treatments
#'    in printouts
#' @param tau.preset  Optional value for manually setting the square root of the
#'    between-study variance, \eqn{\tau^2}
#' @param tol.multiarm  Numeric for the tolerance for consistency of treatment
#'    estimates in multi-arm studies which are consistent by design
#' @param tol.multiarm.se  Numeric for the tolerance for consistency of standard 
#'    errors in multi-arm studies which are consistent by design. This check is 
#'    not conducted if the argument is NULL
#' @param details.chkmultiarm  Logical indicating whether treatment estimates and/
#'    or variances of multi-arm studies with inconsistent results or negative 
#'    multi-arm variances should be printed
#' @param sep.trts  Character used in comparison names as separator between 
#'    treatment labels
#' @param nchar.trts  Numeric defining the minimum number of characters used to
#'    create unique treatment names (see Details of the netmeta function)
#' @param nchar.studlab  Numeric definign the minimum number of characters used 
#'    to create unique study labels
#' @param func.inverse  R function used to calculate the pseudoinverse (defaults
#'    to the Moore-Penrose pseudoinverse) of the Laplacian matrix L
#' @param n1  Number of observations in first treatment group
#' @param n2  Number of observations in second treatment group
#' @param event1  Number of events in first treatment group
#' @param event2  Number of events in second treatment group
#' @param mean1  Mean in first treatment group
#' @param mean2  Mean in second treatment group
#' @param sd1  Standard deviation in first treatment group
#' @param sd2  Standard deviation in second treatment group
#' @param time1  Person time at risk in first treatment group
#' @param time2  Person time at risk in second treatment group
#' @param warn  Logical indicating if warnings should be printed (e.g., if
#'    studies are excluded from meta-analysis due to zero standard errors)
#'
#' @return  List containing studlab, treat1, treat2, TE and seTE as passed to 
#'    the function, as well as TE.common and TE.random, \eqn{n \times n} matrices
#'    (\eqn{n} representing the total number of treatments) containing the 
#'    estimated overall treatment effects from the common/random effects model, 
#'    respectively, along with seTE.common and seTE.random, \eqn{n \times n} 
#'    matrices containing the estimated standard errors for these estimates.
#'    
#' @export
#' ----------------------------------------------------------------------------


netmeta_slim <- function (TE, seTE, treat1, treat2, studlab, data = NULL, 
                          subset = NULL, correlated, sm, 
                          reference.group, baseline.reference = TRUE, 
                          all.treatments = NULL, seq = NULL, 
                          tau.preset = NULL, tol.multiarm = 0.001, 
                          tol.multiarm.se = NULL, details.chkmultiarm = FALSE, 
                          sep.trts = ":", nchar.trts = 666, nchar.studlab = 666, 
                          func.inverse = invmat, n1 = NULL, n2 = NULL, 
                          event1 = NULL, event2 = NULL, mean1 = NULL, 
                          mean2 = NULL, sd1 = NULL, sd2 = NULL, 
                          time1 = NULL, time2 = NULL, warn = TRUE) 
{

  missing.reference.group <- missing(reference.group)
  baseline.reference <- replaceNULL(baseline.reference, TRUE)
  chklogical(baseline.reference)
  if (!is.null(all.treatments)) 
    chklogical(all.treatments)
  method.tau <- "DL"
  if (!is.null(tau.preset)) 
    chknumeric(tau.preset, min = 0, length = 1)
  tol.multiarm <- replaceNULL(tol.multiarm, 0.001)
  chknumeric(tol.multiarm, min = 0, length = 1)
  if (!is.null(tol.multiarm.se)) 
    chknumeric(tol.multiarm.se, min = 0, length = 1)
  details.chkmultiarm <- replaceNULL(details.chkmultiarm, FALSE)
  chklogical(details.chkmultiarm)
  missing.sep.trts <- missing(sep.trts)
  sep.trts <- replaceNULL(sep.trts, ":")
  chkchar(sep.trts, length = 1)
  nchar.studlab <- replaceNULL(nchar.studlab, 666)
  chknumeric(nchar.studlab, length = 1)
  keepdata = TRUE
  chklogical(warn)
  chklogical(baseline.reference)
  missing.nchar.trts <- missing(nchar.trts)
  nchar.trts <- replaceNULL(nchar.trts, 666)
  chknumeric(nchar.trts, min = 1, length = 1)
  nulldata <- is.null(data)
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  if (nulldata) 
    data <- sfsp
  TE <- catch("TE", mc, data, sfsp)
  avail.reference.group.pairwise <- FALSE
  if (inherits(TE, "pairwise")) {
    is.pairwise <- TRUE
    sm <- attr(TE, "sm")
    allstudies <- replaceNULL(attr(TE, "allstudies"), TRUE)
    if (missing.reference.group) {
      reference.group <- attr(TE, "reference.group")
      if (is.null(reference.group)) 
        reference.group <- ""
      else avail.reference.group.pairwise <- TRUE
    }
    keep.all.comparisons <- attr(TE, "keep.all.comparisons")
    if (!is.null(keep.all.comparisons) && !keep.all.comparisons) 
      stop("First argument is a pairwise object created with ", 
           "'keep.all.comparisons = FALSE'.", call. = TRUE)
    if (is.null(attr(TE, "varnames"))) 
      seTE <- TE$seTE
    else seTE <- TE[[attr(TE, "varnames")[2]]]
    treat1 <- TE$treat1
    treat2 <- TE$treat2
    studlab <- TE$studlab
    if (!is.null(TE$n1)) 
      n1 <- TE$n1
    if (!is.null(TE$n2)) 
      n2 <- TE$n2
    if (!is.null(TE$event1)) 
      event1 <- TE$event1
    if (!is.null(TE$event2)) 
      event2 <- TE$event2
    incr1 <- TE$incr1
    incr2 <- TE$incr2
    if (!is.null(incr2)) 
      incr <- incr2
    else incr <- TE$incr
    if (!is.null(TE$mean1)) 
      mean1 <- TE$mean1
    if (!is.null(TE$mean2)) 
      mean2 <- TE$mean2
    if (!is.null(TE$sd1)) 
      sd1 <- TE$sd1
    if (!is.null(TE$sd2)) 
      sd2 <- TE$sd2
    if (!is.null(TE$time1)) 
      time1 <- TE$time1
    if (!is.null(TE$time2)) 
      time2 <- TE$time2
    pairdata <- TE
    data <- TE
    if (is.null(attr(TE, "varnames"))) 
      TE <- TE$TE
    else TE <- TE[[attr(TE, "varnames")[1]]]
  }
  else {
    is.pairwise <- FALSE
    if (missing(sm)) 
      if (!is.null(data) && !is.null(attr(data, "sm"))) 
        sm <- attr(data, "sm")
    else sm <- ""
    seTE <- catch("seTE", mc, data, sfsp)
    treat1 <- catch("treat1", mc, data, sfsp)
    treat2 <- catch("treat2", mc, data, sfsp)
    studlab <- catch("studlab", mc, data, sfsp)
    n1 <- catch("n1", mc, data, sfsp)
    n2 <- catch("n2", mc, data, sfsp)
    event1 <- catch("event1", mc, data, sfsp)
    event2 <- catch("event2", mc, data, sfsp)
    incr1 <- NULL
    incr2 <- catch("incr", mc, data, sfsp)
    mean1 <- catch("mean1", mc, data, sfsp)
    mean2 <- catch("mean2", mc, data, sfsp)
    sd1 <- catch("sd1", mc, data, sfsp)
    sd2 <- catch("sd2", mc, data, sfsp)
    time1 <- catch("time1", mc, data, sfsp)
    time2 <- catch("time2", mc, data, sfsp)
  }
  chknumeric(TE)
  chknumeric(seTE)
  if (!any(!is.na(TE) & !is.na(seTE))) 
    stop("Missing data for estimates (argument 'TE') and ", 
         "standard errors (argument 'seTE') in all studies.\n  ", 
         "No network meta-analysis possible.", call. = FALSE)
  k.Comp <- length(TE)
  if (is.factor(treat1)) 
    treat1 <- as.character(treat1)
  if (is.factor(treat2)) 
    treat2 <- as.character(treat2)
  treat1 <- rmSpace(rmSpace(treat1, end = TRUE))
  treat2 <- rmSpace(rmSpace(treat2, end = TRUE))
  if (length(studlab) == 0) {
    if (warn) 
      warning("No information given for argument 'studlab'. ", 
              "Assuming that comparisons are from independent studies.", 
              call. = FALSE)
    studlab <- seq(along = TE)
  }
  studlab <- as.character(studlab)
  subset <- catch("subset", mc, data, sfsp)
  missing.subset <- is.null(subset)
  correlated <- catch("correlated", mc, data, sfsp)
  if (is.null(correlated)) 
    correlated <- FALSE
  if (!is.logical(correlated)) 
    stop("Argument 'correlated' must be a logical vector.", 
         call. = FALSE)
  if (length(correlated) == 1) 
    correlated <- rep(correlated, length(TE))
  else if (length(correlated) != length(TE)) 
    stop("Different length for arguments 'TE' and 'correlated'.", 
         call. = FALSE)
  if (!is.null(event1) & !is.null(event2)) 
    available.events <- TRUE
  else available.events <- FALSE
  if (!is.null(n1) & !is.null(n2)) 
    available.n <- TRUE
  else available.n <- FALSE
  if (available.events & is.null(incr2)) 
    incr2 <- rep(0, length(event2))
  if (!is.null(mean1) & !is.null(mean2)) 
    available.means <- TRUE
  else available.means <- FALSE
  if (!is.null(sd1) & !is.null(sd2)) 
    available.sds <- TRUE
  else available.sds <- FALSE
  if (!is.null(time1) & !is.null(time2)) 
    available.times <- TRUE
  else available.times <- FALSE
  if (keepdata) {
    if (nulldata & !is.pairwise) 
      data <- data.frame(.studlab = studlab, stringsAsFactors = FALSE)
    else if (nulldata & is.pairwise) {
      data <- pairdata
      data$.studlab <- studlab
    }
    else data$.studlab <- studlab
    data$.order <- seq_along(studlab)
    data$.treat1 <- treat1
    data$.treat2 <- treat2
    data$.TE <- TE
    data$.seTE <- seTE
    data$.correlated <- correlated
    data$.event1 <- event1
    data$.n1 <- n1
    data$.event2 <- event2
    data$.n2 <- n2
    data$.incr1 <- incr1
    data$.incr2 <- incr2
    data$.mean1 <- mean1
    data$.sd1 <- sd1
    data$.mean2 <- mean2
    data$.sd2 <- sd2
    data$.time1 <- time1
    data$.time2 <- time2
    wo <- data$.treat1 > data$.treat2
    if (any(wo)) {
      data$.TE[wo] <- -data$.TE[wo]
      ttreat1 <- data$.treat1
      data$.treat1[wo] <- data$.treat2[wo]
      data$.treat2[wo] <- ttreat1[wo]
      if (isCol(data, ".n1") & isCol(data, ".n2")) {
        tn1 <- data$.n1
        data$.n1[wo] <- data$.n2[wo]
        data$.n2[wo] <- tn1[wo]
      }
      if (isCol(data, ".event1") & isCol(data, ".event2")) {
        tevent1 <- data$.event1
        data$.event1[wo] <- data$.event2[wo]
        data$.event2[wo] <- tevent1[wo]
      }
      if (isCol(data, ".mean1") & isCol(data, ".mean2")) {
        tmean1 <- data$.mean1
        data$.mean1[wo] <- data$.mean2[wo]
        data$.mean2[wo] <- tmean1[wo]
      }
      if (isCol(data, ".sd1") & isCol(data, ".sd2")) {
        tsd1 <- data$.sd1
        data$.sd1[wo] <- data$.sd2[wo]
        data$.sd2[wo] <- tsd1[wo]
      }
      if (isCol(data, ".time1") & isCol(data, ".time2")) {
        ttime1 <- data$.time1
        data$.time1[wo] <- data$.time2[wo]
        data$.time2[wo] <- ttime1[wo]
      }
    }
    if (!missing.subset) {
      if (length(subset) == dim(data)[1]) 
        data$.subset <- subset
      else {
        data$.subset <- FALSE
        data$.subset[subset] <- TRUE
      }
    }
  }
  if (!missing.subset) {
    if ((is.logical(subset) & (sum(subset) > k.Comp)) || 
        (length(subset) > k.Comp)) 
      stop("Length of subset is larger than number of studies.", 
           call. = FALSE)
    TE <- TE[subset]
    seTE <- seTE[subset]
    treat1 <- treat1[subset]
    treat2 <- treat2[subset]
    studlab <- studlab[subset]
    correlated <- correlated[subset]
    if (!is.null(n1)) 
      n1 <- n1[subset]
    if (!is.null(n2)) 
      n2 <- n2[subset]
    if (!is.null(event1)) 
      event1 <- event1[subset]
    if (!is.null(event2)) 
      event2 <- event2[subset]
    if (!is.null(incr1)) 
      incr1 <- incr1[subset]
    if (!is.null(incr2)) 
      incr2 <- incr2[subset]
    if (!is.null(mean1)) 
      mean1 <- mean1[subset]
    if (!is.null(mean2)) 
      mean2 <- mean2[subset]
    if (!is.null(sd1)) 
      sd1 <- sd1[subset]
    if (!is.null(sd2)) 
      sd2 <- sd2[subset]
    if (!is.null(time1)) 
      time1 <- time1[subset]
    if (!is.null(time2)) 
      time2 <- time2[subset]
  }
  labels <- sort(unique(c(treat1, treat2)))
  if (!is.null(seq)) 
    seq <- setseq(seq, labels)
  else {
    seq <- labels
    if (is.numeric(seq)) 
      seq <- as.character(seq)
  }
  sep.trts <- setsep(labels, sep.trts, missing = missing.sep.trts)
  if (any(treat1 == treat2)) 
    stop("Treatments must be different (arguments 'treat1' and 'treat2').", 
         call. = FALSE)
  tabnarms <- table(studlab)
  sel.narms <- !is_wholenumber((1 + sqrt(8 * tabnarms + 1))/2)
  if (sum(sel.narms) == 1) 
    stop("Study '", names(tabnarms)[sel.narms], "' has a wrong number of comparisons.", 
         "\n  Please provide data for all treatment comparisons", 
         " (two-arm: 1; three-arm: 3; four-arm: 6, ...).", 
         call. = FALSE)
  if (sum(sel.narms) > 1) 
    stop("The following studies have a wrong number of comparisons: ", 
         paste(paste0("'", names(tabnarms)[sel.narms], "'"), 
               collapse = ", "), "\n  Please provide data for all treatment comparisons", 
         " (two-arm: 1; three-arm: 3; four-arm: 6, ...).", 
         call. = FALSE)
  n.subnets <- netconnection(treat1, treat2, studlab)$n.subnets
  if (n.subnets > 1) 
    stop("Network consists of ", n.subnets, " separate sub-networks.\n  ", 
         "Use R function 'netconnection' to identify sub-networks.", 
         call. = FALSE)
  excl <- is.na(TE) | is.na(seTE) | seTE <= 0
  if (any(excl)) {
    if (keepdata) {
      if (!missing.subset) {
        data$.excl <- NA
        data$.excl[subset] <- excl
      }
      else data$.excl <- excl
    }
    dat.NAs <- data.frame(studlab = studlab[excl], treat1 = treat1[excl], 
                          treat2 = treat2[excl], TE = format(round(TE[excl], 
                                                                   4)), seTE = format(round(seTE[excl], 4)), stringsAsFactors = FALSE)
    studlab <- studlab[!excl]
    treat1 <- treat1[!excl]
    treat2 <- treat2[!excl]
    TE <- TE[!excl]
    seTE <- seTE[!excl]
    correlated <- correlated[!excl]
    if (!is.null(n1)) 
      n1 <- n1[!excl]
    if (!is.null(n2)) 
      n2 <- n2[!excl]
    if (!is.null(event1)) 
      event1 <- event1[!excl]
    if (!is.null(event2)) 
      event2 <- event2[!excl]
    if (!is.null(incr1)) 
      incr1 <- incr1[!excl]
    if (!is.null(incr2)) 
      incr2 <- incr2[!excl]
    if (!is.null(mean1)) 
      mean1 <- mean1[!excl]
    if (!is.null(mean2)) 
      mean2 <- mean2[!excl]
    if (!is.null(sd1)) 
      sd1 <- sd1[!excl]
    if (!is.null(sd2)) 
      sd2 <- sd2[!excl]
    if (!is.null(time1)) 
      time1 <- time1[!excl]
    if (!is.null(time2)) 
      time2 <- time2[!excl]
    seq <- seq[seq %in% unique(c(treat1, treat2))]
    labels <- labels[labels %in% unique(c(treat1, treat2))]
  }
  tabnarms <- table(studlab)
  sel.narms <- !is_wholenumber((1 + sqrt(8 * tabnarms + 1))/2)
  if (sum(sel.narms) == 1) 
    stop("After removing comparisons with missing treatment effects", 
         " or standard errors,\n  study '", names(tabnarms)[sel.narms], 
         "' has a wrong number of comparisons.\n  ", "Please check data and consider to ", 
         if (is.pairwise & !allstudies) 
           "\n   (i) ", "remove study from network meta-analysis", 
         if (is.pairwise & !allstudies) 
           " or\n  (ii) use argument 'allstudies = TRUE' in pairwise()", 
         ".", call. = FALSE)
  if (sum(sel.narms) > 1) 
    stop("After removing comparisons with missing treatment effects", 
         " or standard errors,\n  the following studies have", 
         " a wrong number of comparisons: ", paste(paste0("'", 
                                                          names(tabnarms)[sel.narms], "'"), collapse = ", "), 
         "\n  ", "Please check data and consider to ", if (is.pairwise & 
                                                           !allstudies) 
           "\n   (i) ", "remove studies from network meta-analysis", 
         if (is.pairwise & !allstudies) 
           " or\n  (ii) use argument 'allstudies = TRUE' in pairwise()", 
         ".", call. = FALSE)
  if (any(excl)) {
    if (warn) 
      warning("Comparison", if (sum(excl) > 1) 
        "s", " with missing TE / seTE or zero seTE not considered ", 
        "in network meta-analysis.", call. = FALSE)
    if (warn) {
      cat("Comparison", if (sum(excl) > 1) 
        "s", " not considered in network meta-analysis:\n", 
        sep = "")
      prmatrix(dat.NAs, quote = FALSE, right = TRUE, rowlab = rep("", 
                                                                  sum(excl)))
    }
  }
  n.subnets <- netconnection(treat1, treat2, studlab)$n.subnets
  if (n.subnets > 1) 
    stop("After removing comparisons with missing treatment effects", 
         " or standard errors,\n  network consists of ", n.subnets, 
         " separate sub-networks.\n  ", "Please check data and consider to remove studies", 
         " from network meta-analysis.", call. = FALSE)
  wo <- treat1 > treat2
  if (any(wo)) {
    TE[wo] <- -TE[wo]
    ttreat1 <- treat1
    treat1[wo] <- treat2[wo]
    treat2[wo] <- ttreat1[wo]
    if (available.n) {
      tn1 <- n1
      n1[wo] <- n2[wo]
      n2[wo] <- tn1[wo]
    }
    if (available.events) {
      tevent1 <- event1
      event1[wo] <- event2[wo]
      event2[wo] <- tevent1[wo]
    }
    if (available.means) {
      tmean1 <- mean1
      mean1[wo] <- mean2[wo]
      mean2[wo] <- tmean1[wo]
    }
    if (available.sds) {
      tsd1 <- sd1
      sd1[wo] <- sd2[wo]
      sd2[wo] <- tsd1[wo]
    }
    if (available.times) {
      ttime1 <- time1
      time1[wo] <- time2[wo]
      time2[wo] <- ttime1[wo]
    }
  }
  if (missing.reference.group & !avail.reference.group.pairwise) {
    go.on <- TRUE
    i <- 0
    while (go.on) {
      i <- i + 1
      sel.i <- !is.na(TE) & !is.na(seTE) & (treat1 == labels[i] | 
                                              treat2 == labels[i])
      if (sum(sel.i) > 0) {
        go.on <- FALSE
        reference.group <- labels[i]
      }
      else if (i == length(labels)) {
        go.on <- FALSE
        reference.group <- ""
      }
    }
  }
  if (is.null(all.treatments)) 
    if (reference.group == "") 
      all.treatments <- TRUE
  else all.treatments <- FALSE
  if (reference.group != "") 
    reference.group <- setref(reference.group, labels)
  p0 <- prepare(TE, seTE, treat1, treat2, studlab, correlated = correlated, 
                func.inverse = func.inverse)
  W.matrix.common <- as.matrix(p0$W)
  Cov.matrix.common <- as.matrix(p0$Cov)
  dat.c <- p0$data
  chkmultiarm(dat.c$TE, dat.c$seTE, dat.c$treat1, dat.c$treat2, 
      dat.c$studlab, dat.c$correlated, tol.multiarm = tol.multiarm, 
      tol.multiarm.se = tol.multiarm.se, details = details.chkmultiarm)
  tdata <- data.frame(studies = dat.c$studlab, narms = dat.c$narms, 
                      order = dat.c$order, stringsAsFactors = FALSE)
  tdata <- tdata[!duplicated(tdata[, c("studies", "narms")]), 
                 , drop = FALSE]
  studies <- tdata$studies[order(tdata$order)]
  narms <- tdata$narms[order(tdata$order)]
  res.c <- nma_ruecker_slim(dat.c$TE, W.matrix.common, sqrt(1/dat.c$weights), 
                            dat.c$treat1, dat.c$treat2, dat.c$treat1.pos, 
                            dat.c$treat2.pos, dat.c$narms, dat.c$studlab, sm, 
                            dat.c$seTE, 0, sep.trts, method.tau, 
                            func.inverse = func.inverse, Cov0 = p0$Cov)
  if (is.null(tau.preset)) {
    tau <- res.c$tau
  }
  else {
    tau <- tau.preset
  }
  p1 <- prepare(TE, seTE, treat1, treat2, studlab, tau = tau, 
                correlated = correlated, func.inverse = func.inverse)
  W.matrix.random <- as.matrix(p1$W)
  Cov.matrix.random <- as.matrix(p1$Cov)
  dat.r <- p1$data
  res.r <- nma_ruecker_slim(dat.r$TE, W.matrix.random, sqrt(1/dat.r$weights), 
                            dat.r$treat1, dat.r$treat2, dat.r$treat1.pos, 
                            dat.r$treat2.pos, dat.r$narms, dat.r$studlab, sm, 
                            dat.r$seTE, tau, sep.trts, method.tau, 
                            func.inverse = func.inverse, Cov0 = p1$Cov)
  
  res <- list(studlab = res.c$studlab, treat1 = res.c$treat1, 
              treat2 = res.c$treat2, TE = res.c$TE, seTE = res.c$seTE.orig, 
              TE.common = res.c$TE.pooled, seTE.common = res.c$seTE.pooled, 
              TE.random = res.r$TE.pooled, seTE.random = res.r$seTE.pooled)
  
  return(res)
}


#
# Modified version of auxiliary function nma_ruecker from netmeta
#
# Source: https://github.com/guido-s/netmeta/blob/develop/R/nma.ruecker.R
#
nma_ruecker_slim <- function(TE, W, seTE,
                             treat1, treat2,
                             treat1.pos, treat2.pos,
                             narms, studlab,
                             sm = "",
                             seTE.orig, tau.direct = 0, sep.trts = ":",
                             method.tau = "DL",
                             func.inverse,
                             Cov0) {
  
  require(MASS)
  w.pooled <- 1 / seTE^2
  
  m <- length(TE)                        # Number of pairwise comparisons (edges)
  n <- length(unique(c(treat1, treat2))) # Number of treatments (vertices)
  df1 <- 2 * sum(1 / narms)              # Sum of degrees of freedom per study
  df.Q <- df1 - (n - 1)                  # Degrees of freedom for Q test
  
  # Drop Matrix attributes
  W <- as.matrix(W)
  class(W) <- "matrix"
  #
  Cov0 <- as.matrix(Cov0)
  class(Cov0) <- "matrix"
  
  ##
  ## B is the edge-vertex incidence matrix (m x n)
  ##
  B <- createB(treat1.pos, treat2.pos, ncol = n)
  ##
  ## B.full is the full edge-vertex incidence matrix (m x n)
  ##
  B.full <- createB(ncol = n)
  ##
  ## M is the unweighted Laplacian, D its diagonal,
  ## A is the adjacency matrix
  ##
  M <- t(B) %*% B    # unweighted Laplacian matrix
  D <- diag(diag(M)) # diagonal matrix
  A <- D - M         # adjacency matrix (n x n)
  ##
  ## L is the weighted Laplacian (Kirchhoff) matrix (n x n)
  ## Lplus is its Moore-Penrose pseudoinverse
  ##
  L <- t(B) %*% W %*% B
  Lplus <- do.call(func.inverse, list(X = L))
  Lplus[is_zero(Lplus)] <- 0
  #
  L1 <- t(B) %*% ginv(Cov0) %*% B
  Lplus1 <- do.call(func.inverse, list(X = L1))
  Lplus1[is_zero(Lplus1)] <- 0
  #
  # R resistance distance (variance) matrix (n x n)
  #
  R <- matrix(0, nrow = n, ncol = n) 
  for (i in 1:n) {
    for (j in 1:n) {
      R[i, j] <- Lplus1[i, i] + Lplus1[j, j] - 2 * Lplus1[i, j]
    }
  }
  ##
  ## V is the vector of effective variances
  ##
  V <- vector(length = m, mode = "numeric")
  for (i in 1:m) {
    V[i] <- R[treat1.pos[i], treat2.pos[i]]
  }
  ##
  ## G is the matrix B %*% Lplus %*% t(B)
  ## H is the projection matrix (also called "hat matrix")
  ##
  ## Interpretation:
  ## (i)    diag(G) = V                 The effective variances
  ## (ii)   diag(H) = V %*% W = V * w   The leverages
  ## (iii)  sum(diag(H)) = n - 1        Rank of projection
  ## (iv)   mean(diag(H)) = (n - 1) / m Mean leverage = average efficiency
  ##
  G <- B %*% Lplus %*% t(B)
  H <- G %*% W
  ##
  ## Variance-covariance matrix for all comparisons
  ##
  Cov <- B.full %*% Lplus1 %*% t(B.full)
  ##
  ## Resulting effects and variances at numbered edges
  ##
  v <- as.vector(H %*% TE)
  ##
  ## Resulting effects, all edges, as a n x n matrix:
  ##
  all <- matrix(NA, nrow = n, ncol = n)
  ##
  for (i in 1:m) {
    all[treat1.pos[i], treat2.pos[i]] <- v[i]
  }
  ##
  for (i in 1:n) {
    for (j in 1:n) {
      for (k in 1:n) {
        if (!is.na(all[i, k]) & !is.na(all[j, k])) {
          all[i, j] <- all[i, k] - all[j, k]
          all[j, i] <- all[j, k] - all[i, k]
        }
        if (!is.na(all[i, j]) & !is.na(all[k, j])) {
          all[i, k] <- all[i, j] - all[k, j]
          all[k, i] <- all[k, j] - all[i, j]
        }
        if (!is.na(all[i, k]) & !is.na(all[i, j])) {
          all[j, k] <- all[i, k] - all[i, j]
          all[k, j] <- all[i, j] - all[i, k]
        }
      }
    }
  }
  ##
  ## Test of total heterogeneity / inconsistency:
  ##
  Q <- as.vector(t(TE - v) %*% W %*% (TE - v))
  if (df.Q == 0)
    pval.Q <-  NA
  else
    pval.Q <- pchisq(Q, df.Q, lower.tail = FALSE)
  ##
  ## Heterogeneity variance
  ##
  I <- diag(m)
  E <- matrix(0, nrow = m, ncol = m)
  for (i in 1:m)
    for (j in 1:m)
      E[i, j] <- as.numeric(studlab[i] == studlab[j])
  ##
  if (df.Q == 0) {
    tau2 <- NA
    tau <- NA
    
  }
  else {
    tau2 <-
      max(0, (Q - df.Q) / sum(diag((I - H) %*% (B %*% t(B) * E / 2) %*% W)))
    tau <- sqrt(tau2)
    
  }
  
  
  ##
  names.treat <- sort(unique(c(treat1, treat2)))
  ##
  
  TE.pooled <- all
  seTE.pooled <- sqrt(R)
  
  #
  rownames(TE.pooled) <- colnames(TE.pooled) <- names.treat
  rownames(seTE.pooled) <- colnames(seTE.pooled) <- names.treat
  
  res <- list(studlab = studlab,
              treat1 = treat1, treat2 = treat2,
              TE = TE, seTE = seTE,
              seTE.orig = seTE.orig,
              TE.pooled = TE.pooled,
              seTE.pooled = seTE.pooled,
              tau = tau)
  
  res
}

