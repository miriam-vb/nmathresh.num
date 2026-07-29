#' Auxiliary functions and classes to support netmeta_slim
#'
#' Package: netmeta
#' Authors: Gerta Rücker <gerta.ruecker@@uniklinik-freiburg.de>, 
#'          Guido Schwarzer <guido.schwarzer@uniklinik-freiburg.de>
#' Source:  https://github.com/guido-s/netmeta
#' License: GPL (>= 2)
#' 
#' @import Matrix
#' @import MASS


prepare <- function(TE, seTE, treat1, treat2, studlab, tau = 0,
                    correlated = FALSE, func.inverse) {
  if (is.na(tau))
    tau <- 0
  
  data <- data.frame(studlab,
                     treat1, treat2,
                     treat1.pos = NA, treat2.pos = NA,
                     TE, seTE, weights = 1 / (seTE^2 + tau^2), correlated,
                     narms = NA, stringsAsFactors = FALSE)
  #
  # Ordering dataset
  #
  o <- order(data$studlab, data$treat1, data$treat2)
  data <- data[o, ]
  #
  # Adapt numbers to treatment IDs
  #
  names.treat <- sort(unique(c(data$treat1, data$treat2)))
  data$treat1.pos <- match(data$treat1, names.treat)
  data$treat2.pos <- match(data$treat2, names.treat)
  #
  data$order <- o
  
  sl <- unique(data$studlab)
  #
  # List with weight matrices
  #
  W.list <- vector("list", length(sl))
  names(W.list) <- sl
  #
  # List with covariance matrices
  #
  C.list <- vector("list", length(sl))
  names(C.list) <- sl
  #
  # Determining number of arms and adjusting weights of multi-arm studies
  #
  for (s in sl) {
    sel.s <- data$studlab == s
    correlated.s <- unique(data$correlated[sel.s]) 
    #
    if (length(correlated.s) != 1)
      stop("Different values for argument 'correlated' for study '", s, "'.",
           call. = FALSE)
    # Only treatment arms from multi-arm studies can be correlated
    if (correlated.s & sum(sel.s) == 1)
      correlated.s <- FALSE
    #
    res.s <- covar_study(1 / data$weights[sel.s], s, correlated.s, func.inverse)
    #
    W.list[[s]] <- res.s$W
    C.list[[s]] <- res.s$Cov
    #
    data$narms[sel.s] <- res.s$n
    data$weights[sel.s] <- diag(res.s$W)
  }
  #
  res <- list(W = bdiag(W.list), Cov = bdiag(C.list), data = data)
  #
  res
}


covar_study <- function(v, studlab, correlated, func.inverse) {
  m <- length(v)
  n <- (1 + sqrt(8 * m + 1)) / 2
  #
  if (correlated) {
    B <- createB(ncol = n)
    V <- diag(diag(t(B) %*% diag(v, nrow = m) %*% B)) - t(B) %*%
      diag(v, nrow = m) %*% B
    #
    Cov <- matrix(0, nrow = m, ncol = m)
    edges <- matrix(nrow = m, ncol = 2)
    #
    r <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        r <- r + 1
        edges[r, ] <- c(i, j)
      }
    }
    #
    for (p in 1:(m - 1)) {
      i <- edges[p, 1]
      j <- edges[p, 2]
      #
      for (q in (p+1):m) {
        k <- edges[q, 1]
        l <- edges[q, 2] 
        #
        Cov[p, q] <- 0.5 * (V[i, l] - V[i, k] + V[j, k] - V[j, l])
        Cov[q, p] <- 0.5 * (V[i, l] - V[i, k] + V[j, k] - V[j, l])
      }
    }
    #
    for (p in 1:m) {
      i <- edges[p, 1]
      j <- edges[p, 2]
      #
      Cov[p, p] <- V[i, j]
    }
    #
    if (qr(Cov)$rank == n - 1)
      W <- ginv(as.matrix(Cov))
    else {
      if (length(v) > 1)
        W <- diag(1 / v)
      else {
        Cov <- matrix(v)
        W <- 1 / Cov
      }
    }
  }
  else {
    if (length(v) > 1) {
      v <- multiarm(v, studlab, func.inverse)$v
      Cov <- diag(v)
      W <- diag(1 / v)
    }
    else {
      Cov <- matrix(v)
      W <- 1 / Cov
    }
  }
  #
  res <- list(v = v, n = n, m = m, Cov = as.matrix(Cov), W = W)
}


multiarm <- function(r, studlab, func.inverse) {
  ##
  ## Dimension of r and R
  ##
  m <- length(r) # Number of edges
  ##
  k <- (1 + sqrt(8 * m + 1)) / 2 # Number of vertices
  if (!(abs(k - round(k)) < .Machine$double.eps^0.5))
    stop("Wrong number of comparisons in multi-arm study.", call. = FALSE)
  ##
  ## Construct edge-vertex incidence matrix of complete graph of
  ## dimension k
  ##
  B <- createB(ncol = k)
  ##
  ## Distribute the edge variances on a symmetrical matrix R of
  ## dimension k x k
  ##
  R <- diag(diag(t(B) %*% diag(r, nrow = m) %*% B)) -
    t(B) %*% diag(r, nrow = m) %*% B
  ##
  ## Construct pseudoinverse Lt from given variance (resistance) matrix R
  ## using a theorem equivalent to Theorem 7 by Gutman & Xiao
  ## Lt <- -0.5 * (R - (R %*% J + J %*% R) / k + J %*%R %*% J / k^2)
  ##
  Lt <- -0.5 * t(B) %*% B %*% R %*% t(B) %*% B / k^2
  ##
  ## Compute Laplacian matrix L from Lt
  ##
  L <- do.call(func.inverse, list(X = Lt))
  ##
  ## Compute weight matrix W and variance matrix V from Laplacian L
  ## 
  W <- diag(diag(L)) - L
  ##
  ## Replace small negative weights with zeros
  ## (i.e., if an absolute negative weight contributes less than 0.1%)
  ##
  W[W < 0 & (abs(W) / sum(abs(W)[lower.tri(W)])) < 0.001] <- 0
  #
  V <- 1 / W
  ##
  ## Compute original variance vector v from V
  ##
  v <- rep(0, m)
  edge <- 0
  for (i in 1:(k - 1)) {
    for (j in (i + 1):k) {
      edge <- edge + 1
      v[edge] <- V[i, j]
    }
  }
  ##
  ## Result
  ##
  res <- list(k = k, r = r, R = R, Lt = Lt, L = L, W = W, V = V, v = v)
  res
}


createB <- function(pos1, pos2, ncol, aggr = FALSE) {
  
  
  if (!aggr) {
    if (missing(pos1) | missing(pos2)) {
      ##
      ## Create full edge-vertex incidence matrix
      ##
      nrow <- choose(ncol, 2)
      B <- matrix(0, nrow = nrow, ncol = ncol)
      ##
      i <- 0
      ##
      for (pos1.i in 1:(ncol - 1)) {
        for (pos2.i in (pos1.i + 1):ncol) {
          i <- i + 1
          B[i, pos1.i] <-  1
          B[i, pos2.i] <- -1
        }
      }
    }
    else {
      ##
      ## Create edge-vertex incidence matrix
      ##
      nrow <- length(pos1)
      ncol <- length(unique(c(pos1, pos2)))
      ##
      B <- matrix(0, nrow = nrow, ncol = ncol)
      ##
      for (i in 1:nrow) {
        B[i, pos1[i]] <-  1
        B[i, pos2[i]] <- -1
      }
    }
  }
  else {
    nrow <- 0
    ##
    ## Determine number of edges (no. of rows of B)
    ##
    for (i in 1:(ncol - 1)) {
      for (j in (i + 1):ncol) {
        ij.count <- 0
        ## Cycle through every possible edge ij
        ## Search pos1 and pos2 to see if at least one of these
        ## combinations is ij
        for (k in seq_along(pos1)) {
          if (pos1[k] == i & pos2[k] == j) {
            ij.count <- ij.count + 1
          }
          else {
            ij.count <- ij.count
          }
        }
        if (ij.count > 0)
          nrow <- nrow + 1
        else
          nrow <- nrow
      }
    }
    ##
    ## Create aggregate B matrix with dimensions e x n
    ##
    B <- matrix(0, nrow = nrow, ncol = ncol)
    ##
    r <- 0
    ## Cycle through each possible pairwise comparison ij
    for (i in 1:(ncol - 1)) {
      for (j in (i + 1):ncol) {
        ij.count <- 0
        for (k in 1:length(pos1)) {
          ## If there is an edge for that pairwise comparison ...
          if (pos1[k] == i & pos2[k] == j)
            ij.count <- ij.count + 1 # ...then ij.count is no longer = 0 ...
          else
            ij.count <- ij.count
        }
        if (ij.count > 0) {
          ## ...and we add this row to B
          r <- r + 1
          B[r, i] <-  1
          B[r, j] <- -1
        }
      }
    }
  }
  
  B
}

#
# Moore-Penrose Pseudoinverse of a Matrix
# 
invmat <- function(X) {
  n <- nrow(X)
  m <- ncol(X)
  ##
  if (n != m)
    stop("Argument 'X' must be a square matrix", call. = FALSE)
  ##
  J <- matrix(1, nrow = n, ncol = n)
  ##
  res <- solve(X - J / n) + J / n
  ##
  res
}


chkchar <- function(x, length = 0, name = NULL, nchar = NULL, single = FALSE) {
  if (!missing(single) && single)
    length <- 1
  if (is.null(name))
    name <- deparse(substitute(x))
  ##
  if (length && length(x) != length)
    stop("Argument '", name, "' must be a character vector of length ",
         length, ".",
         call. = FALSE)
  ##
  if (length == 1) {
    if (!is.null(nchar) && !(nchar(x) %in% nchar))
      if (length(nchar) == 1 && nchar == 1)
        stop("Argument '", name, "' must be a single character.",
             call. = FALSE)
    else
      stop("Argument '", name, "' must be a character string of length ",
           if (length(nchar) == 2)
             paste0(nchar, collapse = " or ")
           else
             paste0(nchar, collapse = ", "),
           ".",
           call. = FALSE)
  }
  ##
  if (!is.character(x))
    stop("Argument '", name, "' must be a character vector.")
  else {
    if (!is.null(nchar) & any(!(nchar(x) %in% nchar)))
      if (length(nchar) == 1 && nchar == 1)
        stop("Argument '", name, "' must be a vector of single characters.",
             call. = FALSE)
    else
      stop("Argument '", name, "' must be a character vector where ",
           "each element has ",
           if (length(nchar) == 2)
             paste0(nchar, collapse = " or ")
           else
             paste0(nchar, collapse = ", "),
           " characters.",
           call. = FALSE)
  }
}


chklogical <- function(x, name = NULL) {
  ##
  ## Check whether argument is logical
  ##
  if (is.null(name))
    name <- deparse(substitute(x))
  ##
  if (is.numeric(x))
    x <- as.logical(x)
  ##
  if (length(x) !=  1 || !is.logical(x) || is.na(x))
    stop("Argument '", name, "' must be a logical.", call. = FALSE)
  ##
  invisible(NULL)
}


chknumeric <- function(x, min, max, zero = FALSE, length = 0,
                       name = NULL, single = FALSE) {
  if (!missing(single) && single)
    length <- 1
  ##
  ## Check numeric variable
  ##
  if (is.null(name))
    name <- deparse(substitute(x))
  ##
  x <- x[!is.na(x)]
  if (length(x) == 0)
    return(NULL)
  ##
  if (!is.numeric(x))
    stop("Non-numeric value for argument '", name, "'.",
         call. = FALSE)
  ##
  if (length && length(x) != length)
    stop("Argument '", name, "' must be a numeric of length ", length, ".",
         call. = FALSE)
  ##
  if (!missing(min) & missing(max)) {
    if (zero & min == 0 & any(x <= min, na.rm = TRUE))
      stop("Argument '", name, "' must be positive.",
           call. = FALSE)
    else if (any(x < min, na.rm = TRUE))
      stop("Argument '", name, "' must be larger equal ",
           min, ".", call. = FALSE)
  }
  ##
  if (missing(min) & !missing(max)) {
    if (zero & max == 0 & any(x >= max, na.rm = TRUE))
      stop("Argument '", name, "' must be negative.",
           call. = FALSE)
    else if (any(x > max, na.rm = TRUE))
      stop("Argument '", name, "' must be smaller equal ",
           min, ".", call. = FALSE)
  }
  ##
  if ((!missing(min) & !missing(max)) &&
      (any(x < min, na.rm = TRUE) | any(x > max, na.rm = TRUE)))
    stop("Argument '", name, "' must be between ",
         min, " and ", max, ".", call. = FALSE)
  ##
  invisible(NULL)
}


chklength <- function(x, k.all, fun = "", text, name = NULL) {
  ##
  ## Check length of vector
  ##
  if (is.null(name))
    name <- deparse(substitute(x))
  ##
  if (length(x) != k.all) {
    funcs <- c("metabin", "metacont", "metacor",
               "metagen", "metainc", "metamean",
               "metaprop", "metarate",
               "funnel", "forest.meta")
    args <- c("event.e", "n.e", "cor",
              "TE", "event.e", "n",
              "event", "event",
              "TE", "TE")
    ##
    idx <- charmatch(fun, funcs, nomatch = NA)
    if (!is.na(idx))
      argname <- args[idx]
    else
      argname <- fun
    ##
    if (missing(text))
      stop("Arguments '", argname, "' and '", name,
           "' must have the same length.",
           call. = FALSE)
    else
      stop(text, call. = FALSE)
  }
  ##
  invisible(NULL)
}


chkclass <- function(x, class, name = NULL) {
  ##
  ## Check class of R object
  ##
  if (is.null(name))
    name <- deparse(substitute(x))
  ##
  n.class <- length(class)
  if (n.class == 1)
    text.class <- paste0('"', class, '"')
  else if (n.class == 2)
    text.class <- paste0('"', class, '"', collapse = " or ")
  else
    text.class <- paste0(paste0('"', class[-n.class], '"', collapse = ", "),
                         ', or ', '"', class[n.class], '"')
  ##
  if (!inherits(x, class))
    stop("Argument '", name,
         "' must be an object of class \"",
         text.class, "\".", call. = FALSE)
  ##
  invisible(NULL)
}


catch <- function(argname, matchcall, data, encl) {
  #
  # Catch value for argument
  #
  eval(matchcall[[match(argname, names(matchcall))]], data, enclos = encl)
}


replaceNULL <- function(x, replace = NA) {
  if (is.null(x))
    return(replace)
  x
}


rmSpace <- function(x, end = FALSE, pat = " ") {
  
  if (!end) {
    while (any(substring(x, 1, 1) == pat, na.rm = TRUE)) {
      sel <- substring(x, 1, 1) == pat
      x[sel] <- substring(x[sel], 2)
    }
  }
  else {
    last <- nchar(x)
    
    while (any(substring(x, last, last) == pat, na.rm = TRUE)) {
      sel <- substring(x, last, last) == pat
      x[sel] <- substring(x[sel], 1, last[sel] - 1)
      last <- nchar(x)
    }
  }
  
  x
}


setsep <- function(x, sep, type = "treatment label",
                   argname = deparse(substitute(sep)),
                   missing = TRUE) {
  labels <- sort(unique(x))
  #
  if (compmatch(labels, sep)) {
    if (!missing)
      warning("Argument '", argname, "': ",
              "separator '", sep, "' used in at least one ", type, ". ",
              "Trying to use predefined separators: ",
              "':', '-', '_', '/', '+', '.', '|', '*'.",
              call. = FALSE)
    #
    if (!compmatch(labels, ":"))
      sep <- ":"
    else if (!compmatch(labels, "-"))
      sep <- "-"
    else if (!compmatch(labels, "_"))
      sep <- "_"
    else if (!compmatch(labels, "/"))
      sep <- "/"
    else if (!compmatch(labels, "+"))
      sep <- "+"
    else if (!compmatch(labels, "."))
      sep <- "."
    else if (!compmatch(labels, "|"))
      sep <- "|"
    else if (!compmatch(labels, "*"))
      sep <- "*"
    else
      stop("All predefined separators (':', '-', '_', '/', '+', ",
           "'.', '|', '*') are used in at least one ", type, ". ",
           "Please specify a different character that should be ",
           "used as separator (argument '", argname, "').",
           call. = FALSE)
  }
  #
  sep
}


compmatch <- function(x, split) {
  
  if (split %in% c("+", ".", "&", "$", "#", "|", "*", "^"))
    split <- paste0("\\", split)
  
  res <- any(grepl(split, x))
  
  return(res)
}


is_wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  if (is.numeric(x))
    res <- abs(x - round(x)) < tol
  else
    res <- NA
  ##
  return(res)
}


is_zero <- function(x, n = 10) {
  return(abs(x) < n * .Machine$double.eps)
}


setref <- function(reference.group, levs, length = 1,
                   varname = "reference.group", error.text) {
  
  if (missing(error.text)) {
    text.start <- paste0("Argument '", varname, "'")
    text.within <- paste0("argument '", varname, "'")
  }
  else {
    text.start <- paste0(toupper(substring(error.text, 1, 1)),
                         substring(error.text, 2))
    text.within <- error.text
  }
  
  
  if (length && length(reference.group) != length)
    stop(text.start,
         if (length == 1)
           " must be a numeric or a character string"
         else
           paste(" must be a numeric of character vector of length", length),
         ".",
         call. = FALSE)
  ##
  if (is.numeric(reference.group)) {
    if (any(is.na(reference.group)))
      stop("Missing value not allowed in ", text.within, ".",
           call. = FALSE)
    if (!all(reference.group %in% seq_len(length(levs))))
      stop(paste0(text.start, " must ",
                  if (length == 1) "be any of the " else "contain ",
                  "integers from 1 to ",
                  length(levs), "."),
           call. = FALSE)
    res <- levs[reference.group]
  }
  else if (is.character(reference.group)) {
    if (any(is.na(reference.group)))
      stop("Missing value not allowed in ", text.within, ".",
           call. = FALSE)
    ##
    if (length(unique(levs)) == length(unique(tolower(levs))))
      idx <- charmatch(tolower(reference.group), tolower(levs), nomatch = NA)
    else {
      idx1 <- charmatch(reference.group, levs, nomatch = NA)
      idx2 <- charmatch(tolower(reference.group), tolower(levs), nomatch = NA)
      if (anyNA(idx1) & !anyNA(idx2))
        idx <- idx2
      else
        idx <- idx1
    }
    ##
    if (anyNA(idx) || any(idx == 0))
      stop("Admissible values for ", text.within, ":\n  ",
           paste(paste0("'", levs, "'"), collapse = " - "),
           "\n  (unmatched value", if (sum(is.na(idx)) > 1) "s",
           ": ",
           paste(paste0("'", reference.group[is.na(idx)], "'"),
                 collapse = " - "),
           ")",
           call. = FALSE)
    res <- levs[idx]
  }
  
  res
}


chkmultiarm <- function(TE, seTE, treat1, treat2, studlab, correlated,
                        tol.multiarm = 0.001,
                        tol.multiarm.se = NULL,
                        details = FALSE, debug = FALSE) {

  #
  # Ordering dataset (if necessary)
  #
  o <- order(studlab, treat1, treat2)
  #
  if (any(o != seq_along(studlab))) {
    TE <- TE[o]
    seTE <- seTE[o]
    treat1 <- treat1[o]
    treat2 <- treat2[o]
    studlab <- studlab[o]
    correlated <- correlated[o]
  }
  
  
  tabnarms <- table(studlab)
  sel.multi <- tabnarms > 1
  #
  if (any(sel.multi)) {
    #
    msgdetails <-
      paste0("  - For more details, re-run netmeta() with argument ",
             "details.chkmultiarm = TRUE.\n")
    #
    studlab.multi <- names(tabnarms)[sel.multi]
    #
    # Check duplicate and incomplete comparisons
    #
    dat.duplicate <- dat.incomplete <-
      data.frame(studlab = "", treat1 = "", treat2 = "",
                 stringsAsFactors = FALSE)
    #
    incomplete <- rep_len(NA, sum(sel.multi))
    duplicate <- rep_len(NA, sum(sel.multi))
    #
    s.idx <- 0
    
    #
    # (1) Check for multi-arm studies with incomplete or duplicate treatment
    #     comparisons
    #
    for (s in studlab.multi) {
      s.idx <- s.idx + 1
      sel <- studlab == s
      #
      TE.s <- TE[sel]
      studlab.s <- studlab[sel]
      treat1.s <- treat1[sel]
      treat2.s <- treat2[sel]
      treats.s <- unique(c(treat1.s, treat2.s))
      #
      n <- (1 + sqrt(8 * length(TE.s) + 1)) / 2
      #
      incomplete[s.idx] <- length(treats.s) != n
      duplicate[s.idx] <- any(table(interaction(treat1.s, treat2.s)) > 1)
      #
      if (incomplete[s.idx])
        dat.incomplete <- rbind(dat.incomplete,
                                data.frame(studlab = studlab.s,
                                           treat1 = treat1.s,
                                           treat2 = treat2.s,
                                           stringsAsFactors = FALSE))
      #
      if (duplicate[s.idx])
        dat.duplicate <- rbind(dat.duplicate,
                               data.frame(studlab = studlab.s,
                                          treat1 = treat1.s,
                                          treat2 = treat2.s,
                                          stringsAsFactors = FALSE))
    }
    #
    if (details & any(incomplete)) {
      dat.incomplete <- dat.incomplete[-1, ]
      cat("\nMulti-arm studies with incomplete treatment comparisons:\n\n")
      prmatrix(dat.incomplete, quote = FALSE, right = TRUE,
               rowlab = rep("", dim(dat.incomplete)[1]))
      cat("\n")
    }
    #
    if (details & any(duplicate)) {
      dat.duplicate <- dat.duplicate[-1, ]
      cat("\nMulti-arm studies with duplicate treatment comparisons:\n\n")
      prmatrix(dat.duplicate, quote = FALSE, right = TRUE,
               rowlab = rep("", dim(dat.duplicate)[1]))
      cat("\n")
    }
    #
    if (any(incomplete)) {
      studlabs <- unique(dat.incomplete$studlab[-1])
      #
      if (length(studlabs) == 1)
        errmsg.incomplete <-
          paste0("  - Study '", studlabs,
                 "' has an incomplete set of comparisons.\n")
      else
        errmsg.incomplete <-
          paste0("  - Studies with incomplete set of comparisons: ",
                 paste(paste0("'", studlabs, "'"), collapse = ", "),
                 "\n")
    }
    else
      errmsg.incomplete <- ""
    #
    if (any(duplicate)) {
      studlabs <- unique(dat.duplicate$studlab[-1])
      #
      if (length(studlabs) == 1)
        errmsg.duplicate <-
          paste0("  - Duplicate comparison in study '", studlabs, "'.\n")
      else
        errmsg.duplicate <-
          paste0("  - Studies with duplicate comparisons: ",
                 paste(paste0("'", studlabs, "'"), collapse = ", "),
                 "\n")
    }
    else
      errmsg.duplicate <- ""
    #
    if (any(incomplete) | any(duplicate))
      stop("Problem",
           if ((sum(incomplete) + sum(duplicate)) > 1) "s",
           " in multi-arm studies!\n",
           errmsg.incomplete, errmsg.duplicate,
           if (!details) msgdetails, call. = FALSE)
    
    #
    # (2) Check for (i) consistency of TE and varTE or (2) negative or zero
    #     treatment arm variance (zero variance only results in a warning)
    #
    dat.TE <- data.frame(studlab = "", treat1 = "", treat2 = "",
                         TE = NA, resid = NA,
                         stringsAsFactors = FALSE)
    #
    dat.varTE <- data.frame(studlab = "", treat1 = "", treat2 = "",
                            varTE = NA, resid.var = NA,
                            seTE = NA, resid.se = NA,
                            stringsAsFactors = FALSE)
    #
    dat.negative <- data.frame(studlab = "", treat = "", var.treat = NA,
                               stringsAsFactors = FALSE)
    #
    dat.zero <- data.frame(studlab = "", treat = "", var.treat = NA,
                           stringsAsFactors = FALSE)
    #
    inconsistent.TE <- inconsistent.varTE <-
      zero.sigma2 <- negative.sigma2 <- rep_len(NA, sum(sel.multi))
    #
    s.idx <- 0
    #
    for (s in studlab.multi) {
      s.idx <- s.idx + 1
      sel <- studlab == s
      #
      TE.s <- TE[sel]
      seTE.s <- seTE[sel]
      varTE.s <- seTE.s^2
      studlab.s <- studlab[sel]
      treat1.s <- treat1[sel]
      treat2.s <- treat2[sel]
      #
      correlated.s <- unique(correlated[sel])
      #
      n <- (1 + sqrt(8 * length(TE.s) + 1)) / 2
      #
      treats.s <- unique(c(treat1.s, treat2.s))
      studlab.s.arms <- rep_len(s, n)
      #
      # Create full edge-vertex incidence matrix
      #
      B <- createB(ncol = n)
      #
      # Check treatment estimates
      #
      TE.diff <- TE.s - B %*% as.vector(ginv(B) %*% TE.s)
      if (debug) {
        cat("*** TE.diff = TE - B %*% as.vector(ginv(B) %*% TE) ***\n")
        print(data.frame(TE = TE.s,
                         TE.calc =  B %*% as.vector(ginv(B) %*% TE.s),
                         TE.diff))
      }
      #
      inconsistent.TE[s.idx] <- any(abs(TE.diff) > tol.multiarm)
      #
      if (any(abs(TE.diff) > tol.multiarm))
        dat.TE <- rbind(dat.TE,
                        data.frame(studlab = studlab.s,
                                   treat1 = treat1.s,
                                   treat2 = treat2.s,
                                   TE = round(TE.s, 8),
                                   resid = round(TE.diff, 8),
                                   stringsAsFactors = FALSE))
      #
      # Check standard errors
      #
      if (correlated.s) {
        inconsistent.varTE[s.idx] <- 0
        negative.sigma2[s.idx] <- 0
        zero.sigma2[s.idx] <- 0
      }
      else {
        A <- abs(B)
        #
        sigma2 <- as.vector(ginv(A) %*% varTE.s)
        #
        varTE.diff <- varTE.s - A %*% sigma2
        #
        if (debug) {
          cat("*** varTE.diff = varTE - A %*% sigma2 ***\n")
          cat("*** with sigma2 = as.vector(ginv(A) %*% varTE) ***\n")
          print(data.frame(varTE = varTE.s,
                           varTE.calc = A %*% sigma2,
                           varTE.diff))
          print(data.frame(sigma2))
        }
        #
        if (!is.null(tol.multiarm.se))
          inconsistent.varTE[s.idx] <- any(abs(varTE.diff) > tol.multiarm.se^2)
        else
          inconsistent.varTE <- rep_len(FALSE, length(inconsistent.varTE))
        #
        is.negative <- sigma2 < 0
        negative.sigma2[s.idx] <- any(is.negative)
        zero.sigma2[s.idx] <- any(is_zero(sigma2[!is.negative]))
        #
        if (inconsistent.varTE[s.idx])
          dat.varTE <- rbind(dat.varTE,
                             data.frame(studlab = studlab.s,
                                        treat1 = treat1.s,
                                        treat2 = treat2.s,
                                        varTE = round(varTE.s, 8),
                                        resid.var = round(varTE.diff, 8),
                                        seTE = round(sqrt(varTE.s), 8),
                                        resid.se = sign(varTE.diff) *
                                          round(sqrt(abs(varTE.diff)), 8),
                                        stringsAsFactors = FALSE))
        #
        if (negative.sigma2[s.idx])
          dat.negative <- rbind(dat.negative,
                                data.frame(studlab = studlab.s.arms,
                                           treat = treats.s,
                                           var.treat = sigma2,
                                           stringsAsFactors = FALSE))
        #
        if (zero.sigma2[s.idx])
          dat.zero <- rbind(dat.zero,
                            data.frame(studlab = studlab.s.arms,
                                       treat = treats.s,
                                       var.treat = round(sigma2, 8),
                                       stringsAsFactors = FALSE))
      }
    }
    #
    iTE <- sum(inconsistent.TE)
    ivarTE <- sum(inconsistent.varTE)
    inconsistent <- iTE > 0 | ivarTE > 0
    #
    inegative.sigma2 <- sum(negative.sigma2)
    izero.sigma2 <- sum(zero.sigma2)
    zero <- izero.sigma2 > 0
    negative <- inegative.sigma2 > 0
    #
    studlab.inconsistent.TE <- character()
    studlab.inconsistent.varTE <- character()
    studlab.zero.sigma2 <- character()
    studlab.negative.sigma2 <- character()
    if (iTE > 0)
      studlab.inconsistent.TE <- studlab.multi[inconsistent.TE]
    if (ivarTE > 0)
      studlab.inconsistent.varTE <- studlab.multi[inconsistent.varTE]
    if (zero)
      studlab.zero.sigma2 <- studlab.multi[zero.sigma2]
    if (negative)
      studlab.negative.sigma2 <- studlab.multi[negative.sigma2]
    #
    studlab.inconsistent <-
      unique(c(studlab.inconsistent.TE, studlab.inconsistent.varTE,
               studlab.zero.sigma2, studlab.negative.sigma2))
    
    #
    # Print information on deviations from consistency assumption in
    # multi-arm studies
    #
    if (details & (inconsistent | zero | negative)) {
      if (length(dat.TE$studlab) > 1) {
        dat.TE <- dat.TE[-1, ]
        dat.TE$TE <- format(dat.TE$TE)
        dat.TE$resid <- format(dat.TE$resid)
        cat("\nMulti-arm studies with inconsistent treatment effects:\n\n")
        prmatrix(dat.TE, quote = FALSE, right = TRUE,
                 rowlab = rep("", dim(dat.TE)[1]))
        cat("\n")
      }
      if (length(dat.varTE$studlab) > 1 & ivarTE > 0) {
        dat.varTE <- dat.varTE[-1, ]
        dat.varTE$varTE <- format(dat.varTE$varTE)
        dat.varTE$resid.var <- format(dat.varTE$resid.var)
        dat.varTE$seTE <- format(dat.varTE$seTE)
        dat.varTE$resid.se <- format(dat.varTE$resid.se)
        cat("\nMulti-arm studies with inconsistent",
            "variances / standard errors:\n\n")
        prmatrix(dat.varTE, quote = FALSE, right = TRUE,
                 rowlab = rep("", dim(dat.varTE)[1]))
        cat("\n")
      }
      #
      # Negative variance
      #
      if (length(dat.negative$studlab) > 1 & negative) {
        dat.negative <- dat.negative[-1, ]
        cat("\nMulti-arm studies with negative treatment arm variance:\n\n")
        prmatrix(dat.negative, quote = FALSE, right = TRUE,
                 rowlab = rep("", dim(dat.negative)[1]))
        cat("\n")
      }
      #
      # Zero variance
      #
      if (length(dat.zero$studlab) > 1 & zero) {
        dat.zero <- dat.zero[-1, ]
        cat("\nMulti-arm studies with zero treatment arm variance:\n\n")
        prmatrix(dat.zero, quote = FALSE, right = TRUE,
                 rowlab = rep("", dim(dat.zero)[1]))
        cat("\n")
        #
        warning(
          paste0("Note, a zero treatment arm variance has been calculated ",
                 "for the following multi-arm ",
                 if (izero.sigma2 == 1) "study" else "studies",
                 ": ",
                 paste(paste0("'", studlab.zero.sigma2, "'"), collapse = ", ")),
          call. = FALSE)
      }
      #
      cat("Legend:\n")
      if (inconsistent)
        cat(" resid - residual deviation (observed minus expected)\n")
      if (iTE > 0)
        cat(" TE", if (ivarTE > 0) "   ", " - treatment estimate\n", sep = "")
      if (ivarTE > 0) {
        cat(" varTE - variance of treatment estimate\n")
        cat(" seTE  - standard error of treatment estimate\n")
      }
      if (negative | zero)
        cat(" var.treat - treatment arm variance\n")
      cat("\n")
    }
    
    
    #
    # Generate error message
    #
    if (inconsistent | negative) {
      #
      if (iTE > 0)
        msgTE <-
          paste0("  ",
                 if (iTE == 1) "- Study " else "- Studies ",
                 "with inconsistent treatment estimates: ",
                 paste(paste0("'", studlab.inconsistent.TE, "'"),
                       collapse = ", "),
                 "\n")
      else
        msgTE <- ""
      #
      if (ivarTE > 0)
        msgvarTE <-
          paste0("  ",
                 if (ivarTE == 1) "- Study " else "- Studies ",
                 "with inconsistent variances: ",
                 paste(paste0("'", studlab.inconsistent.varTE, "'"),
                       collapse = ", "),
                 "\n")
      else
        msgvarTE <- ""
      #
      if (negative)
        msgsigma2 <-
          paste0("  ",
                 if (inegative.sigma2 == 1) "- Study " else "- Studies ",
                 "with negative treatment arm variance: ",
                 paste(paste0("'", studlab.negative.sigma2, "'"),
                       collapse = ", "),
                 "\n",
                 "    Potential solutions:\n",
                 "    1. Use argument 'func.inverse' to specify a different ",
                 "function for matrix inversion;\n",
                 "    2. Use argument 'correlated' to identify studies with ",
                 "correlated treatment arms, e.g., due to body-split design;\n",
                 "    3. Fix data errors.\n")
      else
        msgsigma2 <- ""
      #
      errmsg <-
        paste0("Problem",
               if ((iTE + ivarTE + inegative.sigma2 + izero.sigma2) > 1) "s",
               " in multi-arm studies!\n",
               msgTE, msgvarTE, msgsigma2,
               "  - Please check original data used as input to netmeta().\n",
               if (!details) msgdetails,
               if (inconsistent)
                 paste0("  - Argument",
                        if (iTE & ivarTE) "s",
                        if (iTE) " 'tol.multiarm'",
                        if (iTE & ivarTE) " and",
                        if (ivarTE) " 'tol.multiarm.se'",
                        " in netmeta() can be used to",
                        if (iTE & ivarTE) "\n   ",
                        " relax consistency",
                        if (!(iTE & ivarTE)) "\n   ",
                        " assumption for multi-arm studies (if appropriate)."))
      #
      stop(errmsg, call. = FALSE)
    }
  }
  
  invisible(NULL)
}

#' @rdname netconnection
#' @method netconnection default
#' @export

netconnection.default <- function(data = NULL, treat1, treat2, studlab = NULL,
                                  subset = NULL,
                                  sep.trts = ":",
                                  nchar.trts = 666,
                                  title = "", details.disconnected = FALSE,
                                  warn = FALSE, ...) {
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  
  nulldata <- is.null(data)
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  ##
  ## Catch treat1
  ##
  if (!nulldata & !is.data.frame(data)) {
    if (!missing(treat1) & !missing(treat2) & !missing(studlab))
      stop("Argument 'data' must be a data frame.")
    else if (!missing(treat1) & !missing(treat2)) {
      treat1 <- catch("data", mc, sfsp, sfsp)
      treat2 <- catch("treat1", mc, sfsp, sfsp)
      studlab <- catch("treat2", mc, sfsp, sfsp)
    }
    else if (!missing(treat1)) {
      treat1 <- catch("data", mc, sfsp, sfsp)
      treat2 <- catch("treat1", mc, sfsp, sfsp)
    }
  }
  else {
    ##
    if (nulldata)
      data <- sfsp
    ##
    treat1 <- catch("treat1", mc, data, sfsp)
    treat2 <- catch("treat2", mc, data, sfsp)
    studlab <- catch("studlab", mc, data, sfsp)
    subset <- catch("subset", mc, data, sfsp)
  }
  ##
  if (length(studlab) != 0)
    studlab <- as.character(studlab)
  else {
    if (warn)
      warning("No information given for argument 'studlab'. ",
              "Assuming that comparisons are from independent studies.")
    studlab <- as.character(seq_along(treat1))
  }
  ##
  chknumeric(nchar.trts, min = 1, length = 1)
  ##
  chklogical(details.disconnected)
  chklogical(warn)
  
  
  ##
  ##
  ## (2) Check length of essential variables
  ##
  ##
  
  fun <- "netconnection"
  ##
  k.All <- length(treat1)
  ##
  missing.subset <- is.null(subset)
  ##
  chklength(treat2, k.All, fun,
            text = paste0("Arguments 'treat1' and 'treat2' ",
                          "must have the same length."))
  chklength(studlab, k.All, fun,
            text = paste0("Arguments 'treat1' and 'studlab' ",
                          "must have the same length."))
  ##
  if (is.factor(treat1))
    treat1 <- as.character(treat1)
  if (is.factor(treat2))
    treat2 <- as.character(treat2)
  
  
  ##
  ##
  ## (3) Use subset for analysis
  ##
  ##
  
  if (!missing.subset) {
    if ((is.logical(subset) & (sum(subset) > k.All)) ||
        (length(subset) > k.All))
      stop("Length of subset is larger than number of studies.")
    ##
    treat1 <- treat1[subset]
    treat2 <- treat2[subset]
    studlab <- studlab[subset]
  }
  
  
  ##
  ##
  ## (4) Additional checks
  ##
  ##
  
  if (any(treat1 == treat2))
    stop("Treatments must be different (arguments 'treat1' and 'treat2').")
  ##
  ## Check for correct number of comparisons
  ##
  tabnarms <- table(studlab)
  sel.narms <- !is_wholenumber((1 + sqrt(8 * tabnarms + 1)) / 2)
  ##
  if (sum(sel.narms) == 1)
    stop("Study '", names(tabnarms)[sel.narms],
         "' has a wrong number of comparisons.",
         "\n  Please provide data for all treatment comparisons ",
         "(two-arm: 1; three-arm: 3; four-arm: 6, ...).")
  if (sum(sel.narms) > 1)
    stop("The following studies have a wrong number of comparisons: ",
         paste(paste0("'", names(tabnarms)[sel.narms], "'"),
               collapse = ", "),
         "\n  Please provide data for all treatment comparisons ",
         "(two-arm: 1; three-arm: 3; four-arm: 6, ...).")
  ##
  labels <- sort(unique(c(treat1, treat2)))
  ##
  if (compmatch(labels, sep.trts)) {
    if (!missing(sep.trts))
      warning("Separator '", sep.trts,
              "' used in at least one treatment label. ",
              "Try to use predefined separators: ",
              "':', '-', '_', '/', '+', '.', '|', '*'.",
              call. = FALSE)
    ##
    if (!compmatch(labels, ":"))
      sep.trts <- ":"
    else if (!compmatch(labels, "-"))
      sep.trts <- "-"
    else if (!compmatch(labels, "_"))
      sep.trts <- "_"
    else if (!compmatch(labels, "/"))
      sep.trts <- "/"
    else if (!compmatch(labels, "+"))
      sep.trts <- "+"
    else if (!compmatch(labels, "."))
      sep.trts <- "-"
    else if (!compmatch(labels, "|"))
      sep.trts <- "|"
    else if (!compmatch(labels, "*"))
      sep.trts <- "*"
    else
      stop("All predefined separators (':', '-', '_', '/', '+', ",
           "'.', '|', '*') are used in at least one treatment label.",
           "\n   Please specify a different character that should be ",
           "used as separator (argument 'sep.trts').",
           call. = FALSE)
  }
  
  
  ##
  ##
  ## (5) Determine (sub)network(s)
  ##
  ##
  
  treats <- as.factor(c(as.character(treat1), as.character(treat2)))
  trts <- levels(treats)
  ##
  n <- length(trts)   # Number of treatments
  m <- length(treat1) # Number of comparisons
  ##
  ## Edge-vertex incidence matrix
  ##
  treat1.pos <- treats[1:m]
  treat2.pos <- treats[(m + 1):(2 * m)]
  B <- createB(treat1.pos, treat2.pos, ncol = n)
  ##
  rownames(B) <- studlab
  colnames(B) <- trts
  ##
  L.mult <- t(B) %*% B             # Laplacian matrix with multiplicity
  A <- diag(diag(L.mult)) - L.mult # Adjacency matrix
  D <- netdistance(A)              # Distance matrix
  L <- diag(rowSums(A)) - A        # Laplacian matrix without multiplicity
  ##
  n.subsets <- as.integer(table(round(eigen(L)$values, 10) == 0)[2])
  ##
  ## Block diagonal matrix in case of sub-networks
  ##
  maxdist <- nrow(D)
  D2 <- D
  D2[is.infinite(D2)] <- maxdist
  o <- hclust(dist(D2))$order
  ##
  D <- D[o, o]
  A <- A[o, o]
  L <- L[o, o]
  ##
  A.loop <- A
  D.loop <- D
  id.treats <- character(0)
  id.subnets <- numeric(0)
  more.subnets <- TRUE
  subnet.i <- 0
  dat.subnet <- data.frame()
  ##
  while (more.subnets) {
    subnet.i <- subnet.i + 1
    n.i <- seq_len(nrow(D.loop))
    #
    next.subnet <- min(max(n.i) + 1, n.i[is.infinite(D.loop[, 1])])
    sel.i <- seq_len(next.subnet - 1)
    #
    id.treats <- c(id.treats, rownames(D.loop)[sel.i])
    id.subnets <- c(id.subnets, rep(subnet.i, length(sel.i)))
    ##
    A.i <- A.loop[sel.i, sel.i]
    A.i <- A.i[order(rownames(A.i)), order(colnames(A.i))]
    ##
    for (row.i in 1:(ncol(A.i) - 1)) {
      for (col.i in 2:ncol(A.i)) {
        if (col.i > row.i) {
          if (A.i[row.i, col.i] > 0) {
            trt1.i <- rownames(A.i)[row.i]
            trt2.i <- rownames(A.i)[col.i]
            #
            dat.subnet <-
              rbind(
                dat.subnet,
                data.frame(subnet = subnet.i,
                           comparison = paste(trt1.i, trt2.i, sep = sep.trts),
                           treat1 = trt1.i, treat2 = trt2.i))
          }
        }
      }
    }
    #
    A.loop <- A.loop[-sel.i, -sel.i]
    D.loop <- D.loop[-sel.i, -sel.i]
    #
    more.subnets <- nrow(D.loop) > 0
  }
  #
  # Order matrices within subnets by treatments
  # (data frame 'dat.subnet' is already sorted by subnetwork and treatments)
  #
  seq <- NULL
  for (i in unique(dat.subnet$subnet)) {
    dat.i <- dat.subnet[dat.subnet$subnet == i, ]
    seq <- c(seq, unique(c(dat.i$treat1, dat.i$treat2)))
  }
  #
  A <- A[seq, seq]
  D <- D[seq, seq]
  L <- L[seq, seq]
  #
  # Add subnetwork number to comparisons
  #
  subnet <- rep(NA, length(treat1))
  #
  for (i in seq_along(id.treats))
    subnet[treat1 == id.treats[i]] <- id.subnets[i]
  ##
  comparisons <- dat.subnet$comparison
  subnet.comparisons <- dat.subnet$subnet
  #
  designs <- designs(treat1, treat2, studlab)
  
  res <- list(treat1 = treat1,
              treat2 = treat2,
              studlab = studlab,
              design = designs$design,
              subnet = subnet,
              #
              k = length(unique(studlab)),
              m = m,
              n = n,
              n.subnets = n.subsets,
              d = length(unique(designs$design)),
              #
              seq = seq,
              #
              D.matrix = D,
              A.matrix = A,
              L.matrix = L,
              ##
              designs = unique(sort(designs$design)),
              comparisons = comparisons,
              subnet.comparisons = subnet.comparisons,
              ##
              nchar.trts = nchar.trts,
              ##
              title = title,
              ##
              details.disconnected = details.disconnected,
              ##
              warn = warn,
              call = match.call(),
              version = packageDescription("netmeta")$Version
  )
  
  class(res) <- "netconnection"
  
  res
}


#' @rdname netconnection
#' @method netconnection pairwise
#' @export

netconnection.pairwise <- function(data,
                                   treat1, treat2, studlab = NULL,
                                   subset = NULL,
                                   sep.trts = ":",
                                   drop.NA = TRUE,
                                   nchar.trts = 666,
                                   title = "", details.disconnected = FALSE,
                                   warn = FALSE,
                                   ...) {
  
  #
  #
  # (1) Check arguments
  #
  #
  
  chkclass(data, "pairwise")
  #
  chklogical(drop.NA)
  chklogical(warn)
  #
  # Arguments 'treat1', 'treat2' and 'studlab' ignored
  #
  if (warn) {
    if (!missing(treat1))
      warning("Argument 'treat1' ignored as argument 'data' is an ",
              "object created with pairwise().",
              call. = FALSE)
    #
    if (!missing(treat2))
      warning("Argument 'treat2' ignored as argument 'data' is an ",
              "object created with pairwise().",
              call. = FALSE)
    #
    if (!missing(studlab))
      warning("Argument 'studlab' ignored as argument 'data' is an ",
              "object created with pairwise().",
              call. = FALSE)
  }
  #
  treat1 <- data$treat1
  treat2 <- data$treat2
  studlab <- data$studlab
  #
  if (is.factor(treat1))
    treat1 <- as.character(treat1)
  if (is.factor(treat2))
    treat2 <- as.character(treat2)
  #
  missing.subset <- missing(subset)
  #
  if (!missing.subset) {
    sfsp <- sys.frame(sys.parent())
    mc <- match.call()
    subset <- catch("subset", mc, data, sfsp)
    #
    k.All <- length(treat1)
    if  (length(subset) > k.All)
      stop("Length of subset is larger than number of studies.")
    #
    #if ((is.logical(subset) & (sum(subset) > k.All)) ||
    #    (length(subset) > k.All))
    #  stop("Length of subset is larger than number of studies.")
    #
    if (is.numeric(subset)) {
      if (any(is.na(subset)))
        stop("No missing values allowed in argument 'subset'.")
      if (length(subset) != length(unique(subset)))
        stop("Duplicate values in argument 'subset'.")
      if (any(subset > k.All | subset <= 0))
        stop("Numerical values in argument 'subset' must be between 1 and ",
             k.All, ".")
      #
      subset1 <- rep_len(FALSE, k.All)
      subset1[subset] <- TRUE
      subset <- subset1
    }
  }
  else
    subset <- rep_len(TRUE, length(treat1))
  #
  if (drop.NA) {
    if (is.null(attr(data, "varnames"))) {
      TE <- data$TE
      seTE <- data$seTE
    }
    else {
      TE <- data[[attr(data, "varnames")[1]]]
      seTE <- data[[attr(data, "varnames")[2]]]
    }
    #
    subset <- subset & (!is.na(TE) & !(is.na(seTE) | seTE == 0))
  }
  #
  treat1 <- treat1[subset]
  treat2 <- treat2[subset]
  studlab <- studlab[subset]
  #
  chknumeric(nchar.trts, min = 1, length = 1)
  #
  chklogical(details.disconnected)
  
  
  #
  #
  # (2) Call netconnection.default()
  #
  #
  
  res <- netconnection(treat1 = treat1, treat2 = treat2,
                       studlab = studlab,
                       sep.trts = sep.trts, nchar.trts = nchar.trts,
                       title = title,
                       details.disconnected = details.disconnected,
                       warn = warn,
                       ...)
  #
  res
}


#' @rdname netconnection
#' @method netconnection netmeta
#' @export

netconnection.netmeta <- function(data,
                                  sep.trts = data$sep.trts,
                                  nchar.trts = data$nchar.trts,
                                  title = data$title,
                                  details.disconnected = FALSE,
                                  warn = FALSE, ...) {
  
  chkclass(data, "netmeta")
  
  res <- netconnection(treat1 = data$treat1, treat2 = data$treat2,
                       studlab = data$studlab,
                       sep.trts = sep.trts, nchar.trts = nchar.trts,
                       title = title)
  #
  res
}


#' @rdname netconnection
#' @method netconnection netcomb
#' @export

netconnection.netcomb <- function(data,
                                  sep.trts = data$sep.trts,
                                  nchar.trts,
                                  title = data$title,
                                  details.disconnected = FALSE,
                                  warn = FALSE, ...) {
  
  chkclass(data, "netcomb")
  
  if (inherits(data, "discomb")) {
    if (missing(nchar.trts))
      nchar.trts <- data$nchar.comps
    #
    return(netconnection(treat1 = data$treat1, treat2 = data$treat2,
                         studlab = data$studlab,
                         sep.trts = sep.trts, nchar.trts = nchar.trts,
                         title = title))
  }
  else {
    if (missing(nchar.trts))
      nchar.trts <- data$nchar.trts
    #
    return(netconnection(treat1 = data$treat1, treat2 = data$treat2,
                         studlab = data$studlab,
                         sep.trts = sep.trts, nchar.trts = nchar.trts,
                         title = title))
  }
}


#' @rdname netconnection
#' @method print netconnection
#' @export

print.netconnection <- function(x,
                                digits = max(4, .Options$digits - 3),
                                nchar.trts = x$nchar.trts,
                                details = FALSE,
                                details.disconnected = x$details.disconnected,
                                ...) {
  
  chkclass(x, "netconnection")
  ##
  if (is.null(nchar.trts))
    nchar.trts <- 666
  
  
  chknumeric(digits, length = 1)
  chknumeric(nchar.trts, min = 1, length = 1)
  chklogical(details)
  details.disconnected <- replaceNULL(details.disconnected, FALSE)
  chklogical(details.disconnected)
  
  
  matitle(x)
  ##
  cat("Number of studies: k = ", x$k, "\n", sep = "")
  cat("Number of pairwise comparisons: m = ", x$m, "\n", sep = "")
  cat("Number of treatments: n = ", x$n, "\n", sep = "")
  if (!is.null(x$d))
    cat("Number of designs: d = ", x$d, "\n", sep = "")
  ##
  cat("Number of networks: ", x$n.subnets, "\n", sep = "")
  ##
  if (x$n.subnets > 1) {
    f <- function(x) length(unique(x))
    d <- as.data.frame(x)
    k.subset <- tapply(d$studlab, d$subnet, f)
    ##
    m <- as.matrix(
      data.frame(subnetwork = names(k.subset),
                 k = as.vector(k.subset),
                 m = as.vector(tapply(d$studlab, d$subnet, length)),
                 n = as.vector(tapply(c(d$treat1, d$treat2),
                                      c(d$subnet, d$subnet), f)))
    )
    rownames(m) <- rep("", nrow(m))
    ##
    cat("\nDetails on subnetworks: \n")
    prmatrix(m, quote = FALSE, right = TRUE)
    ##
    if (details.disconnected) {
      cat("\n")
      for (i in seq_len(x$n.subnets)) {
        d.i <- subset(d, d$subnet == i)
        cat("Subnetwork ", i, ":\n", sep = "")
        print(sort(unique(c(d.i$treat1, d.i$treat2))))
      }
    }
  }
  
  
  if (details) {
    cat("\nDistance matrix:\n")
    
    D <- round(x$D.matrix, digits = digits)
    D[is.infinite(D)] <- "."
    ##
    if (x$n.subnets == 1)
      diag(D) <- "."
    ##
    rownames(D) <- treats(rownames(D), nchar.trts)
    colnames(D) <- treats(colnames(D), nchar.trts)
    ##
    prmatrix(D, quote = FALSE, right = TRUE)
    
    diff.rownames <- rownames(x$D.matrix) != rownames(D)
    if (any(diff.rownames)) {
      abbr <- rownames(D)
      full <- rownames(x$D.matrix)
      ##
      tmat <- data.frame(abbr, full)
      names(tmat) <- c("Abbreviation", "Treatment name")
      tmat <- tmat[diff.rownames, ]
      tmat <- tmat[order(tmat$Abbreviation), ]
      ##
      cat("\nLegend:\n")
      prmatrix(tmat, quote = FALSE, right = TRUE,
               rowlab = rep("", length(abbr)))
    }
  }
  
  invisible(NULL)
}


#' @rdname netconnection
#' @export

netconnection <- function(data, ...)
  UseMethod("netconnection")


#' @rdname netdistance
#' @method netdistance default
#' @export

netdistance.default <- function(x, ...) {
  
  # Calculate distance matrix D of adjacency matrix A based on
  # distance algorithm by Mueller et al. (1987) using triangle
  # inequality
  
  chkclass(x, "matrix")
  #
  A <- x
  
  # Starting value for D is sign(A), with 0 replaced by Inf
  #
  n <- nrow(A)
  D <- sign(A)
  #
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (D[i, j] == 0) {
        D[i, j] <- Inf
        D[j, i] <- Inf
      }
    }
  }
  #
  for (d in 1:(n - 1)) {
    for (i in 1:n) {
      for (j in 1:n) {
        if (D[i, j] == d) {
          for (k in 1:n) {
            akj <- D[k, i] + d # = D[k, i] + D[i, j]
            D[k, j] <- min(D[k, j], akj)
          }
        }
      }
    }
  }
  #
  maxdist <- nrow(D)
  D2 <- D
  D2[is.infinite(D2)] <- maxdist
  attr(D, "order") <- hclust(dist(D2))$order
  
  class(D) <- c("netdistance", class(D))
  #
  D
}


#' @rdname netdistance
#' @method netdistance netmeta
#' @export

netdistance.netmeta <- function(x, sort = gs("sort.distance"), ...) {
  
  chkclass(x, "netmeta")
  chklogical(sort)
  
  A <- x$A.matrix
  #
  if (sort) {
    seq <- netconnection(x$treat1, x$treat2)$seq
    A <- A[seq, seq]
  }
  
  res <- netdistance(A)
  #
  if (sort)
    attr(res, "order") <- NULL
  #
  res
}


#' @rdname netdistance
#' @method netdistance netcomb
#' @export

netdistance.netcomb <- function(x, sort = gs("sort.distance"), ...) {
  
  chkclass(x, "netcomb")
  chklogical(sort)
  
  if (inherits(x, "discomb")) {
    A <- x$A.matrix
    #
    if (sort)
      seq <- netconnection(x$treat1, x$treat2)$seq
  }
  else {
    A <- x$x$A.matrix
    #
    if (sort)
      seq <- netconnection(x$x$treat1, x$x$treat2)$seq
  }
  #
  if (sort)
    A <- A[seq, seq]
  
  res <- netdistance(A)
  #
  if (sort)
    attr(res, "order") <- NULL
  #
  res
}


#' @rdname netdistance
#' @method netdistance netconnection
#' @export

netdistance.netconnection <- function(x, ...) {
  
  chkclass(x, "netconnection")
  
  netdistance(x$A.matrix)
}


#' @rdname netdistance
#' @method print netdistance
#' @export

print.netdistance <- function(x, lab.Inf = ".", ...) {
  o <- attr(x, "order")
  #
  if (!is.null(o))
    x <- x[o, o]
  #
  x[is.infinite(x)] <- lab.Inf
  #
  prmatrix(x, quote = FALSE, right = TRUE)
  #
  invisible(NULL)
}


#' @rdname netdistance
#' @export

netdistance <- function(x, ...)
  UseMethod("netdistance")


designs <- function(treat1, treat2, studlab, sep.trts = ":") {
  
  
  id <- seq_along(studlab)
  o <- order(studlab, treat1, treat2)
  ##
  if (any(o != id)) {
    treat1 <- treat1[o]
    treat2 <- treat2[o]
    studlab <- studlab[o]
    id <- id[o]
  }
  
  
  studies <- unique(studlab)
  n.study <- length(unique(studies))
  ##
  designs <- data.frame(studlab = "", design = rep_len("", n.study))
  
  
  for (i in seq_len(n.study)) {
    designs$studlab[i] <- studies[i]
    designs$design[i] <-
      paste(sort(unique(c(treat1[studlab == studies[i]],
                          treat2[studlab == studies[i]]))),
            collapse = sep.trts)
  }
  
  
  dat <- data.frame(studlab, treat1, treat2, o)
  ##
  res <- merge(dat, designs, by = "studlab")
  res <- res[order(res$o), ]
  res$treat1 <- NULL
  res$treat2 <- NULL
  res$o <- NULL
  ##
  res
}

