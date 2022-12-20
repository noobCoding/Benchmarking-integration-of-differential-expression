#Compare x against uniform null distribution y with its range, [y.min, y.max].
ks.test.truncated <- function (x, w.x=NULL, y,..., alternative = c("two.sided", "less", "greater"), exact = NULL, maxCombSize=10000,y.min, y.max) {
  alternative <- match.arg(alternative)
  DNAME <- deparse(substitute(x))
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 1L) 
    stop("not enough 'x' data")
  PVAL <- NULL
  {
    if (is.character(y)) 
      y <- get(y, mode = "function", envir = parent.frame())
    if (!is.function(y)) 
      stop("'y' must be a function or a string naming a valid function")
    METHOD <- "One-sample Kolmogorov-Smirnov test"
    TIES <- FALSE
    if (length(unique(x)) < n) {
      warning("ties should not be present for the Kolmogorov-Smirnov test")
      TIES <- TRUE
    }
    if (is.null(exact)) 
      exact <- (n < 100) && !TIES
    if(is.null(w.x)){
      w.x<-rep(1,length(x))
    }
    names(w.x)=x
    x%<>%as.vector()
    
    # x <- y(sort(x), ...) - (0:(n - 1))/n   
    x.vec<-punif(sort(x),min=y.min,max=y.max)
    
    w.x.vec<-w.x[as.character(sort(x))]%>%as.vector()
    x.substitution<-c(0,cumsum(w.x.vec))/sum(w.x.vec)# 
    
    x<-x.vec-x.substitution[1:(length(x.substitution)-1)]
    STATISTIC <- switch(alternative, two.sided = max(c(x, 
                                                       x.substitution[2:length(x.substitution)]-x.vec)), greater = max(x.substitution[2:length(x.substitution)]-x.vec), less = max(x))
    
    if (exact) {
      PVAL <- 1 - if (alternative == "two.sided")
        result = tryCatch({
          .C(C_pkolmogorov2x, p = as.double(STATISTIC), 
             as.integer(n), PACKAGE = "stats")$p
        }, warning = function(w) {
          warning(w)
        }, error = function(e) {
          .Call(C_pKolmogorov2x, STATISTIC, n)
        }, finally = {
        })
      
      else {
        pkolmogorov1x <- function(x, n) {
          if (x <= 0) 
            return(0)
          if (x >= 1) 
            return(1)
          j <- seq.int(from = 0, to = floor(n * (1 - 
                                                   x)))
          1 - x * sum(exp(lchoose(n, j) + (n - j) * log(1 - 
                                                          x - j/n) + (j - 1) * log(x + j/n)))
        }
        pkolmogorov1x(STATISTIC, n)
      }
    }
    nm_alternative <- switch(alternative, two.sided = "two-sided", 
                             less = "the CDF of x lies below the null hypothesis", 
                             greater = "the CDF of x lies above the null hypothesis")
  }
  names(STATISTIC) <- switch(alternative, two.sided = "D", 
                             greater = "D^+", less = "D^-")
  
  if (is.null(PVAL)) {
    pkstwo <- function(x, tol = 1e-06) {
      if (is.numeric(x)) 
        x <- as.double(x)
      else stop("argument 'x' must be numeric")
      p <- rep(0, length(x))
      p[is.na(x)] <- NA
      IND <- which(!is.na(x) & (x > 0))
      if (length(IND))
        p[IND] <- tryCatch({
          tryRes <- .C(stats:::C_pkstwo, length(x[IND]), p = x[IND], 
                       as.double(tol), PACKAGE = "stats")$p
        }, warning = function(w) {
          warning(w)
        }, error = function(e) {
          tryRes <- .Call(stats:::C_pKS2, p = x[IND], tol)
        }, finally = {
        })
      p
    }
    PVAL <- ifelse(alternative == "two.sided", 1 - pkstwo(sqrt(n) * 
                                                            STATISTIC), exp(-2 * n * STATISTIC^2))
  }
  PVAL <- min(1, max(0, PVAL))
  RVAL <- list(statistic = STATISTIC, p.value = PVAL, alternative = nm_alternative, 
               method = METHOD, data.name = DNAME)
  
  class(RVAL) <- "htest"
  return(RVAL)
}
