#Function to compute weigthed cdfs
w.ecdf=function(q, w){
  n <- length(q)
  if (n < 1)
    stop("'q' must have 1 or more non-missing values")
  if (all(q == sort(q)) == FALSE)
    stop ("'q' must be sorted beforehand")
  if (n != length(w))
    stop ("'q' and 'w' must have the same length")
  
  vals <- unique(q)
  if (anyDuplicated(q)) {
    w <- tapply(w, q, sum)
  }
  fn <- cumsum(w)
  rval <- stats::approxfun(vals, fn,
                           method = "constant", yleft = 0,
                           yright = 1, f = 0, ties = "ordered")
  class(rval) <- c("ecdf", "stepfun", class(rval))
  assign("nobs", n, envir = environment(rval))
  attr(rval, "call") <- sys.call()
  rval
}
