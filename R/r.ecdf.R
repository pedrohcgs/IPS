## Function to compute rearranged (monotone) cdfs
r.ecdf=function(q,wr){
  n <- length(q)
  if (n < 1)
    stop("'q' must have 1 or more non-missing values")
  if (all(q == sort(q)) == FALSE)
    stop ("'q' must be sorted beforehand")
  if (min(wr) < 0)
    stop ("'wr' must be non-negative")
  if (max(wr) > 1)
    stop ("'wr' must be not greater than 1")


  rval <- stats::approxfun(q, wr,
                 method = "constant", yleft = 0,
                 yright = 1, f = 0, ties = "ordered")
  class(rval) <- c("ecdf", "stepfun", class(rval))
  assign("nobs", n, envir = environment(rval))
  attr(rval, "call") <- sys.call()
  rval
}
