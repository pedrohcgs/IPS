# Compute Quantile for weighted distributions
cdfinv <- function(F.ecdf, probs = c(0, 0.25, 0.5, 0.75, 1)) {
  x <- stats::knots(F.ecdf)
  x <- x[base::order(x)]                              #here
  quantiles <- vapply(probs, FUN = function(p) {x[which(p <= F.ecdf(x))[1]]}, FUN.VALUE = 1.0)
  
  nams <- paste(format(round(probs * 100,
                             if (length(probs) > 1) 2 - log10(diff(range(probs))) else 2)),
                "%", sep = "")
  names(quantiles) <- nams
  return(quantiles)
}
