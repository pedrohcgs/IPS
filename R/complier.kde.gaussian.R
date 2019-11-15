complier.kde.gaussian <- function(x, p.eval, bw = stats::bw.nrd, adjust = 1,
                                  whs = NULL, gaussian=T,
                                  ...) {
  n <- length(x)
  if(is.null(whs)) whs <- rep(1, n)
  sd <-  bw(x) * adjust
  if (gaussian==T){
  k.gaussian <- function(x, mean=0, sd=1){
    stats::dnorm(x, mean = mean, sd = sd)
  }
  y <-  base::outer(x, p.eval, k.gaussian, sd = sd)
  
  } else {
    k.Epanechnikov <- function(x, mean=0, sd=1) {
      h <- sqrt(5)*sd
      ifelse((z <- abs(x-mean)) < h, 3/4*(1 - (z/h)^2)/h, 0)
    }
    y <-  base::outer(x, p.eval, k.Epanechnikov, sd = sd)
  }
  
  y <- base::colMeans(whs * y)
  
  return(y)
}
