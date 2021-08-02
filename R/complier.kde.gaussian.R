complier.kde.gaussian <- function(x, p.eval, bw = NULL, adjust = 1,
                                  whs = NULL, gaussian=TRUE,
                                  ...) {
  n <- length(x)
  if(is.null(whs)) whs <- rep(1, n)
  
  if(is.null(bw)) {
    bw <- stats::bw.nrd0(x)
  }
  
  if(bw == "nrd0"){
    bw <- stats::bw.nrd0(x)
  } else if(bw == "nrd"){
    bw <- stats::bw.nrd(x)
  } else   if(bw == "ucv"){
    bw <- stats::bw.ucv(x)
  } else if(bw == "bcv"){
    bw <- stats::bw.bcv(x)
  } else if(bw == "SJ"){
    bw <- stats::bw.SJ(x)
  } 
  
  sd <-  bw(x) * adjust
  if (gaussian==TRUE){
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
