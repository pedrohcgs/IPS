complier.kde.gaussian <- function(x, p.eval, bw = "nrd0", adjust = 1,
                                  whs = NULL, gaussian=TRUE,
                                  ...) {
  n <- length(x)
  if(is.null(whs)) whs <- rep(1, n)
  
 
  if(bw == "nrd0"){
    bw_x <- stats::bw.nrd0(x)
  } else if(bw == "nrd"){
    bw_x <- stats::bw.nrd(x)
  } else   if(bw == "ucv"){
    bw_x <- stats::bw.ucv(x)
  } else if(bw == "bcv"){
    bw_x <- stats::bw.bcv(x)
  } else if(bw == "SJ"){
    bw_x <- stats::bw.SJ(x)
  } 
  
  sd <-  bw_x * adjust
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
