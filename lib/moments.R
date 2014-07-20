
sample_moments <- function(vec, order.max=5, simplify=T) {
  moms <- all.moments(vec, order.max)
  names(moms) <- 0:order.max
  
  cums <- all.cumulants(moms)
  names(cums) <- 0:order.max
  
  pros <- c(mean=mean(vec), var=var(vec), # note: mean and var like mom[1:2]
            skewness=skewness(vec), kurtosis=kurtosis(vec), kurtosis.geary=geary(vec))
  
  if (simplify) {
    ret <- c(moms[-1], cums[-(1:2)], pros)
    names(ret) <- c(sprintf("moment.%i", 1:order.max), sprintf("cumulant.%i", 2:order.max), names(pros))
  } else {
    ret <- list(moments=moms, cumulants=cums, other=pros)
  }
  
  return(ret)
}


unnamed <- function(x, y, to.scale=F) {
  if (to.scale) x   <- scale(x)
  if (to.scale) y   <- scale(y)
  x.y <- lm(x ~ y)$residuals + mean(x)
  y.x <- lm(y ~ x)$residuals + mean(y)
  
  df <- data.frame(
    x.y = sample_moments(x.y), 
    x   = sample_moments(x), 
    y.x = sample_moments(y.x), 
    y   = sample_moments(y)
  )
  
  df$x.y_x = df$x.y - df$x
  df$y.x_y = df$y.x - df$y
  df$x_to_y = df$x.y_x - df$y.x_y
  df$y_to_x = df$y.x_y - df$x.y_x
  
  return(df)
}

