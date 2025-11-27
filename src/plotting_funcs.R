#Function to extract quantiles from MCMC output 
extract <- function(x)
{
  tmp = do.call(rbind, x)
  tmp2 = tmp[ , grep('^a\\[',colnames(tmp))]
  qta = apply(tmp2, 2, quantile, prob=c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1))
  return(qta)
}

#Function to plot posterior bars
post.bar <- function(x, i, h, a, barcolours=c("skyblue","dodgerblue","darkblue","darkgreen"))
{
  rect(xleft = x[2], xright = x[6], ybottom = i - h/5, ytop = i + h/5, border = NA, col = barcolours[1]) # 95% interval rectangle
  segments(x[2], i-h/3.5, x[2], i+h/3.5, lwd = 2, col = barcolours[1]) # horizontal ticks for 95%
  segments(x[6], i-h/3.5, x[6], i+h/3.5, lwd = 2, col = barcolours[1]) # horizontal ticks for 95%
  
  rect(xleft=x[3], xright=x[5], ybottom=i-h/3, ytop=i+h/3, border=NA, col=barcolours[2]) #50% interval
  segments(x[3], i-h/2.5, x[3], i+h/2.5, lwd = 2, col = barcolours[2])   # horizontal ticks for 50%
  segments(x[5], i-h/2.5, x[5], i+h/2.5, lwd = 2, col = barcolours[2])   # horizontal ticks for 50%
  
  points(x[4], i, pch = 16, col = barcolours[3], cex = 2) #posterior median 
  points(a, i, pch = 4, col = barcolours[4], cex = 2, lwd = 2) #simulated (true) arrival time 
}
