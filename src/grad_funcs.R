##Functions to plot the magnitude and direction of gradients

#Extract gradient distribution information
extract_gradinfo <- function(x)
{
  tmp = do.call(rbind, x)
  tmp2 = tmp[ , grep('^nabla\\[',colnames(tmp))]
  return(tmp2)
}

#Determine proportion of distribution which is positive
prop_gthan_zero <- function(data) {
  # Count the number of elements greater than zero
  count_positive <- sum(data > 0)
  # Calculate the proportion
  proportion <- count_positive / length(data)
  return(proportion)
}

#Determine proportion of distribution which is above a specified time threshold
prop_gthan_threshold <- function(data, threshold) {
  # Count the number of elements greater than zero
  count_positive <- sum(data > threshold)
  # Calculate the proportion
  proportion <- count_positive / length(data)
  return(proportion)
}