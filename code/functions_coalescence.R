bin_data = function(x, y, n_breaks, centralize = 'mean'){
  #This function vertically bins a group of data according to the 
  #specified number of bins, and centralizing measure.
  
  #Divide range of data in n_breaks containing the same number of points
  breaks = best_binning(x, y, n_breaks)
  #Find the bin to which each point of x belongs.
  bin_ind = findInterval(x, breaks)
  
  #Preallocate vector of medians
  y_binned = rep(0, n_breaks-1)
  x_binned = rep(0, n_breaks-1)
  sigma_y = rep(0, n_breaks-1)
  sigma_x = rep(0, n_breaks-1)
  #Calculate average of points in each bin
  for (i in seq(n_breaks-1)){
    #Collect y points that fall inside bin i
    bin_y = y[which(bin_ind == i)]
    bin_x = x[which(bin_ind == i)]
    #Calculate its centralized measure
    if (centralize == 'mean'){
      y_binned[i] = mean(bin_y)
      x_binned[i] = mean(bin_x)
      } else if (centralize == 'median'){
        y_binned[i] = median(bin_y)
        x_binned[i] = median(bin_x)
      } else{
        stop("The centralize parameter can only take values 'mean' or 'median' ")
    }
    sigma_y[i] = sd(bin_y)
    sigma_x[i] = sd(bin_x)
  }
  #Create output lists
  out = list(y_binned, sigma_y, x_binned, sigma_x)
  names(out) = c('y_binned', 'sigma_y', 'x_binned', 'sigma_x')
  return(out)
}

best_binning = function(x, y, n_breaks, tol = 10, alpha = 0.2){
  #Function to determine the partition that makes the number of
  #points in each bin as close as possible
  
  #Total number of points
  P = length(x)
  #Number of bins
  B = n_breaks - 1
  #Number of points per bin to reach
  p_opt = floor(P/B)
  #Range of data
  L = max(x) - min(x)
  #Divide range of data in n_breaks equispaced bins
  #Breaks for a uniform partition
  step = max(x)
  b = seq(from = min(x), 
          to = max(x), 
          by = (max(x)-min(x))/B)
  h = diff(b)
  step = L/B
  for (i in seq(n_breaks-2)){#Subtract two because the extremes are fixed
    #Calculate number of points falling inside ith bin with current
    #partition
    p_bin_i = table(findInterval(x, b, all.inside = T))[i]
    #Calculate difference between initial h and optimal h
    dif1 = abs(p_bin_i - p_opt)
    dif2 = 0
    b0 = b[i]
    dif = tol + 1
    #Modify partitioning iteratively until optimal one is reached for ith bin
    while ( dif > tol){
      #Modify ith break position based on the excess or defect of points per bin
      h[i] = h[i]*(1 + alpha*(p_opt - p_bin_i)/(p_opt + p_bin_i))
      #Recalculate vector of breaks
      if (b0 + h[i] > max(x)){
        #Sometimes the push of the new h[i] exceeds the last element of the breaks
        #vector. In that case, I just set that h[i] as the last element of the vector
        #and finish the binning.
        b = c(b[1:i], b0 + h[i])
        dif = tol - 1
      }
      else{
        b = c(b[1:i], seq(from = b0+h[i], 
                          to = max(x),
                          by = (max(x) - (b0+h[i]))/(B-i)))
        dif = abs(p_bin_i-p_opt)
      }
      #Recalulate number of points in ith bin
      p_bin_i = table(findInterval(x, b, all.inside = T))[i]
    }
  }
  return(b)
}
