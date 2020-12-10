bin_data = function(x, y, n_breaks, centralize = 'mean'){
  #Divide range of data in n_breaks equispaced bins
  step = (max(x)-min(x))/(n_breaks - 1)
  breaks = seq(from = min(x), 
               to = max(x), 
               by = step)
  #Find the bin to which each point of x belongs.
  bin_ind = findInterval(x, breaks)
  #Preallocate vector of medians
  y_binned = rep(0, n_breaks-1)
  sigma = rep(0, n_breaks-1)
  #Calculate average of points in each bin
  for (i in seq(n_breaks-1)){
    #Collect y points that fall inside bin i
    bin_data = y[which(bin_ind == i)]
    #Calculate its centralized measure
    if (centralize == 'mean'){
      y_binned[i] = mean(bin_data)
      } else if (centralize == 'median'){
        y_binned[i] = median(bin_data)
      } else{
        stop("The centralize parameter can only take values 'mean' or 'median' ")
    }
    sigma[i] = sd(bin_data)
  }
  #Create vector of mid points
  mids = breaks[-length(breaks)] + step/2
  #Create output lists
  out = list(y_binned, sigma, mids)
  names(out) = c('binned', 'sigma', 'mids')
  return(out)
}
