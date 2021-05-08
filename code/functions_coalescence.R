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

resources_preferences = function(kc_vec, m, c_mat, simulations){
  #Function to determine two matrices. First, the matrix of number of consumers
  #per resource for each value of kc. Second, frequency of each number of
  #preferences for each value of kc.
  
  #Prealocate a matrix of distributions of preferences
  n_preferences = matrix(0, nrow = length(kc_vec),
                         ncol = as.integer(m/3),
                         dimnames = list(kc_vec))
  #Prealocate a matrix of number of consumers per resource.
  n_consumers = matrix(0, nrow = length(kc_vec),
                       ncol = m, dimnames = list(kc_vec))
  for (j in seq(length(kc_vec))){
    #Get simulation indices for which kc = j
    ind = which(simulations$kc == kc_vec[j])
    #Find the rows of all_c that belong to one of these indices.
    coinc = (c_mat$X + 1)%in%ind
    #Mask c_matrix
    c_kc = c_mat[coinc,] 
    #Calculate the number of preferences of each strain
    n_path = rowSums(c_kc[,2:(length(c_kc)-1)])
    #store in the matrix of preferences
    n_preferences[j,] = table(n_path)
    #Calculate number of rows
    n_sim = unique(c_kc$X)
    #Preallocate matrix of n_consumers per community
    n_cons_community = matrix(0, nrow = length(n_sim),
                              ncol = m)
    #Iterate over rows
    for (i in seq(length(n_sim))){
      #Select rows of one community
      rows = which(c_kc$X == n_sim[i])
      #Subset chunk of interest
      c_ji = c_kc[rows[1]:rows[length(rows)],
                  2:(length(c_kc)-1)]
      #Calculate the number of consumers per resource
      n_consumers_i = sort(as.numeric(colSums(c_ji)), decreasing = T)
      n_cons_community[i,] = n_consumers_i
    }
    #Once I the number of consumers per resource for all communities, I average accross axis = 0
    n_consumers[j,] = colMeans(n_cons_community)
  }
  return(list('n_pref' = n_preferences, 'n_cons' = n_consumers))
}

facilitation_cycle = function(l){
  return(l*(1-l)/(1-l*(1-l)))
}

competition_cycle = function(l){
  return((1-l)/(1-l*(1-l)))
}
