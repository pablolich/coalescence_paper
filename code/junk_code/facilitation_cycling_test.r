test = function(l, n){
  return((l*(1-l))^n)
}

n_thresh = function(L, l){
  return(log(L)/(log(l)+log(1-l)))
}

#Initialize the matrix
l_vec = 0.1*seq(10)
n_vec = seq(from= 0, to = 4, by = 0.005)
l_n = length(n_vec)
data = matrix(data = 0, nrow = l_n, ncol = 10)
#Populate the matrix
for (j in seq(length(l_vec))){
  for (i in seq(length(n_vec))){
    data[i, j] = test(l_vec[j], n_vec[i])
  }
}


plot(n_vec, cumsum(data[,5]/cumsum(data[,5])[l_n]), type = 'l', lwd = 3)
points(n_vec, cumsum(data[,1]/cumsum(data[,1])[l_n]), type = 'l', lwd = 3, col = 'purple')
points(n_vec, cumsum(data[,2]/cumsum(data[,2])[l_n]), type = 'l', lwd = 3, col = 'orange')
points(n_vec, cumsum(data[,3]/cumsum(data[,3])[l_n]), type = 'l', lwd = 3, col = 'red')
points(n_vec, cumsum(data[,4]/cumsum(data[,4])[l_n]), type = 'l', lwd = 3, col = 'green')
points(n_vec, cumsum(data[,6]/cumsum(data[,6])[l_n]), type = 'l', lwd = 3, lty = 'dotted')
points(n_vec, cumsum(data[,7]/cumsum(data[,7])[l_n]), type = 'l', lwd = 3, lty = 'dotted')
points(n_vec, cumsum(data[,8]/cumsum(data[,8])[l_n]), type = 'l', lwd = 3, lty = 'dotted')
points(n_vec, cumsum(data[,9]/cumsum(data[,9])[l_n]), type = 'l', lwd = 3,  lty = 'dotted')


n_t = n_thresh(L = 1e-4, l = l_vec[1:9])
plot(l_vec[1:9], n_t, pch = 20, cex = 2, col = 'red')
