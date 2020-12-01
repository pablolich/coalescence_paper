jaccard = function(x, y){
  #Compute the Jaccard index of binary vectors x and y.
  
  #Get index of presences for each vector
  pre_x = which(x == 1)
  pre_y = which(y == 1)
  #Calculate number of elements in the intersection of presences
  m11 = length(intersect(pre_x, pre_y))
  #Get index of absences for each vector
  ab_x = which(x == 0)
  ab_y = which(y == 0)
  #Number of elements present in x and absent in y
  m10 = length(intersect(pre_x, ab_y))
  #Number of elements absent in x and present in y
  m01 = length(intersect(ab_x, pre_y))
  return(m11/(m10 + m01 + m11))
}

lechon = function(x, y){
  #Get index of presences for each vector
  pre_x = which(x == 1)
  pre_y = which(y == 1)
  #Calculate number of elements in the intersection of presences
  m11 = length(intersect(pre_x, pre_y))
  return(m11)
}

lechon_distribution = function(m, k){
  #Computes the probability of sharing k preferences
  
  #Calculate number of vectors with i preferences, with i in {1, ..., m}
  n_pref = choose(m, 1:m)
  #Calculate the probability of having at least k preferences (k or less)
  p_cum = cumsum(n_pref)
  #Note that I exclude the case of 0 preferences, but add 0 explicitly
  #so that the vector has the same length as the binomial vector, which
  #allows for sharing 0 metabolites
  p_cum = c(0, p_cum)
  #Calculate probability to have at least k preferences (k or more)
  p_k = sum(n_pref) - p_cum[k+1]
  #Calculate probability of sharing k metabolites
  p_share_k = dbinom(x = k, prob = (k/m)^2, size = m)
  #Final probability is probability of having at least k preferences
  #multiplied by probability of sharing k metabolites
  return((p_k)^2*p_share_k)
}

library('gtools')
#Given a vector of indices of lenght m
m = 9
v = 1:m
#Generate permutation binary vectors
perm = permutations(n = 2, r = m, repeats.allowed = T) - 1
#Get rid of the element with all 0
perm = perm[-1,]
#Get indices of an upper diagonal matrix
ind = combinations(n = nrow(perm), r = 2, repeats.allowed = T)
#Prealocate vector of jaccard indices
jaccard_v = rep(0, nrow(ind))
#Calculate jaccard index for all combinations
for (i in seq(nrow(ind))){
  jaccard_v[i] = lechon(perm[ind[i,1],], perm[ind[i,2],])
}
plot(table(jaccard_v)/sum(table(jaccard_v)))
#Get theoretical probabilities
p_theory = rep(0, m)
for (i in 0:m){
  p_theory[i+1] = lechon_distribution(9, i)
}

points(0:9, p_theory/sum(p_theory),col = 'red')


#Get only the vector that have n_r strategies
n_r = 1:m
#Compute vector of row sums
r_sum = rowSums(perm)
#Preallocate matrix of probabilities per number of strategies
p_n_r = matrix(0, nrow = m, ncol = m)
for (i in n_r){
  #Get a vector of booleans where TRUE represents preference vectors
  #with i strategies
  bool_nr = r_sum == i
  #Select rows that have i strategies
  perm_n_r = perm[bool_nr,]
  #Sum columns to view probability distribution
  if (is.integer(nrow(perm_n_r))){
    p_n_r[i,] = colSums(perm_n_r)
  }
  else{
    p_n_r[i,] = perm_n_r
  }
}
