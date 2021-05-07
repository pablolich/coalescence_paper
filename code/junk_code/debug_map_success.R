##################################################################################
#Heatmap of coalescence success as a function of kc and kf.
##################################################################################
#Initialize data frame
comb = permutations(n = 3,  r = 3, repeats.allowed = T)
l_vec = unique(coal_results$l1)
kc_vec = unique(coal_results$kc1)
kf_vec = unique(coal_results$kf1)
success_map = data.frame(l = l_vec[comb[,1]], 
                         kc = kc_vec[comb[,2]], 
                         kf = kf_vec[comb[,3]],
                         n_success = 0)
#Get values of leakage
l_values = unique(coal_results$l1)
coal_results_l = coal_results[coal_results$l1 == l,]
for (l in l_values){
  #Get number of coalescence events
  n_coal = length(coal_results_l$kc1)
  for (i in seq(n_coal)){
    #Get similarity of the outcome community
    S = coal_results_l$S[i]
    #Store constants of the winning community
    if (S > 0){
      #Winning community is 1
      kc = coal_results_l$kc1[i]
      kf = coal_results_l$kf1[i]
      #Find that kc-kf coordinate in the map
      point = as.numeric(rownames(success_map[(success_map$l == l &
                                                 success_map$kc == kc &
                                                 success_map$kf == kf),]))
      #Add 1 to the n_succes of that point
      success_map$n_success[point] = success_map$n_success[point] + 1
    }else if (S < 0){
      #Winning community is 2
      kc = coal_results_l$kc2[i]
      kf = coal_results_l$kf2[i]
      point = as.numeric(rownames(success_map[(success_map$l == l &
                                                 success_map$kc == kc &
                                                 success_map$kf == kf),]))
      #Add 1 to the n_succes of that point
      success_map$n_success[point] = success_map$n_success[point] + 1
    }else{
      #There is no winning community
      next
    }
  }
}

#Plot this in the heat map
for (i in l_vec){
  success_l = success_map[success_map$l == i,]
  success_plot = ggplot(data = success_l,
                        aes(x = kc, y = kf,
                            fill = n_success))+
    geom_tile()
  print(success_plot)
}
