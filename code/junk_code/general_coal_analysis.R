#Load results
if (length(args) >= 1){
  coalescence_results = read.csv(paste('../data/recursive_diff_coalescence_results', args[1], '.csv', 
                                          sep = ''))
} else{
  coalescence_results = read.csv('../data/recursive_diff_coalescence_results.csv')
}

# #Specify wether you are cohesive or codependant
# coalescence_results['role'] = 0
# for (i in seq(length(coalescence_results$l))){
#   if (coalescence_results$l[i] > 0.5){
#     coalescence_results$role[i] = 1
#   }
# }

# #Calculate cycling factors
 fc = facilitation_cycle(coalescence_results$l1)
 cc = competition_cycle(coalescence_results$l1)
# #recursive_results$Ftot = recursive_results$Ftot/fc
# coalescence_results$coh = coalescence_results$Fav-#/fc - 
#   coalescence_results$Cav#/cc
# mean_rec = aggregate(coh~l+n_coal+role,
#                      coalescence_results,
#                      FUN = mean)
# 
# std_rec = aggregate(coh~l+n_coal+role,
#                     coalescence_resultsc,
#                     FUN = sd)
# 
# mean_rec['err'] = std_rec$coh/sqrt(100)


#################################################################################################
#Final plot with all posble coalescence events between different leakages
l_vec = unique(coalescence_results$l1)
lij = expand.grid(l_vec, l_vec)

#Preallocate matrices of correlaton
l_mat = matrix(0, nrow = length(l_vec), ncol = length(l_vec),
               dimnames = list(as.character(l_vec), as.character(l_vec)))
for (li in seq(length(lij$Var1))){
  #Get data of coalescence experiments for these values of l
  data_i = coalescence_results[(coalescence_results$l1 == lij$Var1[li]) & 
                                 (coalescence_results$l2 == lij$Var2[li]),] 
  #Calculate correlation between cohesion and similarity
  coh = data_i$F1 - data_i$C1 - data_i$F2 + data_i$C2
  #Get indices for the matrix
  i = which(l_vec == lij$Var1[li])
  j = which(l_vec == lij$Var2[li])
  l_mat[i,j] = cor(data_i$S, coh, method = 'spearman')
}

long_data = melt(l_mat)
#Plot 

ggplot(data = long_data, aes(x = Var1, y = Var2))+
  geom_tile(aes(fill = value))+
  scale_y_continuous(trans = "reverse")+
  scale_x_continuous(position = "top") 


coalescence_results[(coalescence_results$l1 == 0.05) &
                    (coalescence_results$l2 == 0.05) &
                    (coalescence_results$n_sim == 0),]


##################################################################################################
#Plot number of extinctions as a function of leakage and cohesion
# coalescence_results = coalescence_results[(coalescence_results$ext1 == 0) | 
#                                           (coalescence_results$ext2 == 0),]
fcc = facilitation_cycle(coal_results$l1)
ccc = competition_cycle(coal_results$l1)

ggplot(data = coal_results, 
       aes(x = Fav1 - Cav1, 
           y = ext1,
           color = as.factor(l1)))+
  geom_point()+
  geom_smooth(method = "lm", fill = NA,
              aes(group = as.factor(l1)),
              color = 'black')

################################################################################################
#Plot average facilitation and facilitation of extincts as a function of leakage

#Get average facilitation of loosing community
nrow = length(coal_results$l1)
F_looser_av = rep(0, nrow)
for (i in seq(nrow)){
  if (coal_results$looser[i] == 1){
    F_looser_av[i] = coal_results$F1[i]
  }
  else if (coal_results$looser[i] == -1){
    F_looser_av[i] = coal_results$F2[i]
  }
  else{
    F_looser_av[i] = NA
  }
    
}
coal_results['F_looser'] = F_looser_av
#Average along leakage
Fext_av = aggregate(Fext ~ l1, #Average facilitation of species from the loosing community that go extinct
                    coal_results,
                    FUN = mean)

F_av = aggregate(Fav ~ l1, #Avergae facilitation of the mix community
                 coal_results,
                 FUN = mean)

F_looser_av = aggregate(F_looser ~ l1, #Average facilitation of the loosing community
                        coal_results,
                        FUN = mean)

fcc = new_fac(unique(coal_results$l1))

#Join both dataframes
dff = data.frame(rbind(as.matrix(F_av), as.matrix(Fext_av), as.matrix(F_looser_av)))
dff['type'] = c(rep('Mix', length(unique(Fext_av$l1))), 
               rep('Extinct', length(unique(Fext_av$l1))),
               rep('Looser', length(unique(Fext_av$l1))))

dff['dfcc'] = rep(fcc, 3)

ggplot(data = dff,
       aes(x = as.character(l1), 
           y = Fav/dfcc,
           fill = as.character(type)))+
  geom_bar(stat='identity', position=position_dodge())+
  theme(aspect.ratio = 1)+
  labs(x = 'l', y = expression(F[av]),
       fill = 'Group of species')+
  scale_fill_grey()

################################################################################################
#Plot average competition and competition of extincts as a function of leakage

#Get average competition of loosing community
nrow = length(coal_results$l1)
C_looser_av = rep(0, nrow)
for (i in seq(nrow)){
  if (coal_results$looser[i] == 1){
    C_looser_av[i] = coal_results$C1[i]
  }
  else if (coal_results$looser[i] == -1){
    C_looser_av[i] = coal_results$C2[i]
  }
  else{
    C_looser_av[i] = NA
  }
  
}
coal_results['C_looser'] = C_looser_av
#Average along leakage
Cext_av = aggregate(Cext ~ l1, #Average facilitation of species from the loosing community that go extinct
                    coal_results,
                    FUN = mean)

C_av = aggregate(Cav ~ l1, #Avergae facilitation of the mix community
                 coal_results,
                 FUN = mean)

C_looser_av = aggregate(C_looser ~ l1, #Average facilitation of the loosing community
                        coal_results,
                        FUN = mean)

ccc = new_comp(unique(coal_results$l1))

#Join both dataframes
df = data.frame(rbind(as.matrix(C_av), as.matrix(Cext_av), as.matrix(C_looser_av)))
df['type'] = c(rep('Mix', length(unique(Cext_av$l1))), 
               rep('Extinct', length(unique(Cext_av$l1))),
               rep('Looser', length(unique(Cext_av$l1))))

df['dccc'] = rep(ccc, 3)
df['Fav'] = dff$Fav

ggplot(data = df,
       aes(x = as.character(l1), 
           y = Cav - Fav,
           fill = as.character(type)))+
  geom_bar(stat='identity', position=position_dodge())+
  theme(aspect.ratio = 1)+
  labs(x = 'l', y = expression(C[av] - F[av]),
       fill = 'Group of species')+
  scale_fill_grey()


#########################################################################################################

l = seq(from=0, to = 1, by = 0.0001)
f_n = new_fac(l)
c_n = new_comp(l)
plot(l, c_n, type = 'l', col = 'red', ylim = c(0,  1.45))
lines(l, f_n, col = 'green')

#BEFORE ASSEMBLY (RANDOM)
mean_f0 = aggregate(F0av~l, simulations, 
          mean)
points(mean_f0$l, mean_f0$F0av, pch = 20, col = 'gray')
mean_c0 = aggregate(C0av~l, simulations, 
                   mean)
points(mean_c0$l, mean_c0$C0av, pch = 18, col = 'gray')


#AFTER ASEMBLY (ASEMBLED)
mean_fa = aggregate(Fav~l, simulations, 
                   mean)
points(mean_fa$l, mean_fa$Fav, pch = 20, col = 'blue')
mean_ca = aggregate(Cav~l, simulations, 
                   mean)
points(mean_ca$l, mean_ca$Cav, pch = 18, col = 'blue')

#AFTER COALESCENCE (COALESCED)

mean_fc = aggregate(Fav~l1, coalescence_results, 
                    mean)
points(mean_fc$l1, mean_fa$Fav, pch = 20, col = 'purple')
mean_cc = aggregate(Cav~l1, coalescence_results, 
                    mean)
points(mean_cc$l1, mean_cc$Cav, pch = 18, col = 'purple')
