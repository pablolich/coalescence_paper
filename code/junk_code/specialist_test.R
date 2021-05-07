rm(list = ls())
setwd("~/Desktop/coalescence_paper/code")

source('functions_coalescence.R')
require(ggplot2)
require(RColorBrewer)
library(ggpubr)

c_mat = read.csv(file = "../data/c_matrices_equal_test_specialist.csv")
simulations = read.csv(file = '../data/simulation_results_equal_test_specialist.csv')

extinct = which(c_mat[,length(c_mat)] == 1)
c_mat_present = c_mat[extinct,]
# n_sim = unique(c_mat_present$X)
# #Preallocate
# av_least = rep(0, length(n_sim))
# av_pop = rep(0, length(n_sim))
# av_pop_before = rep(0, length(n_sim))
# for (i in n_sim){
#   #Get proportion specialists-generalists (spceialists)
#   #Select chunk from before assembly
#   #Get rows of first matrix
#   rows_before = which(c_mat$X == i)
#   c_i_before = c_mat[rows_before[1]:rows_before[length(rows_before)],
#                              2:(length(c_mat)-1)]
#   #Find number of pathways of each species before assembly
#   paths_before = as.numeric(rowSums(c_i))
#   #Compute proportion of sp-gen for the whole community before assembly
#   sp_pop_before = length(which(paths_before<3))
#   gen_pop_before = length(which(paths_before>8))
#   #Create vector of categories before assembly
#   cat_pop_before = c(rep(0, sp_pop_before), rep(1, gen_pop_before))
#   #Calculate proportion of gen-sp before assembly
#   av_pop_before[i+1] = mean(cat_pop)
#   #Get rows of first matrix
#   rows = which(c_mat_present$X == i)
#   #Select chunk
#   c_i = c_mat_present[rows[1]:rows[length(rows)],
#                       2:(length(c_mat_present)-1)]
#   #Find number of pathways of each species
#   paths = as.numeric(rowSums(c_i))
#   #Compute proportion of sp-gen for the whole community
#   sp_pop = length(which(paths<3))
#   gen_pop = length(which(paths>8))
#   #Create vector of categories
#   cat_pop = c(rep(0, sp_pop), rep(1, gen_pop))
#   #Calculate proportion of gen-sp
#   av_pop[i+1] = mean(cat_pop)
#   #Compute the proportion generalists-specialists
#   #Get demand profile of that community
#   demands = colSums(c_i)
#   #Eliminate 0 from the demands vector
#   demands = demands[demands != 0]
#   #Find the least demanded metabolite
#   least_d = names(which.min(demands))
#   #Find who consumes the least demanded metabolite
#   c_least = which(c_i[least_d] == 1)
#   #Find number of pathways of species that consume from the least consumed metabolite
#   paths_least = paths[c_least]
#   #Find how many specialists and generalists are in thi group of species
#   sp = length(which(paths_least<3))
#   gen = length(which(paths_least>8))
#   #Create vector of categorie
#   cat = c(rep(0, sp), rep(1, gen))
#   if (sp == 0 & gen == 0){
#     cat = -1
#   }
#   #Calculate proportion of gen-sp
#   av_least[i+1] = mean(cat)
# }
# #Add this to the simulations dataframe
# simulations$least = av_least
# simulations$pop = av_pop
# simulations$pop_before = av_pop_before
# #Get rid of rows with minus 1
# simulations = simulations[simulations$least != -1,]
#Merge to get values for different kc only
simulations_kc = aggregate(simulations, 
                 by = list(simulations$kc), 
                 FUN = mean)

#################################################################
#check te distributionof preferences before assembly for each kc
#################################################################
kc_vec = unique(simulations$kc)
m = length(c_mat)-2
mats = resources_preferences(kc_vec, m, c_mat, simulations)

#Plot the same thing in ggplot

n_cons_df = melt(mats$n_cons)
n_cons_df = n_cons_df[order(n_cons_df$Var1),]
n_pref_df = melt(mats$n_pref)
n_pref_df = n_pref_df[order(n_pref_df$Var1),]


##Plot n_cons_df
p_cons = ggplot(data = n_cons_df,
       aes(x = Var2, y = value))+
  geom_point(aes(fill = Var1),
             color = 'black',
             pch = 21,
             size = 3)+
  theme(aspect.ratio = 1,
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.position = 'none')+
  #scale_fill_discrete(labels = c(kc_vec[1], kc_vec[4], kc_vec[7]))+
  scale_fill_distiller(palette = 'RdBu')+
  #scale_fill_viridis()+
  labs(x = 'Resources',
       y = 'Number of consumers',
       fill = expression(k[c]))

##Plot n_cons_df
p_pref = ggplot(data = n_pref_df,
       aes(x = Var2, y = value,
           fill = Var1))+
  geom_bar(stat = 'identity',
           position = 'identity',
           alpha = 1)+
  theme(aspect.ratio = 1,
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.position = 'none')+ 
  scale_fill_distiller(palette = 'RdBu')+
  #scale_fill_viridis()+
  labs(x = 'Number of preferences',
       y = 'Frequency',
       fill = expression(k[c]))

###########################################################################################################
#Do the same thing but with the matrix after assembly
###########################################################################################################

mats_after = resources_preferences(kc_vec, m, c_mat_present, simulations)
#Plot the same thing in ggplot

n_cons_df = melt(mats_after$n_cons)
n_cons_df = n_cons_df[order(n_cons_df$Var1),]
n_pref_df = melt(mats_after$n_pref)
n_pref_df = n_pref_df[order(n_pref_df$Var1),]


##Plot n_cons_df
p_cons_after = ggplot(data = n_cons_df,
                aes(x = Var2, y = value))+
  geom_point(aes(fill = Var1),
             color = 'black',
             pch = 21,
             size = 3)+
  theme(aspect.ratio = 1,
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        legend.key.size = unit(2, 'cm'))+ 
  #scale_fill_discrete(labels = c(kc_vec[1], kc_vec[4], kc_vec[7]))+
  scale_fill_distiller(palette = 'RdBu')+
  #scale_fill_viridis()+
  labs(x = 'Resources',
       y = 'Number of consumers',
       fill = expression(k[c]))

##Plot n_cons_df
p_pref_after = ggplot(data = n_pref_df,
                aes(x = Var2, y = value,
                    fill = Var1))+
  geom_bar(stat = 'identity',
           position = 'identity',
           alpha = 1)+
  theme(aspect.ratio = 1,
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        legend.key.size = unit(2, 'cm'))+ 
  scale_fill_distiller(palette = 'RdBu')+
  #scale_fill_viridis()+
  labs(x = 'Number of preferences',
       y = 'Frequency',
       fill = expression(k[c]))


#########################################################################################################
#Arrange all plots in one page
#########################################################################################################
#Store plots in nice figure
ggarrange(p_cons, p_cons_after, 
          p_pref, p_pref_after,
          nrow = 2, ncol = 2) %>%
  ggsave(filename = "../sandbox/n_consumers_n_preferences.pdf",
           width = 15, height = 15)
