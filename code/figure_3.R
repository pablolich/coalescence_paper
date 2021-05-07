setwd("~/Desktop/coalescence_paper/code")
#Load libraries
require(ggplot2)
require(RColorBrewer)
require(viridis)
require(itertools)
require(tidyverse)
library(ggpubr)
library(gridExtra)
library(reshape2)
library(tikzDevice)
#Colour Palettes Based on the Scientific Colour-Maps
library(scico)
source('functions_coalescence.R')

#########################################################################
#FIGURE 3 Results with preferential feeding
#########################################################################

#########################################################################
#A. Community matrix of effective leakage

#Function to reorder rows and colums
plot_matrix <- function(A){
  n <- nrow(A)
  A <- as.data.frame(A)
  A$row <- n:1
  colnames(A) <- c(1:n, "row")
  dA <- A %>% gather(col, value, -row)
  pl <- ggplot(dA, aes(x = as.numeric(row), y = as.numeric(col), fill = value)) + 
    geom_tile(aes(color = value)) + xlab("Metabolic by-product") + 
    ylab("Species") +
    theme(panel.background = element_blank(),
          legend.position = 'none',
          axis.text = element_text(color = 'white'),
          axis.ticks  = element_blank(),
          axis.title.x = element_text(size = 21,
                                      hjust = 0.1,
                                      margin = margin(t = -10)),
          axis.title.y = element_text(size = 21,
                                      margin = margin(r = -10)),
          plot.margin = margin(b = 10),
          aspect.ratio = 1)+
    scale_color_scico(palette = 'lajolla', direction = -1)+
    scale_fill_scico(palette = 'lajolla', direction = -1)
  return(pl)
}

#Effective D  micro-struct
Deff_micro = as.matrix(read.csv('../data/effective_D.csv', header = F))
dimnames(Deff_micro) = list(as.character(seq(60)), as.character(seq(60)))
# calculate eigenvalues/vectors of A and A^t
eA <- eigen(Deff_micro)
etA <- eigen(t(Deff_micro))
# use the first eigenvector of A to order the rows, 
# the first of eigenvector of A^t to order the cols
order_rows <- order(eA$vectors[,1])
order_cols <- order(etA$vectors[,1])
A_reorderd <- Deff_micro[order_rows, order_cols]
#Plot matrix
Deff_micro = plot_matrix(A_reorderd)


#########################################################################
#B. Average community facilitation and competititon

simulations = read.csv('../data/simulation_results.csv')
simulations_pref = simulations[simulations$K != 0.9,]
sim_l = aggregate(simulations_pref,  list(l = simulations_pref$l), FUN = mean)[-1]
my_dat = sim_l %>% select(4, 9, 10, 11)


cols <- c('$\\mathcal{C}$'= "#7E1900","$\\mathcal{F}$"="#1A3399",
          'expression(C[a])'="#AC7825", 'expression(C[b])'="#AC7825")
C_F_pref = ggplot(data = my_dat, aes(x = l))+
  geom_line(aes(y = Cav + Cbav,
                colour = '$\\mathcal{C}$'),
            size = 4)+
  geom_line(data = my_dat, aes(y = Fav,
                               colour = '$\\mathcal{F}$'),
            #color = "#1A3399", 
            size = 4)+
  geom_line(data = my_dat, aes(y = Cbav),
            color = "#AC7825",
            linetype = 'dashed',
            size = 3)+
  geom_line(data = my_dat, aes(y = Cav),
            color = "#AC7825",
            size = 3,
            linetype = 'dashed')+
  scale_colour_manual(name="",values=cols)+
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.key = element_rect(fill = NA, 
                                  color = NA),
        #legend.key.size = unit(0.7,'cm'),
        legend.key.width = unit(10, 'mm'),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.text = element_text(size = 20),
        legend.spacing.x = unit(2, 'mm'),
        legend.direction = 'horizontal',
        legend.position = c(0.55, 0.9),
        axis.title.x = element_text(size = 21, margin = margin(t = -15)),
        axis.title.y = element_text(size = 21, margin = margin(r = -20)),
        axis.text = element_text(size = 20))+
  scale_x_continuous(breaks = c(0.1, 0.9),
                     labels = c('0.1', '0.9'))+
  scale_y_continuous(breaks = c(0.1, 0.7),
                     labels = c('0.1', '0.7'))+
  labs(x = 'Leakage $(l)$', y = 'Interaction')+
  annotate("text", x = 0.7, y = 0.20, label = '$C_a$', size = 6)+
  annotate("text", x = 0.7, y = 0.30, label = '$C_b$', size = 6)
    

#########################################################################
#C. Similarity plot
#Load the results with structure
coal_res_pref = read.csv('../data/coalescence_clean_results_0.00.0depletion_test.csv')
#DATA WRANGLING
#Set parameter grid of analysis
#vec_kc = unique(coal_res_pref$kc)
vec_leak = unique(coal_res_pref$l1)        
vec_points = seq(20)
#Create grid of parameters
grid = expand.grid(list(#kc = vec_kc, 
  leakage = vec_leak))
#Make room for correlation vectors of each case
grid['cor'] = 0
#Create data frame with that grid
df = data.frame(expand.grid(list(point = vec_points, 
                                 #kc = vec_kc, 
                                 leakage = vec_leak)))
#Add columns to dataframe for plotting purposes
df['x'] = 0
df['S'] = 0
df['x_err'] = 0
df['S_err'] = 0
#Get number of total simulations
n_sim = length(grid$leakage)
#Get number of breaks
n_mids = length(vec_points) 
#Iterate over it and generate plotting points
print('Binning data')
for (i in seq(n_sim)){
  print(i)
  coal_res_pref_i = coal_res_pref[#coal_res_pref$kc == grid$kc[i] & 
    coal_res_pref$l1 == grid$leakage[i],]
  #x = coal_res_pref_i$O1 - coal_res_pref_i$O2
  x = coal_res_pref_i$Fav1 + coal_res_pref_i$Cav2 - coal_res_pref_i$Fav2 - coal_res_pref_i$Cav1
  #Normalize x
  x = x/max(abs(x))
  # x = coal_res_pref_i$H1 - coal_res_pref_i$H2
  # x = coal_res_pref_i$r1 - coal_res_pref_i$r2
  y = coal_res_pref_i$S
  grid['cor'][i,] = cor(x, y, method = 'spearman')
  group = bin_data(x, y, n_mids + 1, centralize = 'mean')
  df[((i-1)*n_mids+1):(i*n_mids),]$x = group$x_binned
  df[((i-1)*n_mids+1):(i*n_mids),]$S = group$y_binned
  df[((i-1)*n_mids+1):(i*n_mids),]$x_err = group$sigma_x
  df[((i-1)*n_mids+1):(i*n_mids),]$S_err = group$sigma_y
  # df$x = group$x_binned
  # df$S = group$y_binned
  # df$x_err = group$sigma_x
  # df$S_err = group$sigma_y
}

#Plot
similarity_micro = ggplot(data = df, 
                          aes(x = x, y = S))+
  geom_ribbon(aes(ymin = S - S_err, ymax = S + S_err,
                  group = leakage,
                  fill = as.factor(leakage)), 
              alpha = 0.3
  )+
  geom_line(aes(group = leakage,
                color = as.factor(leakage)),
            size = 7)+
  labs(x = paste("Cohesion difference", "($\\Theta_1 - \\Theta_2$)"),
       y = paste("Similarity to parents", "($S_{1, 2}$)"))+
  theme(aspect.ratio = 1,
        legend.position = "top",
        legend.spacing.y = unit(0,"cm"),
        legend.title = element_text(size = 21, hjust = 0.5, margin = margin(r = 7)),
        legend.text = element_text(size = 21),
        legend.spacing.x = unit(1, 'mm'),
        #legend.key.size = unit(11,'mm'),
        legend.key.width = unit(18, 'mm'),
        legend.key = element_rect(fill = NA, colour = NA),
        legend.direction = 'horizontal',
        legend.background = element_rect(fill = NA),
        legend.margin = margin(b =  0),
        panel.background = element_blank(), 
        panel.border = element_rect(color = 'black', fill=NA),
        axis.title = element_text(size = 21),
        axis.title.y = element_text(margin = margin(t = 0, 
                                                    r = -10, 
                                                    b = 0, 
                                                    l = 0)),
        axis.text = element_text(size = 20))+
  scale_color_scico_d(palette = 'lajolla')+
  scale_fill_scico_d(palette = 'lajolla')+
  scale_x_continuous(expand = c(-0.045,0),
                     breaks = c(-0.4, 0, 0.45),
                     labels = c("-0.4",  "0", "0.5"))+
  scale_y_continuous(expand = c(0,0),
                     breaks = c(-0.5, 0, 0.5),
                     labels = c("-0.5",  "0", "0.5"))+
  labs(colour = 'Leakage', fill = 'Leakage')


#########################################################################
#D. Depletion-cohesion plot
#Average over simulations
simulations_pref = simulations_pref[(simulations_pref$ER<1000),]
simulations_pref = simulations_pref[(simulations_pref$l == 0.1) |
                                      (simulations_pref$l == 0.5) |
                                      (simulations_pref$l == 0.9), ]
depletion = ggplot(data = simulations_pref,
           aes(x = (Fav/l - (Cav/(1-l) + Cbav/l)),
               y = ER,
               fill = as.factor(l)))+
      geom_point(color = 'grey', pch = 21, size = 2)+
      theme(panel.background = element_blank(),
            panel.border = element_rect(color = 'black', fill=NA),
            aspect.ratio = 1,
            legend.title = element_text(size = 21, margin = margin(r = 7)),
            legend.key = element_blank(),
            axis.title = element_text(size = 21),
            legend.position = 'top',
            legend.spacing.x = unit(3,'mm'),
            legend.margin = margin(b =  0),
            legend.direction = 'horizontal',
            legend.text = element_text(size = 21),
            axis.text = element_text(size = 20),
            axis.title.y = element_text(margin = margin(t = 0, 
                                                        r = -10, 
                                                        b = 0, 
                                                        l = 0),
                                        hjust = 0.3))+
        scale_fill_scico_d(palette = 'lajolla', direction = -1)+
      labs(fill = 'Leakage',
           x = 'Intinsic cohesion ($\\hat{\\Theta}$)',
           y = 'Total resource abundance ($R^{\\bigstar}_{tot}$)')+
      scale_x_continuous(expand = c(0.03, 0),
                         breaks = c(-4, 0),
                         limits = c(-4, 0.1),
                         labels = c("-4",  "0"))+
      scale_y_continuous(expand = c(0.03, 0),
                         breaks = c(0, 550),
                         labels = c("0",  "550"))+
      guides(fill = guide_legend(override.aes = list(size=5)))

#########################################################################
#Plot alltogether
lay <- rbind(c(1, 2, 2, 3, 3),
             c(4, 2, 2, 3, 3))
options( tikzLatexPackages = c( getOption( "tikzLatexPackages" ), "\\usepackage{amssymb}"))
tikz("../sandbox/figure_3.tex", 
     width = 15.33, height = 5.99,
     standAlone = TRUE)
grid.arrange(Deff_micro, similarity_micro, depletion, C_F_pref,
             layout_matrix = lay)
dev.off()


