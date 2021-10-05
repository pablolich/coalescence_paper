setwd("~/Desktop/coalescence_paper/code")
#Load libraries
require(ggplot2)
require(RColorBrewer)
require(viridis)
require(itertools)
require(tidyverse)
library(dplyr)
library(ggpubr)
library(gridExtra)
library(reshape2)
library(tikzDevice)
library(ggstar)
#Colour Palettes Based on the Scientific Colour-Maps
library(scico)
source('functions_coalescence.R')

#########################################################################
#FIGURE 3 Results with preferential feeding
#########################################################################

#########################################################################
#A. Community matrix of effective leakage
Deff_micro = as.matrix(read.csv('../data/effective_D.csv', header = F))
dimnames(Deff_micro) = list(as.character(seq(60)), as.character(seq(60)))
melt_Deff_micro = melt(Deff_micro)
Deff_micro = 
  ggplot(data = melt_Deff_micro, 
         aes(x = Var2, y = rev(Var1), fill = value)) +
  geom_tile(aes(col = value))+
  theme(panel.background = element_blank(),
        legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 21,
                                    hjust = 0.8,
                                    margin = margin(t = 0)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 9)),
        plot.margin = margin(b = 10),
        aspect.ratio = 1)+
  scale_color_scico(palette = 'lajolla', direction = -1)+
  scale_fill_scico(palette = 'lajolla', direction = -1)+
  labs(y = ('Species'), x = ('Metabolic by-product'))

#########################################################################
#B. Average community facilitation and competititon

simulations = read.csv('../data/simulation_results.csv')
#Delete simulations under structured scenario
simulations_pref = simulations[simulations$K != 0.9,]# & simulations$kc == 0.01),]
#Delete rows containing NAs
sim_099 = simulations_pref[simulations_pref$l == 0.99,]
#simulations_pref = simulations_pref[!is.na(simulations_pref$CPasc), ]
sim_l = simulations_pref %>% 
  group_by(l) %>%
  summarise(Cav_av = mean(Cav), Cav_se = sd(Cav)/sqrt(max(simulations_pref$n_sim)),
            Fav_av = mean(Fav), Fav_se = sd(Fav)/sqrt(max(simulations_pref$n_sim)))

cols <- c('$\\mathcal{C}$'= "#7E1900","$\\mathcal{F}$"="#1A3399",
          'expression(C[a])'="#AC7825", 'expression(C[b])'="#AC7825")
C_F_pref = 
  ggplot(data = sim_l, aes(x = l))+
  geom_line(aes(y = Cav_av,
                colour = '$\\mathcal{C}$'),
            size = 4)+
  geom_ribbon(data = sim_l, 
              aes(x = l, 
                  ymin = Cav_av - Cav_se,
                  ymax = Cav_av + Cav_se),
              alpha = 0.3, fill = "#7E1900")+
  geom_ribbon(data = sim_l, 
              aes(x = l, 
                  ymin = Fav_av - Fav_se,
                  ymax = Fav_av + Fav_se),
              alpha = 0.3, fill = "#1A3399")+
  geom_line(data = sim_l, aes(y = Fav_av,
                              colour = '$\\mathcal{F}$'),
            size = 4)+
  scale_colour_manual(name="",values=cols)+
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.key = element_rect(fill = NA, 
                                  color = NA),
        #legend.key.size = unit(0.7,'cm'),
        legend.key.size = unit(10, 'mm'),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.text = element_text(size = 20),
        legend.spacing.x = unit(5, 'mm'),
        legend.direction = 'vertical',
        legend.position = c(0.25, 0.70),
        axis.title.x = element_text(size = 21, margin = margin(t = -15)),
        axis.title.y = element_text(size = 21, margin = margin(r = -20)),
        axis.text = element_text(size = 20))+
  scale_x_continuous(breaks = c(0.1, 0.9),
                     labels = c('0.1', '0.9'))+
  scale_y_continuous(breaks = c(0.1, 0.9),
                     labels = c('0.1', '0.9'))+
  labs(x = 'Leakage $(l)$', y = 'Interaction')

#########################################################################
#C. Similarity plot
#Load the results with structure
coal_res_pref = read.csv('../data/coalescence_clean_results_0.00.0depletion_test.csv')
#DATA WRANGLING
#Set parameter grid of analysis
#vec_kc = unique(coal_res_pref$kc)
vec_leak = unique(coal_res_pref$l1)        
vec_points = seq(25)
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
  x = coal_res_pref_i$Cav2- coal_res_pref_i$Cav1 + coal_res_pref_i$Fav1  - coal_res_pref_i$Fav2
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

#Order according to leakage
leakage = factor(df$leakage, levels = rev(levels(as.factor(df$leakage))))
df$leakage = leakage
#Plot
similarity_micro = 
  ggplot(data = df, 
         aes(x = x, y = S))+
  geom_ribbon(aes(ymin = S - S_err, ymax = S + S_err,
                  group = leakage,
                  fill = leakage), 
              alpha = 0.3
  )+
  geom_line(aes(group = leakage,
                color = leakage),
            size = 7)+
  labs(x = paste("Net competition difference"),
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
                                                    r = -2, 
                                                    b = 0, 
                                                    l = 0)),
        axis.text = element_text(size = 20))+
  scale_color_scico_d(palette = 'lajolla')+
  scale_fill_scico_d(palette = 'lajolla')+
  # scale_x_continuous(expand = c(0,0),
  #                    #limits = c(-0.41, 0.53),
  #                    breaks = c(-0.4, 0, 0.5),
  #                    labels = c("-0.4",  "0", "0.5"))+
  # scale_y_continuous(expand = c(0,0),
  #                    breaks = c(-0.5, 0, 0.5),
  #                    labels = c("-0.5",  "0", "0.5"))+
  labs(colour = 'Leakage', fill = 'Leakage')

#########################################################################
#D. Depletion-cohesion plot
#Average over simulations

simulations = read.csv('../data/simulation_results.csv')
simulations_pref = simulations[(simulations$K != 0.9),]
# sim_l = simulations_pref %>% 
#   group_by(l, kc, kf) %>%
#   summarise(Cav_av = mean(Cav),
#             Fav_av = mean(Fav),
#             r_av = mean(r),
#             ER_av = mean(ER))

simulations_pref = simulations_pref[(simulations_pref$ER<1000),]
simulations_pref_01 = simulations_pref[(simulations_pref$l == 0.1),] 
depletion_01 =
  ggplot(data = simulations_pref_01,
         aes(x = (-Fav + Cav)))+
  geom_point(aes(shape = as.factor(kc),
                 y = r,
                 fill = ER/r), 
             color = 'grey',
             # stroke = 0.1
             size = 2.5)+
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = 'black', fill=NA),
        aspect.ratio = 1,
        legend.title = element_text(size = 21),
        legend.key = element_blank(),
        axis.title = element_text(size = 18),
        legend.position = 'right',
        #legend.spacing.x = unit(0,'mm'),
        legend.margin = margin(b = 30),
        legend.direction = 'vertical',
        legend.box = 'horizontal',
        legend.text = element_text(size = 18),
        axis.text = element_text(size = 20),
        axis.title.y = element_text(margin = margin(t = 0, 
                                                    r = 7, 
                                                    b = 0, 
                                                    l = 0),
                                    hjust = -6,
                                    size = 21),
        axis.title.x = element_text(size = 21))+
  scale_fill_scico(palette = 'lajolla',
                   trans = 'reverse', direction = -1)+
  scale_shape_manual(values = c(21, 22, 24),
                     labels = c('0', '0.5', '0.9'))+
  #scale_color_scico(palette = 'lajolla')+
  labs(fill = '$R^{\\star}_{tot}$',
       color = 'Species\nrichness',
       shape = '$k_c$',
       x = '$\\mathcal{C} - \\mathcal{F}$',
       y = 'Species richness ($r$)')+
  # scale_x_continuous(expand = c(0.03, 0),
  #                    breaks = c(-4, 0.5),
  #                    limits = c(-4, 0.6),
  #                    labels = c("-4",  "0.5"))+
  scale_y_continuous(breaks = c(0, 30, 60),
                     labels = c("0", "30", "60"),
                     limits = c(1, 60))+
  scale_x_continuous(breaks = c(0, 2, 4),
                     limits = c(-0.2, 4.2))+
  annotate('text', x = 2, y = 50, label = '$l = 0.1$',
           size = 7)

legend <- get_legend(depletion_01)

depletion_01 <- depletion_01 + theme(legend.position = "none")
#########################################################################
#E Depletion-cohesion plot

simulations = read.csv('../data/simulation_results.csv')
simulations_pref = simulations[(simulations$K != 0.9),]
sim_l = aggregate(simulations_pref,  list(l = simulations_pref$l), FUN = mean)[-1]
my_dat = sim_l %>% select(4, 9, 10, 11)

simulations_pref = simulations_pref[(simulations_pref$ER<1000),]
simulations_pref_09 = simulations_pref[(simulations_pref$l == 0.9),] 
depletion_09 = 
  ggplot(data = simulations_pref_09,
         aes(x = (-Fav + Cav),
             y = r))+
  geom_point(data = simulations_pref_09,
             aes(shape = as.factor(kc),
                 y = r,
                 fill = ER/r), 
             color = 'grey',
             # stroke = 0.1
             size = 2.5)+
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = 'black', fill=NA),
        aspect.ratio = 1,
        axis.title = element_text(size = 21),
        legend.position = 'none',
        axis.text = element_text(size = 20),
        axis.title.y = element_text(margin = margin(t = 0, 
                                                    r = 0, 
                                                    b = 0, 
                                                    l = 0)))+
  scale_fill_scico(palette = 'lajolla',
                   trans = 'reverse', direction = -1)+
  scale_shape_manual(values = c(21, 22, 24))+
  labs(fill = 'Richness',
       colour = 'Richness',
       x = '$\\mathcal{C} - \\mathcal{F}$',
       y = '$r$')+
  # scale_x_continuous(expand = c(0.04, 0),
  #                    breaks = c(-4, 0.5),
  #                    limits = c(-4, 0.6),
  #                    labels = c("-4",  "0.5"))+
  scale_y_continuous(breaks = c(0, 30, 60),
                     labels = c("0", "30", "60"),
                     limits = c(1, 60))+
  scale_x_continuous(breaks = c(0, 2, 4),
                      limits = c(-0.2, 4.2))+
  annotate('text', x = 2, y = 50, label = '$l = 0.9$',
           size = 7)
#########################################################################
#F. Depletion-cohesion plot

simulations = read.csv('../data/simulation_results.csv')
simulations_pref = simulations[(simulations$K != 0.9),]
sim_l = aggregate(simulations_pref,  list(l = simulations_pref$l), FUN = mean)[-1]
my_dat = sim_l %>% select(4, 9, 10, 11)

simulations_pref = simulations_pref[(simulations_pref$ER<1000),]
simulations_pref_05 = simulations_pref[(simulations_pref$l == 0.5),] 
depletion_05 = 
  ggplot(data = simulations_pref_05,
         aes(x = -Fav + Cav,
             y = r))+
  geom_point(data = simulations_pref_05,
             aes(shape = as.factor(kc),
                 y = r,
                 fill = ER/r), 
             color = 'grey',
             # stroke = 0.1
             size = 2.5)+
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = 'black', fill=NA),
        aspect.ratio = 1,
        axis.title = element_text(size = 21),
        legend.position = 'none',
        axis.text = element_text(size = 20))+
  scale_fill_scico(palette = 'lajolla',
                   trans = 'reverse', direction = -1)+
  scale_shape_manual(values = c(21, 22, 24))+
  labs(fill = 'Richness',
       colour = 'Richness',
       x = '     ',
       y = '     ')+
  # scale_x_continuous(expand = c(0.03, 0),
  #                    breaks = c(-4, 0.5),
  #                    limits = c(-4, 0.6),
  #                    labels = c("-4",  "0.5"))+
  scale_y_continuous(expand = c(0.04, 0),
                     limits = c(0, 61),
                     breaks = c(0, 30, 60),
                     labels = c("0",  "30", "60"))+
  scale_x_continuous(breaks = c(0, 2, 4),
                     limits = c(-0.2, 4.2))+
  annotate('text', x = -2, y = 50, label = '$l = 0.5$',
           size = 7)
#########################################################################
#Plot alltogether
lay <- rbind(c(1, 2, 2, 3, 4),
             c(6, 2, 2, 5, 7))
# lay <- rbind(c(1, 1, 2, 2, 2, 2, 3, 3, 4, 4),
#              c(1, 1, 2, 2, 2, 2, 3, 3, 4, 4),
#              c(6, 6, 2, 2, 2, 2, 5, 5, NA, 7),
#              c(6, 6, 2, 2, 2, 2, 5, 5, NA, 7))
#options( tikzLatexPackages = c( getOption( "tikzLatexPackages" ), "\\usepackage{amssymb}"))
tikz("../sandbox/figure_3_old_cost.tex",
     width = 15.33, height = 5.99,
     standAlone = TRUE)
grid.arrange(Deff_micro, similarity_micro, depletion_05, depletion_09, 
             depletion_01, C_F_pref, legend,
             layout_matrix = lay)
dev.off()


