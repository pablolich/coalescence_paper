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
library(ggstar)
#Colour Palettes Based on the Scientific Colour-Maps
library(scico)
source('functions_coalescence.R')


#D. Depletion-cohesion plot
#Average over simulations

simulations = read.csv('../data/simulation_results_samraat_suggestions_3.csv')
simulations_pref = simulations[(simulations$K == 0.9),]
sim_l = aggregate(simulations_pref,  list(l = simulations_pref$l), FUN = mean)[-1]
my_dat = sim_l %>% select(4, 9, 10, 11)

simulations_pref = simulations_pref[(simulations_pref$ER<1000),]
simulations_pref_01 = simulations_pref[(simulations_pref$l == 0.1),] 
depletion_01_st = ggplot(data = simulations_pref_01,
                      aes(x = (Fav/l - (Cav/(1-l) + Cbav/l))))+
  geom_point(aes(shape = as.factor(kc),
                 y = ER/r,
                 fill = r), 
             size = 2.5,
             color = 'grey',
             stroke = 0.1)+
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = 'black', fill=NA),
        aspect.ratio = 1,
        legend.title = element_text(size = 15, hjust = 0.5),
        legend.key = element_blank(),
        axis.title = element_text(size = 21),
        legend.position = 'right',
        legend.spacing.x = unit(0,'mm'),
        legend.margin = margin(t = -50),
        legend.direction = 'vertical',
        legend.box = 'horizontal',
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        plot.title = element_text(size = 21))+
        # axis.title.y = element_text(margin = margin(t = 0, 
        #                                             r = 7, 
        #                                             b = 0, 
        #                                             l = 0),
        #                            hjust = -0.3))+
  scale_fill_scico(palette = 'lajolla', direction = -1)+
  scale_shape_manual(values = c(21, 22, 24),
                     labels = c('0', '0.5', '0.9'))+
  #scale_color_scico(palette = 'lajolla', direction = -1)+
  labs(fill = 'Species\nrichness',
       color = 'Species\nrichness',
       shape = '$k_c$',
       title = '$l = 0.1$',
       x = '   ',
       y = '$R^{\\bigstar}_{tot}/r$')+
  scale_x_continuous(expand = c(0.03, 0),
                     breaks = c(-0.5, 0),
                     limits = c(-0.5, 0.13),
                     labels = c("-0.5",  "0"))+
  scale_y_continuous(expand = c(0.03, 0),
                     limits = c(0, 20),
                     breaks = c(0, 20),
                     labels = c("0",  "20"))
  # annotate('text', x = -2.9, y = 3.86, label = '$l = 0.1$',
  #          size = 10)

legend <- get_legend(depletion_01_st)

depletion_01_st <- depletion_01_st + theme(legend.position="none")
#########################################################################
#E Depletion-cohesion plot

simulations = read.csv('../data/simulation_results_samraat_suggestions_3.csv')
simulations_pref = simulations[(simulations$K == 0.9),]
sim_l = aggregate(simulations_pref,  list(l = simulations_pref$l), FUN = mean)[-1]
my_dat = sim_l %>% select(4, 9, 10, 11)

simulations_pref = simulations_pref[(simulations_pref$ER<1000),]
simulations_pref_09 = simulations_pref[(simulations_pref$l == 0.9),] 
depletion_09_st = ggplot(data = simulations_pref_09,
                      aes(x = (Fav/l - (Cav/(1-l) + Cbav/l)),
                          y = ER/r))+
  geom_point(data = simulations_pref_09,
             aes(shape = as.factor(kc),
                 y = ER/r,
                 fill = r), 
             size = 2.5,
             color = 'grey',
             stroke = 0.1)+
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = 'black', fill=NA),
        aspect.ratio = 1,
        axis.title = element_text(size = 21),
        legend.position = 'none',
        axis.text = element_text(size = 15),
        plot.title = element_text(size = 21))+
        # axis.title.y = element_text(margin = margin(t = 0, 
        #                                             r = -10, 
        #                                             b = 0, 
        #                                             l = 0)))+
  scale_fill_scico(palette = 'lajolla', direction = -1)+
  scale_color_scico(palette = 'lajolla', direction = -1)+
  scale_shape_manual(values = c(21, 22, 24))+
  labs(fill = 'Richness',
       colour = 'Richness',
       title = '$l = 0.9$',
       x = '    ',
       y = '    ')+
scale_x_continuous(expand = c(0.03, 0),
                   breaks = c(-0.5, 0),
                   limits = c(-0.5, 0.13),
                   labels = c("-0.5",  "0"))+
  scale_y_continuous(expand = c(0.03, 0),
                     limits = c(0, 20),
                     breaks = c(0, 20),
                     labels = c("0",  "20"))
  # annotate('text', x = -2.9, y = 6.43, label = '$l = 0.9$',
  #          size = 10)
#########################################################################
#F. Depletion-cohesion plot

simulations = read.csv('../data/simulation_results_samraat_suggestions_3.csv')
simulations_pref = simulations[(simulations$K == 0.9),]
sim_l = aggregate(simulations_pref,  list(l = simulations_pref$l), FUN = mean)[-1]
my_dat = sim_l %>% select(4, 9, 10, 11)

simulations_pref = simulations_pref[(simulations_pref$ER<1000),]
simulations_pref_05 = simulations_pref[(simulations_pref$l == 0.5),] 
depletion_05_st = ggplot(data = simulations_pref_05,
                      aes(x = (Fav/l - (Cav/(1-l) + Cbav/l)),
                          y = ER/r))+
  geom_point(data = simulations_pref_05,
             aes(shape = as.factor(kc),
                 y = ER/r,
                 fill = r), 
             size = 2.5,
             color = 'grey',
             stroke = 0.1)+
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = 'black', fill=NA),
        aspect.ratio = 1,
        axis.title = element_text(size = 21),
        legend.position = 'none',
        axis.text = element_text(size = 15),
        plot.title = element_text(size = 21))+
        # axis.title.y = element_text(margin = margin(t = 0, 
        #                                             r = -10, 
        #                                             b = 0, 
        #                                             l = 0)))+
  scale_fill_scico(palette = 'lajolla', direction = -1)+
  #scale_color_scico(palette = 'lajolla', direction = -1)+
  scale_shape_manual(values = c(21, 22, 24))+
  labs(fill = 'Richness',
       colour = 'Richness',
       x = '$\\hat{\\Theta}$',
       title = '$l = 0.5$',
       y = '     ')+
  scale_x_continuous(expand = c(0.03, 0),
                   breaks = c(-0.5, 0),
                   limits = c(-0.5, 0.13),
                   labels = c("-0.5",  "0"))+
  scale_y_continuous(expand = c(0.03, 0),
                     limits = c(0, 20),
                     breaks = c(0, 20),
                     labels = c("0",  "20"))
  # annotate('text', x = -2.9, y = 4.29, label = '$l = 0.5$',
  #          size = 10)

#Plot alltogether
lay <- rbind(c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4))
options( tikzLatexPackages = c( getOption( "tikzLatexPackages" ), "\\usepackage{amssymb}"))
tikz("../sandbox/figure_S4.tex",
     width = 9, height = 3,
     standAlone = TRUE)
grid.arrange(depletion_01_st, depletion_05_st, depletion_09_st, legend,
             layout_matrix = lay)
dev.off()
