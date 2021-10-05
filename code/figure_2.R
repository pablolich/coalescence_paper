setwd("~/Desktop/coalescence_paper/code")
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

#########################################################################################################
## 3. Metabolic preferences and matrix no structure
#########################################################################################################
#Load data
#Cmats no struct
c_no_str = as.matrix(read.csv('../data/c_mat0_struct.csv', header = F))
dimnames(c_no_str) = list(as.character(seq(60)), as.character(seq(60)))
#Cmats micro struct
c_micro_str = as.matrix(read.csv('../data/c_mat2.csv', header = F))
dimnames(c_micro_str) = list(as.character(seq(60)), as.character(seq(60)))
#Cmats macro struct
c_macro_str = as.matrix(read.csv('../data/c_mat2_struct.csv', header = F))
dimnames(c_macro_str) = list(as.character(seq(60)), as.character(seq(60)))
#Dmats no struct
D_no_str = as.matrix(read.csv('../data/D0_struct.csv', header = F))
dimnames(D_no_str) = list(as.character(seq(60)), as.character(seq(60)))
#Dmats micro struct
D_micro_str = as.matrix(read.csv('../data/D_mat2.csv', header = F))
dimnames(D_micro_str) = list(as.character(seq(60)), as.character(seq(60)))
#Dmats macro struct
D_macro_str = as.matrix(read.csv('../data/D2_struct.csv', header = F))
dimnames(D_macro_str) = list(as.character(seq(60)), as.character(seq(60)))
#Plot

#Cmats no struct
#########################################################################################################
melt_c_no_str = melt(c_no_str)
tikz("../sandbox/figure_2A.tex", 
     width = 5, height = 5,
     standAlone = TRUE)
c_mat_no_str = ggplot(data = melt_c_no_str, 
                      aes(x = Var2, y = rev(Var1), fill = value)) +
  geom_tile()+
  theme(panel.background = element_blank(),
        legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 23),
        plot.title = element_text(size = 23, hjust = 0.5, 
                                  margin = margin(b = -6)),
        plot.margin = margin(b = 20),
        aspect.ratio = 1)+
  scale_color_scico(palette = 'lajolla')+
  scale_fill_scico(palette = 'lajolla')+
  labs(title = paste('$C(k_c = 0, K_c = 0)$'), x = 'Resource type', y = 'Species')
print(c_mat_no_str)
dev.off()
#Cmats micro struct
#########################################################################################################
melt_c_micro_str = melt(c_micro_str)
tikz("../sandbox/figure_2B.tex", 
     width = 5, height = 5,
     standAlone = TRUE)
c_micro_str = 
  ggplot(data = melt_c_micro_str, 
                     aes(x = Var2, y = rev(Var1), fill = value)) +
  geom_tile(aes(col = value))+
  theme(panel.background = element_blank(),
        legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 23),
        plot.title = element_text(size = 23, hjust = 0.5, 
                                  margin = margin(b = -6)),
        plot.margin = margin(b = 20),
        aspect.ratio = 1)+
  scale_color_scico(palette = 'lajolla')+
  scale_fill_scico(palette = 'lajolla')+
  labs(title = paste('$C(k_c \\neq 0, K_c = 0)$'),
       x = 'Resource type', y = 'Species')
print(c_micro_str)
dev.off()

#D no struct
#########################################################################################################
melt_D_no_str = melt(D_no_str)
tikz("../sandbox/figure_2D.tex", 
     width = 5, height = 5,
     standAlone = TRUE)
D_no_str = ggplot(data = melt_D_no_str, 
                  aes(x = Var2, y = rev(Var1), fill = value)) +
  geom_tile(aes(col = value))+
  theme(panel.background = element_blank(),
        legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 23),
        plot.title = element_text(size = 23, hjust = 0.5, 
                                  margin = margin(b = -6)),
        aspect.ratio = 1)+
  scale_color_scico(palette = 'lajolla', direction = -1)+
  scale_fill_scico(palette = 'lajolla', direction = -1)+
  labs(title = paste('$D(k_f = 0, K_f = 0)$'), 
       x = 'Metabolic by-product', y = 'Consumed resource')
print(D_no_str)
dev.off()
#D micro struct
#########################################################################################################
melt_D_micro_str = melt(D_micro_str)
tikz("../sandbox/figure_2E.tex", 
     width = 5, height = 5,
     standAlone = TRUE)
D_micro_str = 
  ggplot(data = melt_D_micro_str, 
                     aes(x = Var2, y = rev(Var1), fill = value)) +
  geom_tile(aes(col = value))+
  theme(panel.background = element_blank(),
        legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 23),
        plot.title = element_text(size = 23, hjust = 0.5, 
                                  margin = margin(b = -6)),
        aspect.ratio = 1)+
  scale_color_scico(palette = 'lajolla', direction = -1)+
  scale_fill_scico(palette = 'lajolla', direction = -1)+
  labs(title = paste('$D(k_f \\neq 0, K_f = 0)$'), 
       x = 'Metabolic by-product', y = 'Consumed resource')
print(D_micro_str)
dev.off()


#Put together
#########################################################################################################

pdf('../sandbox/figure_2.pdf', width = 12, height = 8.5)
lay <- rbind(c(1, 2, 3, 4))
grid.arrange(c_mat_no_str, D_no_str, c_micro_str, D_micro_str,
          layout_matrix = lay)
dev.off()

#D macro struct
#########################################################################################################
melt_D_macro_str = melt(D_macro_str)
tikz("../sandbox/figure_2F.tex", 
     width = 5, height = 5,
     standAlone = TRUE)
D_macro_str = ggplot(data = melt_D_macro_str, 
                     aes(x = Var2, y = rev(Var1))) +
  geom_tile(aes(fill = value, col = value),
  )+
  theme(panel.background = element_blank(),
        legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 23, colour = 'white'),
        plot.title = element_text(size = 23, hjust = 0.5, 
                                  margin = margin(b = -6)),
        aspect.ratio = 1)+
  scale_color_scico(palette = 'lajolla', direction = -1)+
  scale_fill_scico(palette = 'lajolla', direction = -1)+
  labs(title = paste('$D(k_f = 0, K_f \\neq 0)$'))
print(D_macro_str)
dev.off()
