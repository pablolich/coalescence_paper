setwd("~/Desktop/coalescence_paper/code")

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
c_mat_no_str = ggplot(data = melt_c_no_str, 
                      aes(x = Var2, y = rev(Var1), fill = value)) +
  geom_tile()+
  theme(panel.background = element_blank(),
        legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20, hjust = 0.5, 
                                  margin = margin(b = -6)),
        plot.margin = margin(b = 20),
        aspect.ratio = 1)+
  scale_color_scico(palette = 'lajolla')+
  scale_fill_scico(palette = 'lajolla')+
  labs(title = expression(paste('C(', k[c], ' = 0, ', K[c], ' = 0)')), x = 'Resource type', y = 'Species')
#Cmats micro struct
#########################################################################################################
melt_c_micro_str = melt(c_micro_str)
c_micro_str = ggplot(data = melt_c_micro_str, 
                     aes(x = Var2, y = rev(Var1), fill = value)) +
  geom_tile(aes(col = value))+
  theme(panel.background = element_blank(),
        legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 20, color = 'white'),
        plot.title = element_text(size = 20, hjust = 0.5, 
                                  margin = margin(b = -6)),
        plot.margin = margin(b = 20),
        aspect.ratio = 1)+
  scale_color_scico(palette = 'lajolla')+
  scale_fill_scico(palette = 'lajolla')+
  labs(title = expression(paste('C(', k[c]!=0, ', ', K[c],' = ', '0)')))

#Cmats macro struct
#########################################################################################################
melt_c_macro_str = melt(c_macro_str)
c_macro_str = ggplot(data = melt_c_macro_str, 
                     aes(x = Var2, y = rev(Var1), fill = value)) +
  geom_tile(aes(col = value))+
  theme(panel.background = element_blank(),
        legend.position = 'none',
        axis.text = element_blank(),
        axis.title = element_text(size = 20, colour = 'white'),
        plot.title = element_text(size = 20, hjust = 0.5, 
                                  margin = margin(b = -6)),
        plot.margin = margin(b = 20),
        axis.ticks = element_blank(),
        aspect.ratio = 1)+
  scale_color_scico(palette = 'lajolla')+
  scale_fill_scico(palette = 'lajolla')+
  labs(title = expression(paste('C(', k[c],' = ', '0', ', ', K[c]!=0, ')')))
#D no struct
#########################################################################################################
melt_D_no_str = melt(D_no_str)
D_no_str = ggplot(data = melt_D_no_str, 
                  aes(x = Var2, y = rev(Var1), fill = value)) +
  geom_tile(aes(col = value))+
  theme(panel.background = element_blank(),
        legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20, hjust = 0.5, 
                                  margin = margin(b = -6)),
        aspect.ratio = 1)+
  scale_color_scico(palette = 'lajolla', direction = -1)+
  scale_fill_scico(palette = 'lajolla', direction = -1)+
  labs(title = expression(paste('D(', k[f], ' = 0,  ', K[f], ' = 0)')), 
       x = 'Metabolic by-product', y = 'Consumed resource')

#D micro struct
#########################################################################################################
melt_D_micro_str = melt(D_micro_str)
D_micro_str = ggplot(data = melt_D_micro_str, 
                     aes(x = Var2, y = rev(Var1), fill = value)) +
  geom_tile(aes(col = value))+
  theme(panel.background = element_blank(),
        legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 20, colour = 'white'),
        plot.title = element_text(size = 20, hjust = 0.5, 
                                  margin = margin(b = -6)),
        aspect.ratio = 1)+
  scale_color_scico(palette = 'lajolla', direction = -1)+
  scale_fill_scico(palette = 'lajolla', direction = -1)+
  labs(title = expression(paste('D(', k[f]!=0, ', ', K[f],' = ', '0)')))

#D macro struct
#########################################################################################################
melt_D_macro_str = melt(D_macro_str)
D_macro_str = ggplot(data = melt_D_macro_str, 
                     aes(x = Var2, y = rev(Var1))) +
  geom_tile(aes(fill = value, col = value),
  )+
  theme(panel.background = element_blank(),
        legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 20, colour = 'white'),
        plot.title = element_text(size = 20, hjust = 0.5, 
                                  margin = margin(b = -6)),
        aspect.ratio = 1)+
  scale_color_scico(palette = 'lajolla', direction = -1)+
  scale_fill_scico(palette = 'lajolla', direction = -1)+
  labs(title = expression(paste('D(', k[f],' = ', '0', ', ', K[f]!=0, ')')))

#Put together
#########################################################################################################

pdf('../sandbox/sampling_scenarios.pdf', width = 12, height = 8.5)
ggarrange(c_mat_no_str, c_micro_str, c_macro_str,
          D_no_str, D_micro_str, D_macro_str,
          labels = c('A', 'B', 'C',
                     'D', 'E', 'F'), 
          label.y = 1, 
          label.x = 0.05,
          font.label = list(size = 20))
dev.off()