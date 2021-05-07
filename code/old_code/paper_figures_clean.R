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
#FIGURE 4 Results with guild structure
#########################################################################

#########################################################################
#A. Community matrix of effective leakage
Deff_macro = as.matrix(read.csv('../data/effective_D_struct.csv', header = F))
dimnames(Deff_macro) = list(as.character(seq(60)), as.character(seq(60)))
melt_Deff_macro = melt(Deff_macro)
Deff_macro = ggplot(data = melt_Deff_macro, 
                    aes(x = Var2, y = rev(Var1), fill = value)) +
  geom_tile(aes(col = value))+
  theme(panel.background = element_blank(),
        legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 21,
                                    hjust = 0.8,
                                    margin = margin(t = 9)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 9)),
        plot.margin = margin(b = 10),
        aspect.ratio = 1)+
  scale_color_scico(palette = 'lajolla', direction = -1)+
  scale_fill_scico(palette = 'lajolla', direction = -1)+
  labs(y = ('Species'), x = ('Metabolic by-product'))

#########################################################################
#B. Average community facilitation and competititon
simulations = read.csv('../data/simulation_results_clean.csv')
simulations_guilds = simulations[simulations$K == 0.9,]
sim_l = aggregate(simulations_guilds,  list(l = simulations_guilds$l), FUN = mean)[-1]
my_dat = sim_l %>% select(5, 10, 11, 12)

cols <- c('$\\mathcal{C}$'= "#7E1900","$\\mathcal{F}$"="#1A3399",
          'expression(C[a])'="#AC7825", 'expression(C[b])'="#AC7825")
C_F_macro = ggplot(data = my_dat, aes(x = l))+
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
        legend.key.size = unit(0.7,'cm'),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.text = element_text(size = 20),
        legend.spacing.y = unit(-5, 'mm'),
        legend.position = c(0.45, 0.80),
        axis.title = element_text(size = 21),
        axis.title.y = element_text(margin = margin(r = -20)),
        axis.text = element_text(size = 20))+
  scale_x_continuous(breaks = c(0.1, 0.9),
                     labels = c('0.1', '0.9'))+
  scale_y_continuous(breaks = c(0, 0.25),
                     labels = c('0', '0.25'))+
  labs(x = 'Leakage $(l)$', y = 'Interaction')

#########################################################################
#C. Similarity plot
#Load the results with structure
coal_res_guilds = read.csv('../data/coalescence_clean_results_0.90.9clean.csv')
#DATA WRANGLING
#Set parameter grid of analysis
#vec_kc = unique(coal_res_guilds$kc)
vec_leak = unique(coal_res_guilds$l1)        
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
  coal_res_guilds_i = coal_res_guilds[#coal_res_guilds$kc == grid$kc[i] & 
    coal_res_guilds$l1 == grid$leakage[i],]
  #x = coal_res_guilds_i$O1 - coal_res_guilds_i$O2
  x = coal_res_guilds_i$Fav1 + coal_res_guilds_i$Cav2 - coal_res_guilds_i$Fav2 - coal_res_guilds_i$Cav1
  #Normalize x
  x = x/max(abs(x))
  # x = coal_res_guilds_i$H1 - coal_res_guilds_i$H2
  # x = coal_res_guilds_i$r1 - coal_res_guilds_i$r2
  y = coal_res_guilds_i$S
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
similarity_macro = ggplot(data = df, 
                          aes(x = x, y = S))+
  geom_ribbon(aes(ymin = S - S_err, ymax = S + S_err,
                  group = leakage,
                  fill = as.factor(leakage)), 
              alpha = 0.3
  )+
  geom_line(aes(group = leakage,
                color = as.factor(leakage)),
            size = 8)+
  labs(x = paste("Cohesion difference", "($\\Theta_1 - \\Theta_2$)"),
       y = paste("Similarity to parents", "($S_{1, 2}$)"))+
  theme(aspect.ratio = 1,
        legend.position = "top",
        legend.spacing.y = unit(0,"cm"),
        legend.title = element_text(size = 18, hjust = 0.5),
        legend.text = element_text(size = 15),
        legend.spacing.x = unit(1, 'mm'),
        legend.key.size = unit(11,'mm'),
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
        axis.title.x = element_text(hjust = 0.8),
        axis.text = element_text(size = 20))+
  scale_color_scico_d(palette = 'lajolla')+
  scale_fill_scico_d(palette = 'lajolla')+
  scale_x_continuous(expand = c(-0.1,-0.1),
                     breaks = c(-0.3, 0, 0.3),
                     labels = c("-0.3",  "0", "0.3"))+
  scale_y_continuous(expand = c(0,0),
                     breaks = c(-0.3, 0, 0.3),
                     labels = c("-0.3",  "0", "0.3"))+
  labs(colour = 'Leakage', fill = 'Leakage')
  # guides(colour = guide_legend(title.position = "bottom"),
  #        fill = guide_legend(title.position = "bottom", ))

#########################################################################
#D. Histogram of facilitation links being disrupted
#Average facilitation of species from the loosing community that go extinct
Fext_av = aggregate(F_loos_ext ~ l1, 
                    coal_res,
                    FUN = mean)
#Average facilitation of species from the loosing community 
Flooser = aggregate(F_looser ~ l1, 
                    coal_res,
                    FUN = mean)
#Avergae facilitation of the mix community
F_av = aggregate(Fmix ~ l1, 
                 coal_res,
                 FUN = mean)

df_fac = data.frame(rbind(as.matrix(F_av), as.matrix(Fext_av), as.matrix(Flooser)))
df_fac['type'] = c(rep('Mix', length(unique(Fext_av$l1))), 
               rep('Extinct', length(unique(Fext_av$l1))),
               rep('Looser', length(unique(Fext_av$l1))))
coextinctions = ggplot(data = df_fac,
       aes(x = l1, 
           y = Fmix/l1,
           fill = type))+
  geom_bar(stat='identity', position=position_dodge())+
  theme(aspect.ratio = 1,
        legend.position = 'top',
        legend.title = element_text(size = 21, hjust = 0.5),
        legend.text = element_text(size = 21),
        legend.key.size = unit(8,'mm'),
        legend.key = element_rect(fill = NA, colour = NA),
        legend.direction = 'horizontal',
        legend.background = element_rect(fill = alpha('grey', 0.4)),
        panel.background = element_blank(), 
        panel.border = element_rect(color = 'black', fill=NA),
        axis.title = element_text(size = 21),
        axis.text = element_text(size = 20))+
  
  scale_x_continuous(breaks = c(0.1, 0.5, 0.9),
                     labels = c('0.1', '0.5', '0.9'))+
  labs(x = paste('Leakage', '$(l)$'), y = paste("Average intrinsic facilitation", "($\\mathcal{F}/l$)"),
       fill = 'Group of species')+
  scale_fill_grey()
  # guides(colour = guide_legend(title.position = "bottom"),
  #        fill = guide_legend(title.position = "bottom", ))

lay <- rbind(c(1, 2, 2, 3, 3),
             c(4, 2, 2, 3, 3))
tikz("../sandbox/figure_4.tex", width = 15.33, height = 5.99)
grid.arrange(Deff_macro, similarity_macro, coextinctions, C_F_macro,
                   layout_matrix = lay)
dev.off()


# #tikz("../sandbox/tikz_example.tex", width = 8, height = 4)
# pdf('../sandbox/micro_strucure_results.pdf', height = 3, width = 10)
# ggarrange(Deff_micro,  similarity_micro, C_F_micro,
#           widths = c(0.285, 0.30, 0.32),
#           ncol = 3, nrow = 1)
# dev.off()
# dev.off()

