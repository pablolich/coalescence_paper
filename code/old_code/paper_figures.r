setwd("~/Desktop/coalescence_paper/code")
#Load libraries
require(ggplot2)
require(RColorBrewer)
require(viridis)
require(itertools)
require(tidyverse)
library(ggpubr)
library(reshape2)
library(tikzDevice)
#Colour Palettes Based on the Scientific Colour-Maps
library(scico)
source('functions_coalescence.R')
#Create a vector of arguments
#args = commandArgs(trailingOnly=TRUE)
#Load data
if (length(args) >= 1){
  coal_results = read.csv(paste('../data/coalescence_results_', args[1], '.csv', 
                                sep = ''))
} else{
  coal_results = read.csv('../data/coalescence_results.csv')
}
coal_results = read.csv('../data/coalescence_clean_results_0.00.0depletion_test.csv')

#DATA WRANGLING
#Set parameter grid of analysis
#vec_kc = unique(coal_results$kc)
vec_leak = unique(coal_results$l1)        
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
df['xC'] = 0
df['xF'] = 0
df['x'] = 0
df['S'] = 0
df['x_errC'] = 0
df['x_errF'] = 0
df['x'] = 0
df['S_err'] = 0
#Get number of total simulations
n_sim = length(grid$leakage)
#Get number of breaks
n_mids = length(vec_points) 
#Iterate over it and generate plotting points
par(mfrow=c(3,3), mar = c(3.8, 3.8, 0, 2))
print('Binning data')
for (i in seq(n_sim)){
  print(i)
  coal_results_i = coal_results[#coal_results$kc == grid$kc[i] & 
    coal_results$l1 == grid$leakage[i],]
  #x = coal_results_i$O1 - coal_results_i$O2
  xC = coal_results_i$Cav2  - coal_results_i$Cav1
  xF = coal_results_i$Fav1 - coal_results_i$Fav2
  x = coal_results_i$Cav2  - coal_results_i$Cav1 + 
      coal_results_i$Fav1 - coal_results_i$Fav2
  plot(x, coal_results_i$S, 
       pch = 21, cex = 0.01,
       xlab = '',
       ylab = '',
       col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
  title(ylab = expression(S[1][',2']),
        xlab = expression(Theta[1]-Theta[2]),
        line=2.1, cex.lab=1.5)
  text(x = 0, y = 0.7, 
       label = paste('l = ', as.character(grid$leakage[i])),
       cex = 2)
  plot(xF, coal_results_i$S,
       pch = 21, cex = 0.01,
       xlab = '',
       ylab = '',
       col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
  title(ylab = expression(S[1][',2']),
        xlab = expression(F[1]-F[2]),
        line=2.1, cex.lab=1.5)
  plot(xC, coal_results_i$S,
       pch = 21, cex = 0.01,
       xlab = '',
       ylab = '',
       col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
  title(ylab = expression(S[1][',2']),
        xlab = expression((C[1]-C[2])),
        line=2.1, cex.lab=1.5)
  
  # #Normalize x
  # # xC = xC/max(abs(xC))
  # # xF = xF/max(abs(xF))
  # # x = x/max(abs(x))
  # # x = coal_results_i$H1 - coal_results_i$H2
  # # x = coal_results_i$r1 - coal_results_i$r2
  # y = coal_results_i$S
  # #grid['cor'][i,] = cor(x, y, method = 'spearman')
  # groupC = bin_data(xC, y, n_mids + 1, centralize = 'mean')
  # groupF = bin_data(xF, y, n_mids + 1, centralize = 'mean')
  # group = bin_data(x, y, n_mids + 1, centralize = 'mean')
  # df[((i-1)*n_mids+1):(i*n_mids),]$xC = groupC$x_binned
  # df[((i-1)*n_mids+1):(i*n_mids),]$xF = groupF$x_binned
  # df[((i-1)*n_mids+1):(i*n_mids),]$x = group$x_binned
  # df[((i-1)*n_mids+1):(i*n_mids),]$S = group$y_binned
  # df[((i-1)*n_mids+1):(i*n_mids),]$x_errC = groupC$sigma_x
  # df[((i-1)*n_mids+1):(i*n_mids),]$x_errF = groupF$sigma_x
  # df[((i-1)*n_mids+1):(i*n_mids),]$x_err = group$sigma_x
  # df[((i-1)*n_mids+1):(i*n_mids),]$S_err = group$sigma_y
  # # df$x = group$x_binned
  # # df$S = group$y_binned
  # # df$x_err = group$sigma_x
  # # df$S_err = group$sigma_y
}

#########################################################################################################
## 2. SIMILARITY MACROSTRUCT
#########################################################################################################
ggplot(data = df, 
       aes(x = xF, y = S))+
  geom_ribbon(aes(ymin = S - S_err, ymax = S + S_err,
                  group = leakage,
                  fill = as.factor(leakage)), 
              alpha = 0.9
  )+
  geom_line(aes(group = leakage,
                color = as.factor(leakage)),
            size = 1.5)+
  labs(x = expression(paste("Parent cohesion difference ", "(", Theta[1]-Theta[2], ")")),
       y = expression(paste("Similarity to parents ", "(", S[1][","][2], ")")))+
  theme(aspect.ratio = 1,
        legend.position = c(0.5, 0.13),
        legend.spacing.y = unit(0,"cm"),
        legend.title = element_text(size = 15, hjust = 0.5),
        legend.text = element_text(size = 14),
        legend.spacing.x = unit(1, 'mm'),
        legend.key.size = unit(8,'mm'),
        legend.key = element_rect(fill = NA, colour = NA),
        legend.direction = 'horizontal',
        legend.background = element_rect(fill = NA),
        panel.background = element_blank(), 
        panel.border = element_rect(color = 'black', fill=NA),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15))+
  scale_color_scico_d(palette = 'lajolla')+
  scale_fill_scico_d(palette = 'lajolla')+
  scale_x_continuous(expand = c(-0.06,0),
                     breaks = c(-0.35, 0, 0.35),
                     labels = c("-0.35",  "0", "0.35"))+
  scale_y_continuous(expand = c(0,0),
                     breaks = c(-0.3, 0, 0.3),
                     labels = c("-0.3",  "0", "0.3"))+
  labs(colour = 'Leakage', fill = 'Leakage')+
  guides(colour = guide_legend(title.position = "bottom"),
         fill = guide_legend(title.position = "bottom", ))

#########################################################################################################
## 1. EFFECTIVE METABOLISM MICROSTRUCT
#########################################################################################################

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
          axis.title.x = element_text(size = 19,
                                      margin = margin(t = 3)),
          axis.title.y = element_text(size = 19,
                                      margin = margin(r = 3)),
          plot.margin = margin(t = 0, b = 5),
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
Deff_micro = plot_matrix(A_reorderd)

#########################################################################################################
## 2. EFFECTIVE METABOLISM MACROSTRUCT
#########################################################################################################

#Effective D  macro-struct
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
        axis.title.x = element_text(size = 19,
                                    margin = margin(t = 10)),
        axis.title.y = element_text(size = 19,
                                    margin = margin(r = 17)),
        plot.margin = margin(b = 10),
        aspect.ratio = 1)+
  scale_color_scico(palette = 'lajolla', direction = -1)+
  scale_fill_scico(palette = 'lajolla', direction = -1)+
  labs(y = ('Species'), x = ('Metabolic by-product'))

#########################################################################################################
## 2. AVERAGE COMPETITION AND FACILITATION MACRO
#########################################################################################################

args = 'gen_struct'
#Load data
if (length(args) >= 1){ 
  simulations = read.csv(paste('../data/simulation_results_', args[1], '.csv', 
                               sep = ''))
} else{
  simulations = read.csv('../data/simulation_results.csv')
}

sim_l = aggregate(simulations,  list(l = simulations$l), FUN = mean)[-1]
my_dat = sim_l %>% select(6, 11, 12, 13)

cols <- c('Competition'= "#7E1900","Facilitation"="#1A3399",
          'expression(C[a])'="#AC7825", 'expression(C[b])'="#AC7825")
C_F_macro = ggplot(data = my_dat, aes(x = l))+
  geom_line(aes(y = Cav + Cbav,
                               colour = 'Competition'),
            size = 1.3)+
  geom_line(data = my_dat, aes(y = Fav,
                               colour = 'Facilitation'),
            #color = "#1A3399", 
            size = 1.3)+
  geom_line(data = my_dat, aes(y = Cbav),
            color = "#AC7825",
            linetype = 'dashed',
            size = 1.1)+
  geom_line(data = my_dat, aes(y = Cav),
            color = "#AC7825",
            size = 1.1,
            linetype = 'dashed')+
  scale_colour_manual(name="",values=cols)+
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.key = element_rect(fill = NA, 
                                  color = NA),
        legend.key.size = unit(0.7,'cm'),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.text = element_text(size = 14),
        legend.spacing.y = unit(-5, 'mm'),
        legend.position = c(0.45, 0.86),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 14))+
  scale_x_continuous(breaks = c(0.1, 0.5, 0.9),
                     labels = c('0.1', '0.5', '0.9'))+
  scale_y_continuous(breaks = c(0, 0.13,0.25),
                     labels = c('0','0.13','0.25'))+
  labs(x = 'Leakage', y = 'Interaction')


#########################################################################################################
## 2. AVERAGE COMPETITION AND FACILITATION MICRO
#########################################################################################################

#Load data
args = 'no_structure'
if (length(args) >= 1){ 
  simulations = read.csv(paste('../data/simulation_results_', args[1], '.csv', 
                               sep = ''))
} else{
  simulations = read.csv('../data/simulation_results.csv')
}

sim_l = aggregate(simulations,  list(l = simulations$l), FUN = mean)[-1]
my_dat = sim_l %>% select(4, 16, 17, 18)

cols <- c('Competition'= "#7E1900","Facilitation"="#1A3399",
          'expression(C[a])'="#AC7825", 'expression(C[b])'="#AC7825")
C_F_micro = ggplot(data = my_dat, aes(x = l))+
  geom_line(aes(y = Cav + Cbav,
                colour = 'Competition'),
            size = 1.3)+
  geom_line(data = my_dat, aes(y = Ftot,
                               colour = 'Facilitation'),
            #color = "#1A3399", 
            size = 1.3)+
  geom_line(data = my_dat, aes(y = Cbav),
            color = "#AC7825",
            linetype = 'dashed',
            size = 1.1)+
  geom_line(data = my_dat, aes(y = Cav),
            color = "#AC7825",
            size = 1.1,
            linetype = 'dashed')+
  scale_colour_manual(name="",values=cols)+
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.key = element_rect(fill = NA, colour = NA),
        legend.key.size = unit(0.7,'cm'),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.text = element_text(size = 15),
        legend.position = c(0.5, 0.18),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15))+
  scale_x_continuous(breaks = c(0.1, 0.5, 0.9),
                     labels = c('0.1', '0.5', '0.9'))+
  scale_y_continuous(breaks = c(0, 0.4, 0.8),
                     labels = c('0','0.4', '0.8'),
                     limits = c(0, 0.83))+
  labs(x = 'Leakage', y = 'Interaction')




#########################################################################################################
## 2. SIMILARITY MICROSTRUCT
########################################################################################################

similarity_micro = ggplot(data = df, 
                           aes(x = x, y = S))+
  geom_ribbon(aes(ymin = S - S_err, ymax = S + S_err,
                  group = leakage,
                  fill = as.factor(leakage)), 
              alpha = 0.3
  )+
  geom_line(aes(group = leakage,
                color = as.factor(leakage)),
            size = 1.5)+
  labs(x = expression(paste("Parent cohesion difference ", "(", Theta[1]-Theta[2], ")")),
       y = expression(paste("Similarity to parents ", "(", S[1][","][2], ")")))+
  theme(aspect.ratio = 1,
        legend.position = c(0.2, 0.70),
        legend.title = element_text(size = 15, hjust = 0.5),
        legend.text = element_text(size = 15),
        legend.key.size = unit(8,'mm'),
        legend.background = element_rect(fill = NA),
        panel.background = element_blank(), 
        panel.border = element_rect(color = 'black', fill=NA),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15))+
  scale_color_scico_d(palette = 'lajolla')+
  scale_fill_scico_d(palette = 'lajolla')+
  scale_x_continuous(breaks = c(-0.47, 0, 0.5),
                     labels = c("-0.5",  "0", "0.5"),
                     expand = c(-0.05,0))+
  scale_y_continuous(breaks = c(-0.5,  0,  0.5),
                     labels = c("-0.5", "0",  "0.5"),
                     expand = c(0.02,0))+
  labs(colour = 'Leakage', fill = 'Leakage')

#########################################################################################################
## 3. FACILITATION BARPLOT
#########################################################################################################

#Load the results with structure
coal_res = read.csv('../data/coalescence_clean_results_0.90.9clean.csv')
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

df = data.frame(rbind(as.matrix(F_av), as.matrix(Fext_av), as.matrix(Flooser)))
df['type'] = c(rep('Mix', length(unique(Fext_av$l1))), 
                rep('Extinct', length(unique(Fext_av$l1))),
                rep('Looser', length(unique(Fext_av$l1))))
tikz("../sandbox/tikz_example.tex", width = 8, height = 4)
ggplot(data = df,
       aes(x = l1, 
           y = Fmix/l1,
           fill = type))+
  geom_bar(stat='identity', position=position_dodge())+
  theme(aspect.ratio = 1,
        legend.position = "top",
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA, color = 'black'))+
  scale_x_continuous(breaks = c(0.1, 0.5, 0.9),
                     labels = c('0.1', '0.5', '0.9'))+
  labs(x = paste('Leakage', '$(l)$'), y = paste("Average intrinsic facilitation", "($\\mathcal{F}/l$)"),
       fill = 'Group of species')+
  scale_fill_grey()
dev.off()

#expression(paste("Average facilitation (", F[av], ")"))

#########################################################################################################
## 3. SIMILARITY MICROSTRUCT
#########################################################################################################

#########################################################################################################
## 1. Put plots together
#########################################################################################################
#FIGURE 1
#########################################################################################################

pdf('../sandbox/micro_strucure_results.pdf', height = 3, width = 10)
ggarrange(Deff_micro,  similarity_micro, C_F_micro,
          widths = c(0.285, 0.30, 0.32),
          ncol = 3, nrow = 1)
dev.off()


#########################################################################################################
#FIGURE 2
#########################################################################################################

pdf('../sandbox/macro_strucure_results.pdf', height = 3, width = 10)
ggarrange(Deff_macro, similarity_macro, C_F_macro,
          widths = c(0.285, 0.30, 0.32),
          ncol = 3, nrow = 1)
dev.off()

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
