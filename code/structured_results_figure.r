#Figure to generate results of structured scenario
require(ggplot2)
library(scico)
library(tikzDevice)
library(reshape2)

setwd("~/Desktop/coalescence_paper/code")


#######################
#First, plot the matrices C, D and their product CD as an example

#Matrix C

#Load Cmats macro struct
c_macro_str = as.matrix(read.csv('../data/c_mat2_struct.csv', header = F))
dimnames(c_macro_str) = list(as.character(seq(60)), as.character(seq(60)))

#Plot Cmats macro struct
#########################################################################################################
melt_c_macro_str = melt(c_macro_str)
tikz("../sandbox/c_macro_str.tex", 
     width = 5, height = 5,
     standAlone = TRUE)
c_macro_str = 
ggplot(data = melt_c_macro_str, 
       aes(x = Var2, y = rev(Var1), fill = value)) +
  geom_tile(aes(col = value))+
  theme(panel.background = element_blank(),
        legend.position = 'none',
        axis.text = element_blank(),
        axis.title = element_text(size = 21),
        plot.title = element_text(size = 21, hjust = 0.5),
        plot.margin = margin(b = 20),
        axis.ticks = element_blank(),
        aspect.ratio = 1)+
  scale_color_scico(palette = 'lajolla')+
  scale_fill_scico(palette = 'lajolla')+
  labs(title = paste('$C(k_c = 0, K_c \\neq 0)$'),
       x = 'Resource type', y = 'Species')
print(c_macro_str)
dev.off()

#Matrix D

#Load Dmats macro struct
D_macro_str = as.matrix(read.csv('../data/D2_struct.csv', header = F))
dimnames(D_macro_str) = list(as.character(seq(60)), as.character(seq(60)))

#Plot D macro struct
#########################################################################################################
melt_D_macro_str = melt(D_macro_str)
tikz("../sandbox/D_macro_str.tex", 
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
        axis.title = element_text(size = 21),
        plot.title = element_text(size = 21, hjust = 0.5),
        aspect.ratio = 1)+
  scale_color_scico(palette = 'lajolla', direction = -1)+
  scale_fill_scico(palette = 'lajolla', direction = -1)+
  labs(title = paste('$D(k_f = 0, K_f \\neq 0)$'), 
       x = 'Metabolic by-product', y = 'Consumed resource')
print(D_macro_str)
dev.off()

#Matrix CD
#Plot Community matrix of effective leakage

#Load
Deff_macro = as.matrix(read.csv('../data/effective_D_struct.csv', header = F))
dimnames(Deff_macro) = list(as.character(seq(60)), as.character(seq(60)))
melt_Deff_macro = melt(Deff_macro)

tikz("../sandbox/Deff_macro.tex", 
     width = 5, height = 5,
     standAlone = TRUE)
#Plot
Deff_macro = ggplot(data = melt_Deff_macro, 
                    aes(x = Var2, y = rev(Var1), fill = value)) +
  geom_tile(aes(col = value))+
  theme(panel.background = element_blank(),
        legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 21),
        axis.title.y = element_text(size = 21),
        aspect.ratio = 1,
        plot.title = element_text(size = 21, hjust = 0.5))+
  scale_color_scico(palette = 'lajolla', direction = -1)+
  scale_fill_scico(palette = 'lajolla', direction = -1)+
  labs(title = '$D_{eff} = CD$', y = 'Species', x = 'Metabolic by-product')
print(Deff_macro)
dev.off()


#######################

#Load data for depletion plots

simulations = read.csv('../data/simulation_results_block_averaging.csv.csv')
simulations_pref = simulations[(simulations$K == 0.9),]
sim_l = aggregate(simulations_pref,  list(l = simulations_pref$l), FUN = mean)[-1]
my_dat = sim_l %>% select(4, 9, 10, 11)

#Second, plot the resource depletion plots for l = 0.1, 0.5, 0.9
simulations_pref = simulations_pref[(simulations_pref$ER<1000),]
simulations_pref_01 = simulations_pref[(simulations_pref$l == 0.1),] 
tikz("../sandbox/depletion_01_st.tex", 
     width = 5, height = 5,
     standAlone = TRUE)
depletion_01_st = 
ggplot(data = simulations_pref_01,
       aes(x = (-Fav + Cav)))+
  geom_point(aes(shape = as.factor(kc),
                 y = r,
                 fill = ER/r), 
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
        axis.text = element_text(size = 21),
        plot.title = element_text(size = 21))+
  scale_fill_scico(palette = 'lajolla')+
  scale_shape_manual(values = c(21, 22, 24),
                     labels = c('0', '0.5', '0.9'))+
  #scale_color_scico(palette = 'lajolla', direction = -1)+
  labs(fill = '$R^{\\star}_{tot}/r',
       color = '$R^{\\star}_{tot}/r',
       shape = '$k_c$',
       title = '$l = 0.1$',
       x = '   ',
       y = 'Species richness $(r)$')
  # scale_x_continuous(expand = c(0.03, 0),
  #                    breaks = c(-0.5, 0),
  #                    limits = c(-0.5, 0.13),
  #                    labels = c("-0.5",  "0"))+
  # scale_y_continuous(expand = c(0.03, 0),
  #                    limits = c(0, 20),
  #                    breaks = c(0, 20),
  #                    labels = c("0",  "20"))
# annotate('text', x = -2.9, y = 3.86, label = '$l = 0.1$',
#          size = 10)

legend <- get_legend(depletion_01_st)


depletion_01_st <- depletion_01_st + theme(legend.position="none")

print(depletion_01_st)
dev.off()

tikz("../sandbox/legend.tex", 
     width = 5, height = 5,
     standAlone = TRUE)
print(legend)
dev.off()

  #########################################################################
tikz("../sandbox/depletion05.tex", 
     width = 5, height = 5,
     standAlone = TRUE)
  #F. Depletion-cohesion plot for l = 0.5
  simulations_pref_05 = simulations_pref[(simulations_pref$l == 0.5),] 
  depletion_05_st = 
    ggplot(data = simulations_pref_05,
                           aes(x = (-Fav +Cav)))+
    geom_point(data = simulations_pref_05,
               aes(shape = as.factor(kc),
                   y = r,
                   fill = ER/r), 
               size = 2.5,
               color = 'grey',
               stroke = 0.1)+
    theme(panel.background = element_blank(),
          panel.border = element_rect(color = 'black', fill=NA),
          aspect.ratio = 1,
          axis.title = element_text(size = 21),
          legend.position = 'none',
          axis.text = element_text(size = 21),
          plot.title = element_text(size = 21))+
    # axis.title.y = element_text(margin = margin(t = 0, 
    #                                             r = -10, 
    #                                             b = 0, 
    #                                             l = 0)))+
    scale_fill_scico(palette = 'lajolla')+
    #scale_color_scico(palette = 'lajolla', direction = -1)+
    scale_shape_manual(values = c(21, 22, 24))+
    labs(fill = 'Richness',
         colour = 'Richness',
         x = '$\\mathcal{C} - \\mathcal{F}$',
         title = '$l = 0.5$',
         y = '     ')
    # scale_x_continuous(expand = c(0.03, 0),
    #                    breaks = c(-0.5, 0),
    #                    limits = c(-0.5, 0.13),
    #                    labels = c("-0.5",  "0"))+
    # scale_y_continuous(expand = c(0.03, 0),
    #                    limits = c(0, 20),
    #                    breaks = c(0, 20),
    #                    labels = c("0",  "20"))
  # annotate('text', x = -2.9, y = 4.29, label = '$l = 0.5$',
  #          size = 10)
    print(depletion_05_st)
    dev.off()
  #########################################################################
  #E Depletion-cohesion plot for l = 0.9
    tikz("../sandbox/depletion09.tex", 
         width = 5, height = 5,
         standAlone = TRUE)
  simulations_pref_09 = simulations_pref[(simulations_pref$l == 0.9),] 
  depletion_09_st = 
  ggplot(data = simulations_pref_09,
         aes(x = (-Fav + Cav)))+
    geom_point(data = simulations_pref_09,
               aes(shape = as.factor(kc),
                   y = r,
                   fill = ER/r), 
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
    scale_fill_scico(palette = 'lajolla')+
    scale_color_scico(palette = 'lajolla')+
    scale_shape_manual(values = c(21, 22, 24))+
    labs(fill = 'Richness',
         colour = 'Richness',
         title = '$l = 0.9$',
         x = '    ',
         y = '    ')
  # scale_x_continuous(expand = c(0.03, 0),
  #                    breaks = c(-0.5, 0),
  #                    limits = c(-0.5, 0.13),
  #                    labels = c("-0.5",  "0"))+
  # scale_y_continuous(expand = c(0.03, 0),
  #                    limits = c(0, 20),
  #                    breaks = c(0, 20),
  #                    labels = c("0",  "20"))
  # annotate('text', x = -2.9, y = 6.43, label = '$l = 0.9$',
  #          size = 10)
  
  print(depletion_09_st)
  dev.off()
  #########################################################################
  #Similarity plot
  
  #Load the results with structure
  coal_res_guilds = read.csv('../data/coalescence_clean_results_0.9block_averaging.csv')
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
    x = coal_res_guilds_i$Fav2  - coal_res_guilds_i$Fav1+ coal_res_guilds_i$Cav1 - coal_res_guilds_i$Cav2
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
  
  tikz("../sandbox/similarity_plot.tex", 
       width = 5, height = 5,
       standAlone = TRUE)
  #Plot
  similarity_macro = 
    ggplot(data = df, 
                            aes(x = x, y = S))+
    geom_ribbon(aes(ymin = S - S_err, ymax = S + S_err,
                    group = leakage,
                    fill = as.factor(leakage)), 
                alpha = 0.3
    )+
    geom_line(aes(group = leakage,
                  color = as.factor(leakage)),
              size = 8)+
    labs(x = paste("Net competition difference"),
         y = paste("Similarity to parents", "($S_{1, 2}$)"))+
    theme(aspect.ratio = 1,
          legend.position = "top",
          legend.spacing.y = unit(0,"cm"),
          legend.title = element_text(size = 21, hjust = 0.5, margin = margin(r = 7)),
          legend.text = element_text(size = 17),
          legend.spacing.x = unit(1, 'mm'),
          #legend.key.size = unit(11,'mm'),
          legend.key.width = unit(15, 'mm'),
          legend.key = element_rect(fill = NA, colour = NA),
          legend.direction = 'horizontal',
          legend.background = element_rect(fill = NA),
          legend.margin = margin(b =  0),
          panel.background = element_blank(), 
          panel.border = element_rect(color = 'black', fill=NA),
          axis.title = element_text(size = 21),
          axis.text = element_text(size = 20))+
    scale_color_scico_d(palette = 'lajolla')+
    scale_fill_scico_d(palette = 'lajolla')+
    # scale_x_continuous(expand = c(-0.1,-0.1),
    #                    breaks = c(-0.3, 0, 0.3),
    #                    labels = c("-0.3",  "0", "0.3"))+
    # scale_y_continuous(expand = c(0,0),
    #                    breaks = c(-0.3, 0, 0.3),
    #                    labels = c("-0.3",  "0", "0.3"))+
    labs(colour = 'Leakage', fill = 'Leakage')
    
    
    print(similarity_macro)
    dev.off()
    ###############################################################################################\
    ###Plot all together
    
  