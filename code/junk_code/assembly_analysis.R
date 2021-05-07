setwd("~/Desktop/coalescence_paper/code")
#Load libraries
require(ggplot2)
require(RColorBrewer)
require(viridis)
require(reshape2)
library(tidyverse)
library(ggpubr)
#Create a vector of arguments
#args = commandArgs(trailingOnly=TRUE)
#Load data
if (length(args) >= 1){ 
  simulations = read.csv(paste('../data/simulation_results_', args[1], '.csv', 
                               sep = ''))
} else{
  simulations = read.csv('../data/simulation_results.csv')
}
#Sample rows from this huge dataset
rows = sample(1:nrow(simulations), 1e4)
#simulations = simulations[rows,]
#Specify what to group by
groups = list(kc = simulations$kc,
              #kf = simulations$kf,
              #l = simulations$l,
              beta = simulations$beta)
#Perform average of the groups, i.e., average across simulations
sim = aggregate(x = simulations, 
                by = groups,
                FUN = mean)[,(length(groups)+1):(length(simulations)+length(groups))]
#sim = sim[sim$l == 0.3,]
# simulations_l_3 = simulations[simulations$l == 0.3,]
# simulations_l_3['kc_kf'] = simulations_l_3$kc + simulations_l_3$kf

#Plot
if (length(args) >= 1){
  pdf(paste('../sandbox/F_C_diagram', args[1], '.pdf',
            sep = ''))
} else{
  pdf('../sandbox/F_C_diagram.pdf')
}
ggplot(data = simulations, aes(x = Cav, y = Fav, 
                                   fill = as.factor(l),
                                   shape = as.factor(l)))+
  geom_point(colour = 'black',
             stroke = 0.2)+
  # annotate("text", x = 1.6, y = 2.5, 
  #          label = expression(paste(m[j], '=', frac(l[j], paste('(',1-l[j],')')))),
  #          size = 5)+
  labs(x = expression(paste(C[0])),
       y = expression(paste(F[0])))+
  scale_shape_manual(values = c(21, 22, 23, 24, 25, 3, 4))+
  geom_point(data = simulations[table(simulations$r) > 500,], 
             x = simulations[table(simulations$r) > 500,]$Cav, 
             y = simulations[table(simulations$r) > 500,]$Fav, 
             color  = 'black',
             shape = 20,
             size = 0.3)+
  geom_smooth(data = simulations, 
              method=lm,  linetype="dashed",
              se=FALSE, fullrange=TRUE,
              color = 'black', size = 0.5,
              show.legend = FALSE,
              aes(x = C0av, y = F0av,
                  group = as.factor(l)))+
  geom_smooth(aes(x = Cav, y = Fav, 
                  group = as.factor(l)),
              se=FALSE, fullrange = TRUE, 
              method = lm,
              show.legend = FALSE)+
  theme(aspect.ratio = 1)+
  # geom_abline(slope = 0.1/0.9, intercept = 0, 
  #             colour = brewer.pal(n = 5, name = "RdYlBu")[1],
  #             size = 1)+
  # geom_abline(slope = 0.3/0.7, intercept = 0, 
  #             colour = brewer.pal(n = 5, name = "RdYlBu")[2],
  #             size = 1)+
  # geom_abline(slope = 1, intercept = 0, 
  #             colour = brewer.pal(n = 5, name = "RdYlBu")[3],
  #             size = 1)+
  # geom_abline(slope = 0.7/0.3, intercept = 0, 
  #             colour = brewer.pal(n = 5, name = "RdYlBu")[4],
  #             size = 1)+
  # geom_abline(slope = 0.9/0.1, intercept = 0, 
  #             colour = brewer.pal(n = 5, name = "RdYlBu")[5],
  #             size = 1)+
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = 'grey'),
        panel.border  = element_rect(colour = 'black', fill = NA),
        legend.title = element_text(size = 15),
        legend.title.align = 0.5,
        legend.key = element_blank())+
  #scale_fill_viridis()+ 
  xlim(0,3)+
  ylim(0,3) 
dev.off()

#Again
ggplot(data = simulations,
       aes(x = Cav, y = Fav))+
  geom_point(aes(color = as.factor(l)))
#Remove communities with r = 1
matrix_richness =  aggregate(simulations, 
                             by = list(simulations$kc, 
                                       simulations$kf,
                                       simulations$l), 
                             FUN = mean)


#Plot matrix of richness as a function of kc and kf
pdf('../sandbox/kc_kf_harvest_richness.pdf', width = 12, height = 7)
l = unique(matrix_richness$l)
plots = list()
counter = 1
r_rng = range(matrix_richness$r)
r_mid = (r_rng[2]-r_rng[1])/2 + r_rng[1]
for (i in l){
  data = matrix_richness[matrix_richness$l == i,]
  p = ggplot(data = data,
         aes(x = data$kc, y = data$kf,
             fill = r))+
    geom_tile()+
    theme(aspect.ratio = 1,
          panel.border = element_rect(color = 'black', fill = NA),
          axis.title = element_text(size = 20),
          axis.text = element_text(colour = 'black', size = 10),
          plot.title = element_text(size =  17, hjust = 0.5),
          legend.position = 'none',
          plot.margin = margin(0,0,-5,0, "mm"))+
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = c(0, '', 0.5, '',  1),
                       expand = c(0,0))+
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = c(0, '', 0.5, '',  1),
                       expand = c(0,0))+
    scale_fill_distiller(palette = 'RdBu')+
    # scale_fill_gradient2(high="#B2182B", mid="#F7F7F7", low="#2166AC", #colors in the scale
    #                      midpoint=r_mid,    #same midpoint for plots (mean of the range)
    #                      breaks=seq(30,50,10), #breaks in the scale bar
    #                      limits=c(floor(r_rng[1]), ceiling(r_rng[2]))) + #same limits for plots
    labs(title  = paste('l = ', i),
         x = expression(k[c]),
         y = expression(k[f]))
  if (counter == 3){
    p = p + 
      theme(legend.position = 'right',
            legend.title = element_text(hjust = 0.25, size = 20),
                  legend.text = element_text(size = 10)) + 
      labs(fill = 'r')
  }
  plots[[counter]] = p
  counter = counter + 1
  
}

h_rng = range(matrix_richness$B + matrix_richness$Fin)
h_mid = (h_rng[2]-h_rng[1])/2 + h_rng[1]
for (i in l){
  data = matrix_richness[matrix_richness$l == i,]
  p = ggplot(data = data,
             aes(x = data$kc, y = data$kf,
                 fill = B + Fin))+
    geom_tile()+
    theme(aspect.ratio = 1,
          panel.border = element_rect(color = 'black', fill = NA),
          axis.title = element_text(size = 20),
          axis.text = element_text(colour = 'black', size = 10),
          plot.title = element_text(size =  17, hjust = 0.5),
          legend.title = element_text(hjust = 0.25, size = 20),
          legend.text = element_text(size = 10),
          plot.margin = margin(0,0,-5,0, "mm"))+
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = c(0, '', 0.5, '',  1),
                       expand = c(0,0))+
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = c(0, '', 0.5, '',  1),
                       expand = c(0,0))+
    # scale_fill_gradient2(high="#B2182B", mid="#F7F7F7", low="#2166AC", #colors in the scale
    #                      midpoint=h_mid,    #same midpoint for plots (mean of the range)
    #                      breaks=seq(3.8,4.2,0.1), #breaks in the scale bar
    #                      limits=c(3.3, 4.3)) + #same limits for plots
    scale_fill_distiller(palette = 'RdBu')+
    labs(title  = paste('l = ', i),
         fill = 'H',
         x = expression(k[c]),
         y = expression(k[f]))
  plots[[counter]] = p
  counter = counter + 1
  
}
#Store plots in nice figure
ggarrange(plots[[1]], plots[[2]], plots[[3]], 
          plots[[4]], plots[[5]], plots[[6]],
          ncol = 3, nrow = 2,
          labels = c('A', 'B', 'C', 'D', 'E', 'F'),  
          label.x = 0, label.y = 0.91,
          font.label = list(size = 20))
dev.off()


#Eliminate rows with NAs
simulations_clean = simulations[!is.na(simulations$Fav),]
#Get rid of outliers
simulations_clean = simulations_clean[simulations_clean$av_ab < 200 ,]
matrix_richness_clean =  aggregate(simulations_clean, 
                                   by = list(simulations_clean$kc, 
                                             #simulations_clean$kf,
                                             simulations_clean$l), 
                                   FUN = mean)

#Plot matrix of cohesion as a function of kc and kf
pdf('../sandbox/kc_kf_cohesion.pdf')

l = unique(matrix_richness$l)
for (i in l){
  data = matrix_richness_clean[matrix_richness_clean$l ==i,]
  p = ggplot(data = data,
             aes(x = kc, y = kf,
                 fill = Fav - Cav))+
    geom_tile()+
    theme(aspect.ratio = 1)+
    scale_fill_distiller(palette = 'RdBu')+
    labs(title  = paste('l = ', i))
  print(p)
}
dev.off()

################################################################
#Barplot of richness depending on leakage.
richness_barplot = select(simulations, l, r, kc, kf)
data = melt(table(richness_barplot))
ggplot(data, 
       aes(x = r, y = value, 
           fill = as.factor(l)))+
  geom_bar(stat = 'identity')

#######################################################################
#Plot all the richness heatmaps
#######################################################################
pdf('../sandbox/kc_kf_richness.pdf', width = 12, height = 17)
l = unique(matrix_richness$l)
plots = list()
plots[]
counter = 1
r_rng = range(matrix_richness$r)
r_mid = (r_rng[2]-r_rng[1])/2 + r_rng[1]
#Initialize dataframe for storing max richness
df = data.frame(matrix(0, nrow = length(l), ncol = 3, dimnames = list(as.character(l), c('kc', 'kf', 'l'))))
for (i in l){
  data = matrix_richness[matrix_richness$l == i,]
  #Find the maximum richness for each plot
  rowmax = data[which.max(data$r),c('kc','kf', 'l')]
  df[as.character(i),] = rowmax
  p = ggplot(data = data,
             aes(x = data$kc, y = data$kf,
                 fill = r))+
    geom_tile()+
    theme(aspect.ratio = 1,
          panel.border = element_rect(color = 'black', fill = NA),
          axis.title = element_text(size = 20),
          axis.text = element_text(colour = 'black', size = 10),
          plot.title = element_text(size =  17, hjust = 0.5),
          legend.title = element_text(hjust = 0.25, size = 20),
          legend.text = element_text(size = 10),
          plot.margin = margin(0,0,-5,0, "mm"))+
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = c(0, '', 0.5, '',  1),
                       expand = c(0,0))+
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = c(0, '', 0.5, '',  1),
                       expand = c(0,0))+
    scale_fill_distiller(palette = 'RdBu')+
    # scale_fill_gradient2(high="#B2182B", mid="#F7F7F7", low="#2166AC", #colors in the scale
    #                      midpoint=r_mid,    #same midpoint for plots (mean of the range)
    #                      breaks=seq(30,50,10), #breaks in the scale bar
    #                      limits=c(floor(r_rng[1]), ceiling(r_rng[2]))) + #same limits for plots
    labs(title  = paste('l = ', i),
         fill = 'r',
         x = expression(k[c]),
         y = expression(k[f]))
  plots[[counter]] = p
  counter = counter + 1
}
#Store plots in nice figure
ggarrange(plots[[1]], plots[[2]], plots[[3]], 
          plots[[4]], plots[[5]], plots[[6]],
          plots[[7]],plots[[8]],plots[[9]],
          #plots[[10]],
          ncol = 3, nrow = 3,
          #labels = c('A', 'B', 'C', 'D', 'E', 'F','G'),  
          label.x = 0, label.y = 0.91,
          font.label = list(size = 20))
dev.off()


################################################################################################
#Plot pathway distribution along kc for different values of l
l_vec = unique(matrix_richness_clean$l)
for (l in l_vec){
  data_l = matrix_richness_clean[matrix_richness_clean$l == l,]
  plot(data_l$kc,  data_l$av_path,main = paste(l))
}


##########################################################################################
#Plot of facilitation competition as a function of leakage

sim_l = aggregate(simulations,  list(l = simulations$l), FUN = mean)

# plot(sim_l$l, sim_l$C0av,  pch = 20, col  = 'pink', cex = 1.2,
#      ylim = c(0,1), xlim = c(0,1))
plot(sim_l$l, sim_l$Cav + sim_l$Cbav, type = 'l', 
     col = 'red', ylim = c(0.05, 0.25), lty = 'dashed',
     lwd = 3, xlab = 'l', ylab = 'Interaction level', xlim = c(0.1, 0.9))
#points(sim_l$l, sim_l$F0av, pch = 20, col = 'green', cex = 1.2)
#points(sim_l$l, sim_l$C0bav,  pch = 20, col  = 'pink', cex = 1.2)
# abline(a = 0, b = 1, lty = 'dashed', col = 'grey')
# abline(a = 1, b = -1, lty = 'dashed', col = 'grey')
# abline(h = 1, lty = 'dashed', col = 'grey')
# points(sim_l$l, sim_l$Cbav, pch = 24, col = 'pink')
# points(sim_l$l, sim_l$Cav, pch = 24, col = 'pink')
lines(sim_l$l, sim_l$Fav, col = 'green', 
      lty = 'dashed', lwd = 3)
legend(0.5, 0.1, legend=c("C", "F"),
       col=c("red", "green"), lty='dashed', lwd = 3)
#points(sim_l$l, sim_l$C0av + sim_l$C0bav, pch = 20, col = 'red')

sim_kf = aggregate(simulations,  list(l = simulations$kf), FUN = mean)

sim_kf

plot(sim_kf$kf, sim_kf$Cbav )
plot(sim_kf$kf, sim_kf$Fav)
