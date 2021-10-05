library(ggplot2)
library(viridis)
library(tikzDevice)
library(cowplot)
setwd("~/Desktop/coalescence_paper/code")


exponential = function(x, beta){
  return(1/beta * exp(-x/beta))
}
simulations = read.csv('../data/simulation_results_new_cost.csv')
#Delete simulations under structured scenario
simulations_pref = simulations[simulations$K != 0.9,]# & simulations$kc == 0.01),]
#Delete rows containing NAs
sim_099 = simulations_pref[simulations_pref$l == 0.99,]


hist(simulations_pref$Cav, freq = F, breaks = 100, col = 'blue', ylim = c(0, 2.5),
     xlab = 'Competition')
hist(simulations_pref$C0av, freq = F, breaks = 100, add = T, col=rgb(1,0,0,0.5))

hist(simulations_pref$Fav, freq = F, breaks = 100, col = 'blue', ylim = c(0, 2.5),
     xlab = 'Facxilitation')
hist(simulations_pref$F0av, freq = F, breaks = 100, add = T, col=rgb(1,0,0,0.5))

hist(simulations_pref$Fav - simulations_pref$Cav, freq = F, breaks = 100, col = 'blue', ylim = c(0, 2.5),
     xlab = 'Cohesion')
hist(simulations_pref$F0av - simulations_pref$C0av, freq = F, breaks = 100, add = T, col=rgb(1,0,0,0.5))

hist(simulations_pref[simulations_pref$kc == 0.01,]$r, 
     xlim = c(0, 60))
hist(simulations_pref[simulations_pref$kc == 0.455,]$r, 
     add = T)  
hist(simulations_pref[simulations_pref$kc == 0.9,]$r, 
     add = T)

#Plot F and C and color by richness###########################################################################

simulations = read.csv('../data/simulation_results_small_assembly.csv')
simulations_l = simulations[simulations$l == 0.9,]
simulations_l['id'] = paste('(', simulations_l$kc, ', ', simulations_l$kf, ')')
before_after_scat = 
  ggplot(data = simulations_l)+
  geom_point(color = 'grey', pch = 3, size = 1.5, stroke = 3,
             aes(x = C0av, y = F0av))+
  geom_point(pch = 21, size = 6.5,color = 'grey', 
             aes(x = Cav, y = Fav, 
                 fill = as.factor(id)))+
  
  theme(aspect.ratio = 1,
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA,
                                  color = NA),
        legend.key.size = unit(0.7, 'cm'),
        legend.text = element_text(size = 30),
        legend.position = c(0.25, 0.8),
        legend.title = element_text(size = 35, hjust = 0.4,
                                    face = 'bold', color = 'black',
                                    margin = margin(t = 20)),
        legend.spacing.y = unit(1, 'mm'),
        # legend.spacing.y = unit(0, 'mm'),
        axis.title.x = element_text(size = 35),
        axis.title.y = element_text(size = 35),
        axis.text = element_text(size = 35))+
  #scale_color_brewer(palette = 'RdYlBu')
  scale_fill_brewer(palette = 'Set1', direction = -1)+
  labs(color = '$(k_c \\ , \\ k_f)$', 
       fill = '$(k_c \\ , \\ k_f)$',
       x  = 'Competition $(\\mathcal{C})$',
       y = 'Facilitation $(\\mathcal{F})$')+
  guides(fill = guide_legend(override.aes = list(size=10),
                             nrow=4,byrow=TRUE,
                             title.position = 'top',
                             title.vjust = 3,
                             title.hjust = 0.5))
tikz("../sandbox/before_after_scatter.tex",
     width = 8, height = 8,
     standAlone = TRUE)
print(before_after_scat)
dev.off()  

#theme(legend.position = 'none')

#Before and after plot#########################################################################################

n_pref = read.csv('../data/n_preferences_before_after_new.csv', row.names = 1)
#Get max number of preferences
m = max(n_pref$n_pref)
#Aggregate according to simulations
n_pref_grouped = n_pref %>% 
  group_by(n_pref, kc, l)  %>% 
  summarise(p0nr_av = mean(N0nr/m), p0nr_sd = sd(N0nr/m),
            pnr_av = mean(Nnr/r), pnr_sd = sd(Nnr/r), 
            rel.ab_av = mean(rel_ab), rel.ab_sd = sd(rel_ab)) %>%
  filter(kc != 0.85)
n_pref_grouped['err'] = 1/exponential(n_pref_grouped$n_pref, 6)*sqrt(n_pref_grouped$pnr_sd^2 + 
                                                                       n_pref_grouped$p0nr_sd^2)
n_pref_grouped['err_'] = 1/exponential(n_pref_grouped$n_pref, 6)* sqrt((n_pref_grouped$pnr_av - n_pref_grouped$p0nr_av)^2*n_pref_grouped$rel.ab_sd^2 + 
                                                                        n_pref_grouped$rel.ab_av^2*(n_pref_grouped$pnr_sd^2 + 
                                                                         n_pref_grouped$p0nr_sd^2) )
# %>%
#   filter(kc != 0.9)


before_after = 
  ggplot(data = n_pref_grouped)+ 
  geom_line(size = 2.5, aes(x = n_pref,
                          y = (pnr_av-p0nr_av)*1/exponential(n_pref, 6),
                          color = as.factor(kc),
                          linetype = as.factor(l)))+
  geom_abline(size = 2, slope = 0, intercept = 0, color = 'black', linetype = 'longdash')+
  scale_color_brewer(palette = 'Set1', direction = -1)+
  scale_fill_brewer(palette = 'Set1', direction = -1)+
  scale_linetype_manual(values = c('solid', 'longdash'), 
                     labels = c('$0.1$', '$0.9$'))+
  theme(aspect.ratio = 0.666667,
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.key = element_rect(fill = NA,
                                  color = NA),
        legend.key.size = unit(1.2, 'cm'),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 30, face = 'bold'),
        legend.direction = 'vertical',
        legend.box = 'horizontal',
        legend.position = 'none',
        legend.spacing.y = unit(0, 'mm'),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 25))+
    scale_x_continuous(expand = c(0, 0),
                       limits = c(0, 32))+
    scale_y_continuous(limits = c(-0.2, 1.1))+
  labs(x = 'Number of preferences', 
       y = 'Change in proportion after assembly',
       linetype = '$l$',
       color = '\\textbf{$k_c$}')
  
  #Get rid of specialists with 1 or 2 preferences
n_pref_grouped_gen = n_pref_grouped %>% filter(n_pref > 2 & kc != 0.85)
inset_weights = 
    ggplot(data = n_pref_grouped_gen)+ 
    geom_line(size = 1.5, aes(x = n_pref, 
                            y = (pnr_av-p0nr_av)*rel.ab_av*1/exponential(n_pref, 6),
                            color = as.factor(kc),
                            linetype = as.factor(l)))+
    geom_abline(size = 0.75, slope = 0, intercept = 0, color = 'black', linetype = 2)+
    #scale_color_viridis_d(option = 'magma')+
    scale_linetype_manual(values = c('solid', 'longdash'), 
                          labels = c('$0.1$', '$0.9$'))+
    theme(aspect.ratio = 0.666667,
          panel.background = element_rect(fill = NA),
          panel.grid = element_blank(),
          plot.background = element_blank(), #element_rect(fill = NA),
          panel.border = element_rect(fill = NA),
          legend.key = element_rect(fill = NA),
          legend.key.size = unit(1.2, 'cm'),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 25, face = 'bold'),
          legend.direction = 'vertical',
          legend.box = 'horizontal',
          legend.position = 'none',
          legend.spacing.y = unit(0, 'mm'),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_text(size = 20))+
    scale_x_continuous(limits = c(3, 60),
                       breaks = c(3, 30, 60))+
scale_color_brewer(palette = 'Set1', direction = -1)+
      scale_fill_brewer(palette = 'Set1', direction = -1)+
  scale_y_continuous(breaks = 0)+
    labs(linetype = '$l$',
         color = '\\textbf{$k_c$}')
  

  tikz("../sandbox/before_after.tex",
       width = 8, height = 5.33333,
       standAlone = TRUE)
  ggdraw() +
    draw_plot(before_after) +
    draw_plot(inset_weights, x = -0.145, y = 0.45, height = 0.5)
  dev.off()
  
##########Histogram of relative abundances###################################################################
  
  rel_ab = 
    ggplot(data = n_pref_grouped)+
    geom_line(size = 2.5, aes(x = n_pref,
                  y = rel.ab_av,
                  color = as.factor(kc)))+
    geom_abline(size = 0.75, slope = 0, intercept = 0, color = 'black', linetype = 2)+
    theme(aspect.ratio = 1,
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA),
          legend.key = element_rect(fill = NA,
                                    color = NA),
          legend.key.size = unit(1.2, 'cm'),
          legend.text = element_text(size = 30),
          legend.title = element_text(size = 35, face = 'bold'),
          legend.direction = 'vertical',
          legend.box = 'horizontal',
          legend.position = c(0.8, .65),
          legend.spacing.y = unit(0, 'mm'),
          axis.title.x = element_text(size = 35),
          axis.title.y = element_text(size = 35),
          axis.text = element_text(size = 35))+
    scale_x_continuous(limits = c(0, 11),
                       expand = c(0, 0),
                       breaks = c(1, 5, 10))+
    scale_color_brewer(palette = 'Set1', direction = -1)+
  labs(x = '$n_r$',
       y = 'Abundance proportion',
       color = '$k_c$')
  tikz("../sandbox/rel_ab.tex",
       width = 8, height = 8,
       standAlone = TRUE)
  print(rel_ab)
  dev.off()  
  
  
#Histogram about richneses colored by kc######################################################################

  hist_rich =ggplot(data = simulations_pref,
                    aes(x = r))+
    geom_histogram(bins = 20,
                   alpha = 0.5,
                   color = 'grey',
                   aes(fill = as.factor(kc)),
                   position = 'identity')+
    theme(aspect.ratio = 0.666667,
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA),
          legend.key = element_rect(fill = NA,
                                    color = NA),
          legend.key.size = unit(0.7, 'cm'),
          legend.text = element_text(size = 20),
          legend.position = c(0.87, 0.80),
          legend.title = element_text(size = 25, hjust = 0.4,
                                      face = 'bold', color = 'black'),
          legend.background = element_rect(fill = NA),
          # legend.spacing.y = unit(0, 'mm'),
          axis.title.x = element_text(size = 25),
          axis.text.x = element_text(size = 25),
          axis.title.y = element_text(size = 20),
          axis.text.y = element_text(size = 20))+
    scale_fill_brewer(palette = 'Set1', direction = -1)+
    scale_y_continuous(expand = c(0, 0))+
    labs(x = 'Species richness',
         y = 'Counts',
         fill = '$k_c$')
  
  tikz("../sandbox/hist_rich.tex",
       width = 8, height = 5.3333333,
       standAlone = TRUE)
  print(hist_rich)
  dev.off()  
