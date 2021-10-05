setwd("~/Desktop/coalescence_paper/code")
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(tikzDevice)
library(RColorBrewer)
coal_rec = read.csv('../data/recursive_results_original.csv',row.names = 1 )
#Create unique id for each group
coal_rec['id'] = paste(as.character(coal_rec$l1), as.character(coal_rec$l2), 
                       sep = '')
div_rows = which(coal_rec$Rav == 0)
coal_rec = coal_rec[-div_rows,]
#Number of groups
groups = unique(coal_rec$id)
n_groups = length(groups)
#Group by id
grouped_results <- coal_rec %>% group_by(n_coal, id) %>% 
                                summarise(Ctav = mean(Cav_i), Ctsd = sd(Cav_i),
                                          Ftav = mean(Fav_i), Ftsd = sd(Fav_i),
                                          Rtav = mean(Rav), Rtsd = sd(Rtav),
                                          ninv_av = mean(n_inv), ninv_sd = sd(n_inv),
                                          ltav = mean(lav),
                                          costav = mean(cost), cost_sd = sd(cost),
                                          r_av = mean(r), r_sd = sd(r)) %>% 
                                filter(!(n_coal == 0 & id == '0.10.9')) %>%
                                arrange(id, n_coal)

interactions = 
  ggplot(data = grouped_results)+
    # geom_ribbon(aes(x = n_coal, ymin = (Ftav - Ftsd)/ltav, ymax = (Ftav + Ftsd)/ltav,
    #                 fill = factor(id, 
    #                                levels = c("0.10.1",  "0.90.1","0.10.9", "0.90.9"))), alpha = 0.1)+
    # geom_ribbon(aes(x = n_coal, ymin = Ctav - Ctsd, ymax = Ctav + Ctsd,
    #                 fill = factor(id, 
    #                               levels = c("0.10.1",  "0.90.1","0.10.9", "0.90.9"))), alpha = 0.1)+
    geom_line(size = 2, aes(x = n_coal, y = Ftav/ltav,
                  color = factor(id, 
                                 levels = c("0.10.1",  "0.90.1","0.10.9", "0.90.9"))))+
    geom_line(size = 2, aes(x = n_coal, y =Ctav,
                  color = factor(id, 
                                 levels = c("0.10.1",  "0.90.1","0.10.9", "0.90.9"))),
            linetype = 'longdash')+
    geom_segment(aes(x = 6, y = 0.4, xend = 8.5, yend = 0.4), size = 2)+
    geom_segment(aes(x = 6, y = 0.33, xend = 8.5, yend = 0.33), linetype = 'longdash', size = 2)+
    annotate(geom="text", x=15, y=0.4, label="Facilitation ($\\mathcal{F}/l$)")+
    annotate(geom="text", x=15, y=0.33, label="Competition ($\\mathcal{C}$)")+
    theme(aspect.ratio = 0.666667,
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA),
          legend.key = element_rect(fill = NA,
                                    color = NA),
          legend.key.width = unit(1.2, 'cm'),
          legend.spacing.x = unit(0, 'cm'),
          legend.text = element_text(size = 5),
          legend.title = element_text(size = 25, hjust = 0.48),
          legend.box = 'vertical',
          legend.direction = 'horizontal',
          legend.position = 'bottom',
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          axis.text = element_text(size = 15))+
    scale_color_brewer(palette = "RdYlBu", name = "$l_P \\rightarrow l_T$", 
                       labels = c("$0.1 \\rightarrow 0.1$", 
                                  "$0.1 \\rightarrow 0.9$", 
                                  "$0.9 \\rightarrow 0.1$", 
                                  "$0.9 \\rightarrow 0.9$"))+
    # scale_fill_brewer(palette = "RdYlGn", name = "l2 --> l1", 
    #                    labels = c("0.1      0.1", "0.1      0.9", "0.9      0.1", "0.9      0.9"))+
    # scale_y_continuous(breaks = c(0.3, 0.8, 1.3))+
    guides(color = guide_legend(override.aes = list(size = 3),
                                nrow = 2,
                                byrow = TRUE,
                                title.position = 'top'))+
    labs(y = "Interaction strength",
         x = "")
   
  legend <- get_legend(interactions)
  
  interactions <- interactions + theme(legend.position="none")

depletion = 
  ggplot()+
  geom_line(size = 2, data = grouped_results, aes(x = n_coal, y = Rtav,
                                        color = id))+
  theme(aspect.ratio = 0.666667,
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.position = "none",
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 13),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 15))+
  scale_color_brewer(palette = "RdYlBu")+
  #scale_y_continuous(breaks = c(0.3, 0.8, 1.3))+
  guides(color = guide_legend(label.position = "left")) +
  labs(y = "Avg. resource concentratrion $(R^{\\star}_{tot})$",
       x = "")

invasions = 
ggplot()+
  geom_line(size = 2, linetype ='longdash', data = grouped_results, 
            aes(x = n_coal, y = ninv_av, color = id))+
  geom_line(size = 2, data = grouped_results, 
            aes(x = n_coal, y = r_av, color = id))+
  geom_segment(aes(x = 6, y = 35, xend = 8.5, yend = 35), size = 2)+
  geom_segment(aes(x = 6, y = 25, xend = 8.5, yend =25), linetype = 'longdash', size = 2)+
  annotate(geom="text", x=15, y=35, label="Species richness")+
  annotate(geom="text", x=16, y=25, label="Succesful invasions")+
  geom_abline(slope = 0, intercept = 60, linetype = 'dashed', size = 1.5)+
  geom_abline(slope = 0, intercept = 0, linetype = 'dashed', size = 1.5)+
  theme(aspect.ratio = 0.666667,
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.position = "none",
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text = element_text(size = 15))+
  scale_color_brewer(palette = "RdYlBu")+
  guides(color = guide_legend(label.position = "left")) +
  labs(y = "Number of species",
       x = "Number of coalescence exposures")+
  scale_y_continuous(limits = c(-1, 61))

cost = 
ggplot()+
  geom_line(size = 2, data = grouped_results, aes(x = n_coal, y = costav,
                                                  color = id))+
  theme(aspect.ratio = 0.666667,
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.position = "none",
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text = element_text(size = 15))+
  scale_color_brewer(palette = "RdYlBu")+
  guides(color = guide_legend(label.position = "left")) +
  labs(y = "Avg maintenace cost",
       x = "")

###############################################################################################################

#Plot abundance of each n-generalist species group along dynamics

rec_n_r_results = read.csv('../data/recursive_n_pref_results.csv')

grouped_results_nr = rec_n_r_results %>% group_by(n_pref, n_coal) %>%
  summarise(p_nr_av = mean(N_n_pref)) %>%
  filter(n_pref < 6)

n_pref_ab = 
  ggplot(data = grouped_results_nr) + 
  geom_line(aes(x = n_coal, 
                y = p_nr_av,
                color = as.factor(n_pref + 1)),
            size = 2) + 
  theme(aspect.ratio = 0.666667,
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.position = c(0.5, 0.5),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size = 15),
        legend.box = 'vertical',
        legend.direction = 'horizontal',
        legend.title = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text = element_text(size = 15))+
  scale_color_brewer(palette = "RdYlBu", direction = -1)+
  labs(y = "Proportional abundance",
       x = "",
       color = '$n_r$')



#Plot alltogether
lay <- rbind(c(NA, 3, NA, 1, 1, 1),
             c(5, 5, 5, 2, 2, 2),
             c(6, 6, 6, 4, 4, 4))
tikz("../sandbox/figure_5.tex",
     width = 8, height = 8.5,
     standAlone = TRUE)
grid.arrange(interactions, depletion, legend, invasions, cost, n_pref_ab,
             layout_matrix = lay)
dev.off()



