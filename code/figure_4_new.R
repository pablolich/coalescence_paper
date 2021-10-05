setwd("~/Desktop/coalescence_paper/code")
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(tikzDevice)
library(RColorBrewer)
coal_rep = read.csv('../data/effect_facilitation_new_cost.csv', row.names = 1 )
coal_rep['id'] = paste(coal_rep$l_a, coal_rep$l_b)

l_B = c(0.1       , 0.10808081, 0.11616162, 0.12424242, 0.13232323,
0.14040404, 0.14848485, 0.15656566, 0.16464646, 0.17272727,
0.18080808, 0.18888889, 0.1969697 , 0.20505051, 0.21313131,
0.22121212, 0.22929293, 0.23737374, 0.24545455, 0.25353535,
0.26161616, 0.26969697, 0.27777778, 0.28585859, 0.29393939,
0.3020202 , 0.31010101, 0.31818182, 0.32626263, 0.33434343,
0.34242424, 0.35050505, 0.35858586, 0.36666667, 0.37474747,
0.38282828, 0.39090909, 0.3989899 , 0.40707071, 0.41515152,
0.42323232, 0.43131313, 0.43939394, 0.44747475, 0.45555556,
0.46363636, 0.47171717, 0.47979798, 0.48787879, 0.4959596 ,
0.5040404 , 0.51212121, 0.52020202, 0.52828283, 0.53636364,
0.54444444, 0.55252525, 0.56060606, 0.56868687, 0.57676768,
0.58484848, 0.59292929, 0.6010101 , 0.60909091, 0.61717172,
0.62525253, 0.63333333, 0.64141414, 0.64949495, 0.65757576,
0.66565657, 0.67373737, 0.68181818, 0.68989899, 0.6979798 ,
0.70606061, 0.71414141, 0.72222222, 0.73030303, 0.73838384,
0.74646465, 0.75454545, 0.76262626, 0.77070707, 0.77878788,
0.78686869, 0.79494949, 0.8030303 , 0.81111111, 0.81919192,
0.82727273, 0.83535354, 0.84343434, 0.85151515, 0.85959596,
0.86767677, 0.87575758, 0.88383838, 0.89191919, 0.9)    
l_A = c(0.1, 0.5, 0.9)
coal_rep['la'] = l_A[coal_rep$l_a+1]
coal_rep['lb'] = l_B[coal_rep$l_b+1]


coal_rep_grouped = coal_rep %>% group_by(la, lb) %>% summarise(S_av = mean(S), S_sd = sd(S),
                                                               F_av = mean(Fav), F_av_sd = sd(Fav)) #%>%
  #filter(la != 0.5)
#Get coefficients of transformation from l to F
coeff = lm(coal_rep_grouped$lb~coal_rep_grouped$F_av)
rep_coal = 
  ggplot(data = coal_rep_grouped,
       aes(x = F_av, y = S_av))+
  geom_ribbon(aes(x = F_av, ymin = S_av - S_sd, ymax = S_av + S_sd,
                  fill = as.factor(la)), alpha = 0.3)+
  geom_line(size = 5.5, aes(group = la,
                colour = as.factor(la)))+
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.position = c(0.20, 0.85),
        legend.key.width = unit(1, 'cm'),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20))+
  scale_color_brewer(palette = 'Set1',
                     name = '$l_A$')+
  scale_fill_brewer( palette = 'Set1',
                    name = '$l_A$')+
  scale_y_continuous(breaks = c(-0.5, 0, 0.5))+
  scale_x_continuous(expand = c(0, 0), 
                     sec.axis = sec_axis(~ .*coeff$coefficients[2] + coeff$coefficients[1],
                                         name = 'Leakage of community B $(l_B)$',
                                         breaks = c(0.2, 0.5, 0.8)))+
  labs(x = "Facilitation of community $B$ $(\\mathcal{F}_B)$",
       y = "Similarity to parents $(S_{A, B})$")
  tikz("../sandbox/figure_4_n.tex",
     width = 8, height = 8.5,
     standAlone = TRUE)
print(rep_coal)
dev.off()

#Plot autocorrelation
autocorr_matrix = read.csv('../data/autocorr.csv', row.names = 1)

