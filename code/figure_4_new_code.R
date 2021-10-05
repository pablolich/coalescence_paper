setwd("~/Desktop/coalescence_paper/code")
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(gridExtra)
dev.off()
coal_rec = read.csv('../data/recursive_results.csv',row.names = 1 )
#Create unique id for each group
coal_rec['id'] = paste(as.character(coal_rec$l1), as.character(coal_rec$l2), 
                       sep = '')
div_rows = which(coal_rec$Rav == 0)
coal_rec = coal_rec[-div_rows,]
#Number of groups
groups = unique(coal_rec$id)
n_groups = length(groups)
#Group by id
grouped_results <- coal_rec %>% group_by(n_coal, id) %>% summarise(Ctav = mean(Cav_i), Ctsd = sd(Cav_i),
                                           Ftav = mean(Fav_i), Ftsd = sd(Fav_i),
                                           Rtav = mean(Rav), Rtsd = sd(Rtav),
                                           ninv_av = mean(n_inv), ninv_sd = sd(n_inv),
                                           ltav = mean(lav)) %>% arrange(id, n_coal)

interactions = ggplot()+
  geom_line(data = grouped_results, aes(x = n_coal, y = Ftav/ltav,
                                        color = id))+
  geom_line(data = grouped_results , aes(x = n_coal, y =Ctav,
                                         color = id),
            linetype = 'longdash')

depletion = ggplot()+
  geom_line(data = grouped_results, aes(x = n_coal, y = Rtav,
                                        color = id))

invasions = ggplot()+
  geom_line(data = grouped_results, aes(x = n_coal, y = ninv_av,
                                        color = id))


#Plot alltogether
lay <- rbind(c(1, 2, 2, 3, 4),
             c(6, 2, 2, 5, 7))
options( tikzLatexPackages = c( getOption( "tikzLatexPackages" ), "\\usepackage{amssymb}"))
tikz("../sandbox/figure_3.tex",
     width = 15.33, height = 5.99,
     standAlone = TRUE)
grid.arrange(Deff_micro, similarity_micro, depletion_05, depletion_09, 
             depletion_01, C_F_pref, legend,
             layout_matrix = lay)
dev.off()
