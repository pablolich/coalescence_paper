setwd("~/Desktop/coalescence_paper/code")
#Load libraries
require(ggplot2)
require(RColorBrewer)
require(viridis)
#Load data
simulations = read.csv('../data/simulation_results.csv')
#Sample rows from this huge dataset
rows = sample(1:nrow(simulations), 1e4)
simulations = simulations[rows,]
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
p = ggplot(data = sim, aes(x = C0, y = F0))+
  geom_point(aes(fill = kc + kf,
                 shape = as.factor(l)), 
             colour = 'grey',
             stroke = 0.1)+
  theme(aspect.ratio = 1)+
  geom_abline(slope = 0.1/0.9, intercept = 0, 
              colour = brewer.pal(n = 5, name = "RdYlBu")[1],
              size = 1)+
  geom_abline(slope = 0.3/0.7, intercept = 0, 
              colour = brewer.pal(n = 5, name = "RdYlBu")[2],
              size = 1)+
  geom_abline(slope = 1, intercept = 0, 
              colour = brewer.pal(n = 5, name = "RdYlBu")[3],
              size = 1)+
  geom_abline(slope = 0.7/0.3, intercept = 0, 
              colour = brewer.pal(n = 5, name = "RdYlBu")[4],
              size = 1)+
  geom_abline(slope = 0.9/0.1, intercept = 0, 
              colour = brewer.pal(n = 5, name = "RdYlBu")[5],
              size = 1)+
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = 'grey'),
        panel.border  = element_rect(colour = 'black', fill = NA),
        legend.title = element_text(size = 15),
        legend.title.align = 0.5,
        legend.key = element_blank())+
  annotate("text", x = 1.6, y = 2.5, 
           label = expression(paste(m[j], '=', frac(l[j], paste('(',1-l[j],')')))),
           size = 5)+
  labs(x = expression(paste(C[0])),
       y = expression(paste(F[0])),
       shape = expression(paste(l[j])),
       fill = expression(paste(k[c] + k[f])))+
  scale_shape_manual(values = c(21, 22, 23, 24, 25))+
  scale_fill_viridis()
#Show
p


##Fit data from sim to an exponentiexal

fit <- nls(sim$C0~I(exp(x*a)), data=sim, start=list(a = 2))
exponential = function(c0, a, x){
  return(c0*exp(a*x))
}