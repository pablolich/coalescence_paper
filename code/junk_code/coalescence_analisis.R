setwd("~/Desktop/coalescence_paper/code")
#Load libraries
require(ggplot2)
require(RColorBrewer)
require(viridis)
require(itertools)
library(ggpubr)
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
  coal_results_i = coal_results[#coal_results$kc == grid$kc[i] & 
                                coal_results$l1 == grid$leakage[i],]
  #x = coal_results_i$O1 - coal_results_i$O2
  x = coal_results_i$Fav1 + coal_results_i$Cav2 - coal_results_i$Fav2 - coal_results_i$Cav1 
  # x = coal_results_i$H1 - coal_results_i$H2
  # x = coal_results_i$r1 - coal_results_i$r2
  y = coal_results_i$S
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

print('Plotting')
if (length(args) >= 1){
  pdf(file = paste('../sandbox/correlations_richness', args[1], '.pdf',
                   sep = ''))
} else{
  pdf(file = '../sandbox/correlations_richness.pdf')
}

for (i in unique(df$leakage)){
    data = df[(df$leakage == i),]
    mytitle = paste('l = ', as.character(i))
    pl = ggplot(data, 
                aes(x = x, y = S))+
      geom_ribbon(aes(ymin = S - S_err, ymax = S + S_err), 
                  alpha = 0.5
                  #width = .05,
                  #size = 0.2,
      )+
      geom_point(size = 1.5)+
      # geom_errorbarh(aes(xmin = x - x_err, xmax = x + x_err),
      #                size = 0.2,
      #                height = 0.02,
      #                color = 'grey')+
      labs(x = expression(Theta),
           y = expression(S(C[1],~C[2])),
           title = mytitle)+
      theme(aspect.ratio = 1,
            legend.position = c(0.2, 0.8),
            panel.background = element_blank(), 
            panel.border = element_rect(color = 'black', fill=NA),
            axis.title = element_text(size = 20),
            axis.text = element_text(size = 15))+
      scale_x_continuous(expand = c(0,0))+
      scale_y_continuous(expand = c(0,0))
    print(pl)
    }
dev.off()

#Plot correlations barplot
if (length(args) >= 1){
  pdf(file = paste('../sandbox/correlation_barplot_richness', args[1], '.pdf'))
} else{
  pdf(file = '../sandbox/correlation_barplot_richness.pdf')
}
ggplot(data = grid, aes(#fill = as.factor(kc), 
                        x = as.factor(leakage), y = cor))+
  geom_bar(position = 'dodge', stat = 'identity')+
  scale_fill_viridis(discrete = TRUE)+
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA),
        panel.grid.major.x = element_line(color = 'grey'),
        aspect.ratio = 1,
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))+
  geom_hline(aes(yintercept = 0), data = grid)+
  labs(x = 'leakage',
       y = 'correlation',
       fill = 'Richness')
dev.off()

##################################################################################
#Heatmap of coalescence success as a function of kc and kf.
##################################################################################
#Initialize data frame
l_vec = unique(coal_results$l1)
kc_vec = unique(coal_results$kc1)
kf_vec = unique(coal_results$kf1)
success_map = expand.grid(l_vec, kc_vec, kf_vec)
colnames(success_map) = c('l', 'kc', 'kf')
success_map['n_success'] = 0
#Get values of leakage
l_values = unique(coal_results$l1)
for (l in l_values){
  #Select coalescence events only with that leakage
  coal_results_l = coal_results[coal_results$l1 == l,]
  #Get number of coalescence events
  n_coal = length(coal_results_l$kc1)
  for (i in seq(n_coal)){
    #Get similarity of the outcome community
    S = coal_results_l$S[i]
    #Store constants of the winning community
    if (S > 0){
      #Winning community is 1
      kc = coal_results_l$kc1[i]
      kf = coal_results_l$kf1[i]
      #Find that kc-kf coordinate in the map
      point = as.numeric(rownames(success_map[(success_map$l == l &
                                                 success_map$kc == kc &
                                                 success_map$kf == kf),]))
      #Add 1 to the n_succes of that point
      success_map$n_success[point] = success_map$n_success[point] + 1
    }else if (S < 0){
      #Winning community is 2
      kc = coal_results_l$kc2[i]
      kf = coal_results_l$kf2[i]
      point = as.numeric(rownames(success_map[(success_map$l == l &
                                                 success_map$kc == kc &
                                                 success_map$kf == kf),]))
      #Add 1 to the n_succes of that point
      success_map$n_success[point] = success_map$n_success[point] + 1
    }else{
      #Both communities win (or they loose the same percent of species)
      kc1 = coal_results_l$kc1[i]
      kf1 = coal_results_l$kf1[i]
      kc2 = coal_results_l$kc2[i]
      kf2 = coal_results_l$kf2[i]
      point1 = as.numeric(rownames(success_map[(success_map$l == l &
                                                 success_map$kc == kc1 &
                                                 success_map$kf == kf1),]))
      point2 = as.numeric(rownames(success_map[(success_map$l == l &
                                                  success_map$kc == kc2 &
                                                  success_map$kf == kf2),]))
      success_map$n_success[point1] = success_map$n_success[point1] + 1
      success_map$n_success[point2] = success_map$n_success[point2] + 1
      next
    }
  }
}

#Plot this in the heat map
pdf('../sandbox/kc_kf_success.pdf', width = 12, height = 15)
plots = list()
plots[]
counter = 1
s_rng = range(success_map$n_success)
s_mid = (s_rng[2]-s_rng[1])/2 + s_rng[1]
for (i in l_vec){
  success_l = success_map[success_map$l == i,]
  p = ggplot(data = success_l,
            aes(x = kc, y = kf,
                fill = n_success))+
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
    #                      midpoint=s_mid,    #same midpoint for plots (mean of the range)
    #                      breaks=seq(200,600,200), #breaks in the scale bar
    #                      limits=c(floor(s_rng[1]-1), ceiling(s_rng[2]))) + #same limits for plots
    scale_fill_distiller(palette = 'RdBu')+
    labs(title  = paste('l = ', i),
         fill = '+',
         x = expression(k[c]),
         y = expression(k[f]))
  plots[[counter]] = p
  counter = counter + 1
}
#Store plots in nice figure
ggarrange(plots[[1]], plots[[2]], plots[[3]],
          # plots[[4]], plots[[5]], plots[[6]],
          # plots[[7]], plots[[8]], plots[[9]],
          # plots[[10]], plots[[11]], plots[[12]],
          # plots[[13]], plots[[14]], plots[[15]],
          # plots[[16]], plots[[17]], plots[[18]],
          # plots[[19]], plots[[20]],
          ncol = 3, nrow = 1,  
          label.x = 0, label.y = 0.91,
          font.label = list(size = 20))
dev.off()



#######################################################################
#Plot all the average abundance heatmaps
#######################################################################
pdf('../sandbox/kc_kf_std_demand.pdf', width = 13, height = 15)
l = unique(matrix_richness_clean$l)
plots = list()
plots[]
counter = 1
r_rng = range(matrix_richness_clean$std_dem)
r_mid = (r_rng[2]-r_rng[1])/2 + r_rng[1]
for (i in l){
  data = matrix_richness_clean[matrix_richness_clean$l == i,]
  p = ggplot(data = data,
             aes(x = data$kc, y = data$kf,
                 fill = std_dem))+
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
    # scale_fill_gradient2(high="#2166AC", mid="#F7F7F7", low="#B2182B", #colors in the scale
    #                      midpoint=r_mid,    #same midpoint for plots (mean of the range)
    #                      breaks=seq(30,50,10), #breaks in the scale bar
    #                      limits=c(floor(r_rng[1]), ceiling(r_rng[2]))) + #same limits for plots
    labs(title  = paste('l = ', i),
         fill = 'dem',
         x = expression(k[c]),
         y = expression(k[f]))
  plots[[counter]] = p
  counter = counter + 1
}
#Store plots in nice figure
ggarrange(plots[[1]], plots[[2]], plots[[3]], 
          plots[[4]], plots[[5]], plots[[6]],
          plots[[7]], plots[[8]], plots[[9]],
          plots[[10]],
          ncol = 3, nrow = 3,
          #labels = c('A', 'B', 'C', 'D', 'E', 'F','G', 'H', 'I', 'J'),  
          label.x = 0, label.y = 0.91,
          font.label = list(size = 20))
dev.off()

################################################################################################

#Calculate average (kc, kf) weighted by number of success for each l
#Initialize storing  matrix
points = matrix(0, nrow = length(l_vec), ncol = 5)
l_vec = unique(success_map$l)
for (i in seq(length(l_vec))){
  data_l = success_map[success_map$l == l_vec[i],]
  mean_kc = sum(data_l$kc*data_l$n_success)/sum(data_l$n_success)
  sem_kc = std(rep(data_l$kc, times = data_l$n_success))/sqrt(sum(data_l$n_success))
  mean_kf = sum(data_l$kf*data_l$n_success)/sum(data_l$n_success)
  sem_kf = std(rep(data_l$kf, times = data_l$n_success))/sqrt(sum(data_l$n_success))
  
  points[i,1] = l_vec[i]
  points[i,2] = mean_kc
  points[i,3] = sem_kc
  points[i,4] = mean_kf
  points[i,5] = sem_kf
}
colnames(points) = c('l', 'kc', 'kc_err', 'kf', 'kf_err')
data = data.frame(points)
ggplot(data, aes(x = kc,  y = kf ))+
  geom_errorbar(aes(ymin = kf - kf_err, ymax = kf + kf_err), 
                width = .001,
                size = 0.2,
                color = 'grey')+
  geom_errorbarh(aes(xmin = kc - kc_err, xmax = kc + kc_err),
                 size = 0.2,
                 height = 0.002,
                 color = 'grey')+
  geom_point(aes(colour = l))+
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA),)+
  scale_color_viridis()

################################################################################################
#Recursive coalescence
################################################################################################
if (length(args) >= 1){
  recursive_results= read.csv(paste('../data/recursive_coalescence_results', args[1], '.csv', 
                                sep = ''))
} else{
  recursive_results = read.csv('../data/coalescence_results.csv')
}

recursive_results = recursive_results[(recursive_results$l != 0.01) & 
                                        (recursive_results$l != 0.5),]

#Specify wether you are cohesive or codependant
recursive_results['role'] = 0
for (i in seq(length(recursive_results$l))){
  if (recursive_results$l[i] > 0.5){
    recursive_results$role[i] = 1
  }
}
# recursive_results = recursive_results[recursive_results$Ftot != 0,]

fc = facilitation_cycle(recursive_results$l)
cc = competition_cycle(recursive_results$l)
recursive_results$cfac = cc
recursive_results$ffac = fc
#recursive_results$Ftot = recursive_results$Ftot/fc
recursive_results$coh = recursive_results$Fav-#/fc - 
                        recursive_results$Cav#/cc
mean_rec = aggregate(coh~l+n_coal+role,
                     recursive_results,
                     FUN = mean)

std_rec = aggregate(coh~l+n_coal+role,
                    recursive_results,
                    FUN = sd)

mean_rec['err'] = std_rec$coh/sqrt(100)

#Check competition
mean_comp = aggregate(Cav ~ n_coal + l,
                     recursive_results,
                     FUN = mean)
#Check competition
mean_fcomp = aggregate(cfac ~ n_coal + l,
                      recursive_results,
                      FUN = mean)
mean_comp$cfac = mean_fcomp$cfac

#Check facilitation
mean_faci = aggregate(Fav ~ n_coal + l,
                      recursive_results,
                      FUN = mean)
#Check faciliitation
mean_ffac = aggregate(ffac ~ n_coal + l,
                     recursive_results,
                     FUN = mean)
mean_faci$ffac = mean_ffac$ffac

#Check coheson
mean_faci$coh = mean_comp$Cav - mean_faci$Fav
#Plot
ggplot(data = mean_faci)+
  geom_line(aes(x = n_coal, y = coh,
                group = as.factor(l),
                color = l),
            size = 0.7)+
  # geom_ribbon(data = mean_rec,
  #             aes(x = n_coal,
  #                 ymin = Ftot - err, ymax = Ftot + err,
  #                 group = l,
  #                 fill = l),
  #             alpha = 0.3)+
   scale_color_viridis()+
  # scale_fill_viridis()+
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA))

ggplot(data = recursive_results)+
  geom_line(aes(x = n_coal, y = Ftot,
                group = l + n_sim, color = as.factor(l + n_sim)),
            alpha = 0.2)

#################################################################################################
#check initial values of F and C as a function of leakage
data_f =  mean_faci[mean_faci$n_coal == 0,]
data_c = mean_comp[mean_comp$n_coal == 0,]

plot(data_c$l, data_c$Cav, pch = 20)
lines(data_c$l, comp(data_c$l))
points(data_f$l, data_f$Fav)
lines(data_f$l, fac(data_c$l))
#################################################################################################
#Plot total richness of communities as function of time

mean_R = aggregate(R~l+n_coal,
                     recursive_results,
                     FUN = mean)
std_R = aggregate(R~l+n_coal,
                    recursive_results,
                    FUN = sd)
mean_R['err'] = std_R$R
ggplot(data = recursive_results)+
  geom_line(aes(x = n_coal, y = R, 
                group = l+n_sim,
                color = as.factor(l)))+
  # geom_ribbon(aes(x = n_coal,
  #                 ymin = R - err, ymax = R + err,
  #                 group = l,
  #                 fill = as.factor(l)),
  #             alpha = 0.3)+
  theme(aspect.ratio = 1)

#################################################################################################
#Plot species richness
mean_r = aggregate(r~l+n_coal,
                   recursive_results,
                   FUN = mean)
std_r = aggregate(r~l+n_coal,
                  recursive_results,
                  FUN = sd)
mean_r['err'] = std_r$r
ggplot(data = recursive_results)+
  geom_line(aes(x = n_coal, y = r, 
                group = l + n_sim,
                color = as.factor(l)))+
  # geom_ribbon(aes(x = n_coal,
  #                 ymin = r - err, ymax = r + err,
  #                 group = l,
  #                 fill = as.factor(l)),
  #             alpha = 0.3)+
  theme(aspect.ratio = 1)
#################################################################################################
#Plot number of extinctions

mean_ext = aggregate(ext~l+n_coal,
                   recursive_results,
                   FUN = mean)
mean_inv = aggregate(inv~l+n_coal,
                     recursive_results,
                     FUN = mean)
std_ext = aggregate(ext~l+n_coal,
                  recursive_results,
                  FUN = sd)
mean_ext['err'] = std_ext$ext
mean_ext['inv'] = mean_inv$inv
ggplot(data = mean_ext)+
  geom_ribbon(aes(x = n_coal,
                  ymin = ext - err, ymax = ext + err,
                  group = l,
                  fill = as.factor(l)),
              alpha = 0.3)+
  geom_line(aes(x = n_coal, y = ext, 
                group = l,
                color = as.factor(l)),
            size = 1.3)+
  geom_line(aes(x = n_coal, y = inv,
                group = l,
                col = as.factor(l)),
            linetype = 'dashed')+
  theme(aspect.ratio = 1)


#################################################################################################
#Plot average number of pathways
mean_path = aggregate(path ~ n_coal+l, 
                      recursive_results,
                      FUN = mean)

ggplot(data = mean_path)+
  geom_line(aes(x = n_coal, 
                y = path,
                group = n_sim + l,
                color = as.factor(l)))

ggplot(data = recursive_results) +
  geom_line(aes(x  = n_coal,
                y = sd_met,
                group = n_sim + l, 
                color = as.factor(l)))



#################################################################################################
#Check results
if (length(args) >= 1){
  recursive_diff_results = read.csv(paste('../data/recursive_diff_coalescence_results', args[1], '.csv', 
                                    sep = ''))
} else{
  recursive_diff_results = read.csv('../data/recursive_diff_coalescence_results.csv')
}

#Plot cohesion-similarity plot

coh = recursive_diff_results$F1 - recursive_diff_results$C1 - 
  recursive_diff_results$F2 + recursive_diff_results$C2


plot(recursive_diff_results$n_coal, )
plot(recursive_diff_results$n_coal, recursive_diff_results$Fav - recursive_diff_results$Cav)
