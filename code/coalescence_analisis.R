setwd("~/Desktop/coalescence_paper/code")
#Load libraries
require(ggplot2)
require(RColorBrewer)
require(viridis)
source('functions_coalescence.R')
#Load data
coal_results = read.csv('../data/coalescence_results.csv')
#Set parameter grid of analysis
vec_rich = unique(coal_results$r1)
vec_leak = unique(coal_results$l1)
vec_points = seq(20)
#Create grid of parameters
grid = expand.grid(list(richness = vec_rich, 
                        leakage = vec_leak))
#Create data frame with that grid
df = data.frame(expand.grid(list(point = vec_points, 
                                 richness = vec_rich, 
                                 leakage = vec_leak)))
#Add columns to dataframe for plotting purposes
df['x'] = 0
df['S'] = 0
df['x_err'] = 0
df['S_err'] = 0
#Get number of total simulations
n_sim = length(grid$richness)
#Get number of breaks
n_mids = length(vec_points) 
#Iterate over it and generate plotting points
for (i in seq(n_sim)){
  coal_results_i = coal_results[coal_results$r1 == grid$richness[i] & 
                                coal_results$l1 == grid$leakage[i],]
  x = coal_results_i$O2 - coal_results_i$O1
  y = coal_results_i$S
  group = bin_data(x, y, n_mids + 1, centralize = 'median')
  df[((i-1)*n_mids+1):(i*n_mids),]$x = group$x_binned
  df[((i-1)*n_mids+1):(i*n_mids),]$S = group$y_binned
  df[((i-1)*n_mids+1):(i*n_mids),]$x_err = group$sigma_x
  df[((i-1)*n_mids+1):(i*n_mids),]$S_err = group$sigma_y
}

#Plot

for (i in unique(df$leakage)){
  for (j in unique(df$richness)){
    data = df[(df$leakage == i) &
              (df$richness == j),]
    mytitle = paste('l = ', as.character(i), ', ',
                    'r = ', as.character(j))
    pl = ggplot(data = data, 
                aes(x = x, y = S,
                    color = as.factor(richness)))+
      geom_point()+
      geom_errorbar(aes(ymin = S-S_err, ymax=S+S_err), 
                    width=.2,
                    position=position_dodge(0.05))+
      geom_errorbarh(aes(xmin = x - x_err ,xmax = x + x_err),
                     width=.2,
                     position=position_dodge(0.05))+
      labs(title = mytitle )+
      theme(aspect.ratio = 1,
            legend.position = 'none')
    print(pl)
  }
}

