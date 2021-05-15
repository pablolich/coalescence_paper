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