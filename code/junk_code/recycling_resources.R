setwd("~/Desktop/coalescence_paper/code")
#Create a vector of arguments
args = commandArgs(trailingOnly=TRUE)
#Load data
if (length(args) >= 1){
  c_mats = read.csv(paste('../data/c_matrices_', args[1], '.csv', 
                               sep = ''))
  abundances = read.csv(paste('../data/abundances_', args[1], '.csv',
                              sep = ''))
  D_mats = read.csv(paste('../data/D_matrices_', args[1], '.csv', 
                          sep = ''))
} else{
    c_mats = read.csv('../data/c_matrices.csv')
    abundances = read.csv('../data/abundances.csv')
}

#Change column names
names(abundances) = c('X', as.character(seq(70)-1))
abundances['X'] = rownames(abundances)
#Melt the data
abundances_m = melt(data = abundances, id.vars = c('X'), variable.name = 'species', value.name = 'stable_state')
#Order data according to simulation number
abundances_m = abundances_m[order(abundances_m$X),]
#Merge abundances with c_mat by row
c_mats_ab= cbind(c_mats, abundances_m$stable_state)
colnames(c_mats_ab)[length(c_mats_ab)] = 'stable_state'
#Eliminate extinct species
c_mats_ab = c_mats_ab[c_mats_ab$stable_state > 1 ,]
#Get richness of each community
richness = table(c_mats_ab$X)
#Get indices of communities with only one member
ind_1 = as.numeric(names(richness[richness == 1]))
#Get rid of this in c_mats and D_mats
c_final = c_mats_ab[c_mats_ab$X != ind_1,]
