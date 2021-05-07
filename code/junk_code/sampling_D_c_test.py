#Testing different cases of the sampling method

from model import *
from functions import *

s = 10; m = 5
c = preference_matrix(m, s, 0.99) 
cum_c = np.cumsum(c, axis = 0)
demands = cum_c[-1,:]
D = crossfeeding_matrix(demands, m, 0.9999)
#Try the cohesion calculator
l = 0.5*np.ones(m)
cohesion_matrix(l, c, D)
D_ = crossfeeding_matrix(demands, m, -0.9999)
crossfeeding_vector = np.sum(D, axis = 0)
cross_ = np.sum(D_, axis = 0)
#Plot
#Get vectors of transparencies
alph_exp = np.exp(5*np.arange(s)/s)
alph = (alph_exp-min(alph_exp))/max(alph_exp-min(alph_exp))
col = 'black'
for i in range(s):
    if i == s-1:
        col = 'red'
        lw = 2
        plt.plot(np.arange(m), cum_c[i,:]*sum(cross_)/sum(demands), 
                 alpha = alph[i], color = col, 
                 linewidth = lw)
    else:
        plt.plot(np.arange(m), cum_c[i,:]*sum(cross_)/sum(demands),
                 alpha = alph[i], color = col)

plt.plot(np.arange(m), crossfeeding_vector, color = 'limegreen', 
         label = r'$k_f \rightarrow 1$')
plt.plot(np.arange(m), cross_, color = 'darkgoldenrod', 
         label = r'$k_f \rightarrow -1$')
plt.title(r'$k_c \rightarrow 1$', fontsize = 25)
plt.legend(fontsize = 15)
plt.xlabel('Metabolite', fontsize = 20)
plt.ylabel(r'$d_i$', fontsize = 20)
plt.show()
