#Plot population and resource concentration time series
#import matplotlib.pylab as plt
#for sp in range(len(z0)):
#    if sp<params['s']:
#        plt.plot(sol.t, sol.y[sp])   
#    else:
#        plt.plot(sol.t, sol.y[sp], linestyle = '--')
#plt.xscale('log')
#plt.show()


#Integrate equations
#Add these to the parameter dictionary
#params['D'] = D
#params['c'] = c
#params['x'] = maint_vec.reshape(s,1)
#params['l'] = df['l'][i]*np.ones(m).reshape(m,1)

#Solve diferential equations
#sol = solve_ivp(lambda t,z: equations(t,z, params),
#                tspan, z0,
#                method = 'BDF', atol = 0.0001 )


#Plog results in python
#plt.scatter(df['C'], df['F'], c = df['kc'], edgecolor = 'k')
#clb = plt.colorbar()
#clb.ax.set_title(r'$k_c$')
#plt.xlabel(r'$C$', fontsize = 20)
#plt.ylabel(r'$F$', fontsize = 20)
#plt.legend()
#plt.show()


##Plot population and resource concentration time series after coalescence
#import matplotlib.pylab as plt
#for sp in range(len(z0)):
#    if sp<params['s']:
#        if sp<s1:
#            plt.plot(sol.t, sol.y[sp], color = 'blue')   
#        else:
#            plt.plot(sol.t, sol.y[sp], color = 'red')   
#            
#    else:
#        plt.plot(sol.t, sol.y[sp], linestyle = '--')
#plt.xscale('log')
#plt.show()
