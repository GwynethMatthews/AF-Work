import numpy as np 
import matplotlib.pyplot as plt

results1 = np.load('Averages_tau_50_charge_storage_on_pr_60.npy')
results2 = np.load('Averages_tau_70_charge_storage_on_pr_80.npy')
results3 = np.load('Averages_tau_90_charge_storage_on_pr_100.npy')
results4 = np.load('Averages_tau_110_charge_storage_on_pr_120.npy')
results5 = np.load('Averages_tau_130_charge_storage_on_pr_140.npy')
#nu, tau, p , if charge stored or not, pace_rate,
# prob of termination if enters AF, 
# probability of entering AF,
# average time in AF (includes values for when AF in not entered i.e. t_AF = 0),
# average time of terminationg for repeats that terminate after entering AF
nu = np.linspace(0.325, 1, 28, endpoint = True)
p = np.array([0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15,0.2,0.25,0.30,0.4,0.5,0.6,0.7,0.8,1])

AF_probs1 = []
AF_probs2 = []
AF_probs3 = []
AF_probs4 = []
AF_probs5 = []

for i in results1:
    a = i[8]-((10*i[4])+250)
    #print(i[8])
    #print(i[4])
    #print(a)
    if a < 0:
        a = 0
    AF_probs1.extend([a]) 
for i in results2:
    a = i[8]-((10*i[4])+250)
    if a < 0:
        a = 0
    AF_probs2.extend([a]) 
for i in results3:
    a = i[8]-((10*i[4])+250)
    if a < 0:
        a = 0
    AF_probs3.extend([a])   
for i in results4:
    a = i[8]-((10*i[4])+250)
    if a < 0:
        a = 0
    AF_probs4.extend([a]) 
for i in results5:
    a = i[8]-((10*i[4])+250)
    if a < 0:
        a = 0
    AF_probs5.extend([a]) 
print(min(AF_probs1))
print(min(AF_probs2))
print(min(AF_probs3))
print(min(AF_probs4))
print(min(AF_probs5))

ax1 = plt.subplot2grid(shape=(5,8), loc=(0,0), colspan=2, rowspan=2)
ax2 = plt.subplot2grid((5,8), (0,3), colspan=2, rowspan=2)
ax3 = plt.subplot2grid((5,8), (0,6), colspan=2, rowspan=2)
ax4 = plt.subplot2grid((5,8), (3,1), colspan=2, rowspan=2)
ax5 = plt.subplot2grid((5,8), (3,5), colspan=2, rowspan=2)

mesh1 = ax1.pcolor(p,nu, np.array(AF_probs1).reshape(27,20))
mesh2 = ax2.pcolor(p,nu, np.array(AF_probs2).reshape(27,20))
mesh3 = ax3.pcolor(p,nu, np.array(AF_probs3).reshape(27,20))
mesh4 = ax4.pcolor(p,nu, np.array(AF_probs4).reshape(27,20))
mesh5 = ax5.pcolor(p,nu, np.array(AF_probs5).reshape(27,20))

#mesh1.set_clim(0,1)
#mesh2.set_clim(0,1)
#mesh3.set_clim(0,1)
#mesh4.set_clim(0,1)
#mesh5.set_clim(0,1)

cbar1 = plt.colorbar(mesh1, ax = ax1)
cbar2 = plt.colorbar(mesh2, ax = ax2)
cbar3 = plt.colorbar(mesh3, ax = ax3)
cbar4 = plt.colorbar(mesh4, ax = ax4)
cbar5 = plt.colorbar(mesh5, ax = ax5)

#ax1.colorbar()
#plt.pcolor(p,nu, np.array(Af_probs).reshape(27,20)[:,1:])
ax1.set_title('tau = 50', fontsize = 20)
ax1.set_xlabel('p',fontsize = 20)
ax1.set_ylabel('nu',fontsize = 20)
ax1.tick_params(axis='both', which='major', labelsize=16)
ax2.set_title('tau = 70', fontsize = 20)
ax2.set_xlabel('p',fontsize = 20)
ax2.set_ylabel('nu',fontsize = 20)
ax2.tick_params(axis='both', which='major', labelsize=16)
ax3.set_title('tau = 90', fontsize = 20)
ax3.set_xlabel('p',fontsize = 20)
ax3.set_ylabel('nu',fontsize = 20)
ax3.tick_params(axis='both', which='major', labelsize=16)
ax4.set_title('tau = 110', fontsize = 20)
ax4.set_xlabel('p',fontsize = 20)
ax4.set_ylabel('nu',fontsize = 20)
ax4.tick_params(axis='both', which='major', labelsize=16)
ax5.set_title('tau = 130', fontsize = 20)
ax5.set_xlabel('p',fontsize = 20)
ax5.set_ylabel('nu',fontsize = 20)
ax5.tick_params(axis='both', which='major', labelsize=16)