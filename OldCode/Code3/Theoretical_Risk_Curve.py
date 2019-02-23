"""Risk curve with theoretical risk curve, Kishan's data and my data"""
import Atrium_class as AC
import numpy as np
import matplotlib.pyplot as plt
import pickle
import scipy.stats as sp
A = AC.Atrium()

L = A.size
delta = A.dysfunctional_prob
tau = A.rp
nus = np.array([0.02, 0.04, 0.06, 0.08, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3,0.5,1])
T_risk = (1 - ((1-(1-nus)**tau)**(delta*L*L)))
results = pickle.load(open("Risk_Curve_keep.p","rb"))
#results = pickle.load(open("data_risk_curve_compare.p","rb"))
#print(results[2])
K_results = pickle.load(open("Kishan_data_risk_curve.p","rb"))
std_values = []
for i in range(len(results[2])):
    std_values.extend([np.std(results[2][i])/np.sqrt(len(results[2][i]))])
#plt.scatter(np.array(results[0]),np.array(results[1]),c='r',marker = 'x',label = 'My Data')
#plt.violinplot(results[2],results[0],showmeans = False,showextrema = False,widths = 0.005,points = 20)
plt.errorbar(np.array(results[0]),np.array(results[1]),yerr = std_values,fmt='x',c = 'r',label = 'My Data')
plt.plot(nus,T_risk, 'g', label = 'Theoretical Data')
#plt.scatter(K_results[0],K_results[1],c='b',marker = 'x', label = "Kishan's Data" )
plt.errorbar(np.array(K_results[0]),np.array(K_results[1]),yerr = K_results[2],fmt='x',c = 'b',label = 'Kishan Data')
plt.legend(fontsize = 20)
plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=20)
plt.title('Risk of AF', fontsize = 24)
plt.xlabel('Fraction of Transverse Connections',fontsize = 20)
plt.ylabel('Time in AF/Risk of AF',fontsize = 20)
plt.show()
ktest1 = sp.ks_2samp(results[1][1:-2],K_results[1])
ktest2 = sp.ks_2samp(results[1][1:-2],T_risk)
print(ktest1)
print(ktest2)
