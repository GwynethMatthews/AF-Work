import Atrium_class as AC
import numpy as np
import matplotlib.pyplot as plt
import pickle
import scipy.stats as sp
results1 = pickle.load(open("risk_curve_0.02-0.1.p","rb"))
results2 = pickle.load(open("risk_curve_0.11-0.15.p","rb"))
results3 = pickle.load(open("risk_curve_0.16-0.20.p","rb"))
results4 = pickle.load(open("risk_curve_0.21-0.25.p","rb"))
results5 = pickle.load(open("risk_curve_0.26-0.30.p","rb"))
results6 = pickle.load(open("risk_curve_0.16.p","rb"))
results7 = pickle.load(open("risk_curve_0.17.p","rb"))
results8 = pickle.load(open("risk_curve_0.18.p","rb"))
results9 = pickle.load(open( "risk_curve_0.16-0.19.p", "rb" ) )
K_results = pickle.load(open("Kishan_data_risk_curve.p","rb"))
L = 200
delta = 0.05
tau = 50
nus = np.array([0.02, 0.04, 0.06, 0.08, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3,0.5,1])
T_risk = (1 - ((1-(1-nus)**tau)**(delta*L*L)))
results_nu = np.concatenate((results1[0],results2[0],results9[0],results3[0][-1:],results4[0],results5[0]))
results_risk = np.concatenate((results1[1],results2[1],results9[1],results3[1][-1:],results4[1],results5[1]))
results_repeats = np.concatenate((results1[2],results2[2],results9[2],results3[2][-1:],results4[2],results5[2]))
plt.scatter(results_nu,results_risk,c='r',marker = 'x', label = "My Data")
plt.scatter(K_results[0],K_results[1],c='b',marker = 'x', label = "Kishan's Data" )
plt.plot(nus,T_risk, 'g', label = 'Theoretical Data')
plt.show()
print(results_repeats)
data = [results_nu,results_risk,results_repeats]
pickle.dump(data,open( "Risk_Curve_keep.p", "wb" ) )
