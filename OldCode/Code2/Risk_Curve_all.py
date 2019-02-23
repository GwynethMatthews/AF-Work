"""Theoretical Risk Curve"""
import Atrium_class as AC
import numpy as np
import matplotlib.pyplot as plt
import pickle
A = AC.Atrium()

L = A.size
delta = A.dysfunctional_prob
tau = A.rp
nus = np.array([0,0.05,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.5,1])
T_risk = (1 - ((1-(1-nus)**tau)**(delta*L*L)))
results = pickle.load(open("data_risk_curve1.p","rb"))
plt.scatter(np.array(results[0]),np.array(results[1]))
plt.plot(nus,T_risk)
nu_K = [0.02,0.04,0.06,0.08,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3]
risk_K = [0.99981,0.99983,0.9998,0.99968,0.99919,0.99772,0.99152,0.96099,0.86184,0.60984,0.29714,0.16381,0.039206,0.017807,0.0056277,0.020737,4.83E-05,4.922E-05,0.00082172,0.0001084,0,0,0,0,9.406E-05]
plt.scatter(nu_K,risk_K)
plt.show()