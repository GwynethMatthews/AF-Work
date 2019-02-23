import Atrium_Final as AF
import numpy as np
import pickle

import matplotlib.pyplot as plt
"""Defaults: L=200,v=0.5,d=0.05,e=0.05,rp=50,tot_time = 10**6
                 ,pace_rate = 220,seed1 = 1,seed2=2,seed3=3"""
def RiskCurve(repeats,nus,seeds):
    """Collects data for and plots the risk curve. Saves data."""
    AF_risks = []
    v_tAFs = []
    for j in range(len(nus)):
        print(j)
        v_tAF = []
        #propagates = 0
        for i in range(repeats):
            #print(i)
            A = AF.SourceSinkModel(hexagonal=True, threshold=1,rp=50, nu_para=j,nu_trans=j, 
                                   p_nonfire=0.05, pace_rate=220,seed_connections=seeds[j][i][0], seed_prop=seeds[j][i][1])
            A.cmp_full()
            
            v_tAF.extend([A.t_AF/A.tot_time])
            print(A.t_AF)
            if A.AF == True:
                print(A.AF)
            #propagates += sum(A.number_of_excitations[A.last_col])

        sample_avg = sum(v_tAF)

        v_tAFs.extend([v_tAF])
        AF_risks.extend([sample_avg])
    data = [nus,AF_risks,v_tAFs]
    #data = np.array(data)
    pickle.dump(data,open("RiskCurveTrialforPoster7.p","wb"))
    return data

repeats = 1
nus_trans = np.linspace(0.33, 1, 68, endpoint = True)
nus_para = 1
s = np.random.randint(0, 2**31, (80,1,2),dtype='int')
data = RiskCurve(repeats,nus_trans,s)

plt.scatter(data[0],data[1])
