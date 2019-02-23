import Atrium2 as AC
"""Defaults: L=200,v=0.5,d=0.05,e=0.05,rp=50,tot_time = 10**6
                 ,pace_rate = 220,seed1 = 1,seed2=2,seed3=3"""
import numpy as np
import matplotlib.pyplot as plt
import pickle

def RiskCurve(repeats,nus_trans):
    """Collects data for and plots the risk curve. Saves data."""
    AF_risks = [] # averaged time in AF
    AF_risks_std = []
    v_tAFs = [] # list of lists of each value of time in AF for each repeat
    #nus = [0,0.001,0.01,0.02,0.05,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.3,0.5,1]
    #nus_K = [0.02,0.04,0.06,0.08,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18
           #,0.19,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.3,1]
    #nus = [0,0.005,0.01,0.02,0.05,0.08,0.09,0.1,0.11,0.12,0.13,
           #0.15,0.17,0.18,0.19,0.2,0.21,0.3,0.5,1]
    nus = nus_trans 
    seed1 = 1     
    seed2 = 2
    seed3 = 3
    seed4 = 4
    for j in nus:
        v_tAF = []
        print(j)
        for i in range(repeats):
            print(i)
            A = AC.Atrium(v_tran=j, seed1=seed1, seed2=seed2, seed3=seed3,seed4 = seed4)
            A.CMP2D_time_AF()
            v_tAF.extend([A.tot_AF/A.tot_time])
            seed1 +=20
            seed2 +=20
            seed3 +=20
            seed4 +=20
        sample_avg = sum(v_tAF)/repeats
        sample_std = np.std(v_tAF)
        v_tAFs.extend([v_tAF])
        AF_risks.extend([sample_avg])
        AF_risks_std.extend([sample_std])
    data = [nus,AF_risks,v_tAFs]
    #pickle.dump(data,open( "risk_curve_0.5-1.p", "wb" ) )
    return data
#nus_trans = [0.02,0.04,0.06]    
#RiskCurve(10,nus_trans)