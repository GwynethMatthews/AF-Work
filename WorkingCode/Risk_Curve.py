import Atrium as AC
import numpy as np
import pickle
"""Defaults: L=200,v=0.5,d=0.05,e=0.05,rp=50,tot_time = 10**6
                 ,pace_rate = 220,seed1 = 1,seed2=2,seed3=3"""
def RiskCurve(repeats,nus_trans,nus_para,seed1,seed2,seed3,seed4):
    """Collects data for and plots the risk curve. Saves data."""
    AF_risks = []
    v_tAFs = []
    nus_para = nus_para
    nus = nus_trans 
    for j in nus:
        v_tAF = []
        for i in range(repeats):
            A = AC.Atrium(hexagonal = False, v_para=nus_para,v_tran_1= j, pace_rate = 1980,seed1=seed1, seed2=seed2, seed3=seed3,seed4 = seed4)
            A.CMP2D_time_AF()
            v_tAF.extend([A.tot_AF/A.tot_time])
            seed1 +=20
            seed2 +=20
            seed3 +=20
            seed4 +=20
        sample_avg = sum(v_tAF)/repeats

        v_tAFs.extend([v_tAF])
        AF_risks.extend([sample_avg])
    data = [nus_para,nus,AF_risks,v_tAFs]
    return data

def RiskCurve2(repeats,dfs):
    """Collects data for and plots the risk curve. Saves data."""
    AF_risks = []
    v_tAFs = []
    seed1 = 1
    seed2 = 2
    seed3 = 3
    seed4 = 4
    for j in dfs:
        print(j)
        v_tAF = []
        for i in range(repeats):
            A = AC.Atrium(hexagonal = False, v_para=1,v_tran_1 = 0.3, d=i, seed1=seed1, seed2=seed2, seed3=seed3,seed4 = seed4)
            A.CMP2D_time_AF()
            v_tAF.extend([A.tot_AF/A.tot_time])
            seed1 +=20
            seed2 +=20
            seed3 +=20
            seed4 +=20
        sample_avg = sum(v_tAF)/repeats

        v_tAFs.extend([v_tAF])
        AF_risks.extend([sample_avg])
    data = [dfs,AF_risks,v_tAFs]
    pickle.dump(data,open( "dysfunc_cells7.p", "wb" ) )
repeats = 50
dfs = np.arange(0,1.05,0.05)
RiskCurve2(repeats,dfs[18:])