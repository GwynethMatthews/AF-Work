import Risk_Curve as RC
"""
Use coarse graining to find the rate of conduction and stuff.
Can't use 9*pace_rate.
Last years masters papers have the phase singularity.
Phase space for Hex and Square.
Target cells on the right to find percolation, NO DYSFUNCTIONAL CELLS,
If excited cells == 0 stop recording. 

Watch the animations to see how AF forms. Look at neighbourhood.
Try assiging a random connection value to each connections in the Hex grid. 
Look at statistics of connections that "how many of my neighbours are excited", measure dynamically, multiple cells
Can there be paroxysmal AF if it's based on connections?

"""

import pickle 
#import numpy as np 
def Risk_Plane(repeats, nus_para):
    data_full = []
    seed1 = 1
    seed2 = 2
    seed3 = 3
    seed4 = 4
    for k in nus_para:
        seed1 += 200
        seed2 += 200
        seed3 += 200
        seed4 += 200
        print(k)
        if k <= 0.2:
            nus_trans = [0.02,0.5,0.7,0.72,0.74,0.76,0.78,0.8,0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1]
        if k > 0.2 and k <= 0.3:
            nus_trans = [0.02,0.5,0.6,0.62,0.64,0.66,0.68,0.7,0.72,0.74,0.76,0.78,0.8,0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1]
        if k > 0.3 and k <= 0.4:
            nus_trans = [0.02,0.5,0.52,0.54,0.56,0.58,0.6,0.62,0.64,0.66,0.68,0.7,0.72,0.74,0.76,0.78,0.8,0.82,0.84,0.86,0.88,0.9,1]
        if k > 0.4 and k <= 0.5:
            nus_trans = [0.02,0.36,0.38,0.4,0.42,0.44,0.46,0.48,0.5,0.52,0.54,0.56,0.58,0.6,0.62,0.64,0.66,0.68,0.7,0.72,0.74,0.76,0.78,0.8,0.82,0.84,0.86,1]
        if k > 0.5 and k <= 0.6:
            nus_trans = [0.02,0.3,0.32,0.34,0.36,0.38,0.4,0.42,0.44,0.46,0.48,0.5,0.52,0.54,0.56,0.58,0.6,0.62,0.64,0.66,0.68,0.7,1]
        if k > 0.6 and k <= 0.7:
            nus_trans = [0.02,0.2,0.22,0.24,0.26,0.28,0.3,0.32,0.34,0.36,0.38,0.4,0.42,0.44,0.46,0.48,0.5,0.52,0.54,0.56,0.58,0.6,1]
        if k > 0.7 and k <= 0.8:
            nus_trans = [0.02,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.32,0.34,0.36,0.38,0.4,0.42,0.44,0.46,0.48,0.5,1]
        if k > 0.8:
            nus_trans = [0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.32,0.34,0.36,0.38,0.4,1]
        #nus_trans = np.arange(0,1.02,0.02)
        data = RC.RiskCurve(repeats,nus_trans,k,seed1,seed2,seed3,seed4)
        data_full.extend([data])
    pickle.dump(data_full,open( "risk_plane15.p", "wb" ) )

repeats = 5
nus_para = [0.9,0.95,1]
Risk_Plane(repeats,nus_para)