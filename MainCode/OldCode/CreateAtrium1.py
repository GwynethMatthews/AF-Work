import numpy as np
import random as rnd
L=5 #size of atrium
d = 0.05 #probability of dysfunction
v = 0.2  # probability of transverse connection down

def CreateAtrium(L,v,d):
    """Creats the atrium. Atrium[i,j] gives a site (row,column). 
    Atrium[i,j[0] gives phase. 
    Atrium[i,j][1] dysfunctionality (False = dysfunctional), 
    Atrium[i,j[2] gives neightbouring sites"""    
    Atrium = np.ndarray((L,L),dtype = list)
    for i in range(L):
        for j in range(L):
            y = rnd.uniform(0,1)
            if d > y:
                Atrium[i,j] = [4, False, []]
            if d <= y:
                Atrium[i,j] = [4, True,[]]
    for i in range(L):
        for j in range(L):
            if j == 0:
                Atrium[i,j][2].extend([(i,j+1)])
            if j == L-1:
                Atrium[i,j][2].extend([(i,j-1)])
            if j>0 and j<L-1:
                Atrium[i,j][2].extend([(i,j-1)])
                Atrium[i,j][2].extend([(i,j+1)])
            x = rnd.uniform(0,1)
            if x <= v:
                if i == L-1:
                    Atrium[i,j][2].extend([(0,j)])
                    Atrium[0,j][2].extend([(i,j)])
                if i != L-1:
                    Atrium[i,j][2].extend([(i+1,j)])
                    Atrium[i+1,j][2].extend([(i,j)])
    return Atrium
Atrium = CreateAtrium(L,v,d)