import numpy as np
import random as rnd
import cProfile
import numpy as np
import random as rnd
L=3 #size of atrium
d = 0.05 #probability of dysfunction
v = 0.2  # probability of transverse connection down

def CreateAtrium(L,v,d):
    """Creats the atrium. Atrium[i,j] gives a site (row,column). 
    Atrium[i,j[0] gives phase. 
    Atrium[i,j][1] dysfunctionality (False = dysfunctional), 
    Atrium[i,j[2] gives neightbouring sites"""    
    Atrium = np.ndarray(L*L,dtype = list)
    Phases = np.ndarray(L*L)
    Phases.fill(4)
    y = np.indices((1, L*L))[1]
    x  = y.reshape(L,L) # The index for that site within the long arrays
    for j in np.arange(L*L):
        z = rnd.uniform(0,1)
        if d > z: # dysfunctional
            Atrium[j] = [j,False,[]]
        if d <= z: # functional
            Atrium[j] = [j,True,[]]
        if j in np.arange(0,L*L,L): # first column
            Atrium[j][2].extend([(j+1)])
        elif j in (np.arange(0,L*L,L)+L-1): # last column
           # print(j)
            Atrium[j][2].extend([(j-1)])
        else: # body columns
            Atrium[j][2].extend([(j-1),(j+1)])
        w = rnd.uniform(0,1)
    for j in np.arange(L*L):
        if w <= v: # transverse connections
            if j in np.arange(L*L-L,L*L):
                Atrium[j][2].extend([j-(L*L-L)])
                Atrium[j-(L*L-L)][2].extend([j])
            else:
                Atrium[j][2].extend([j+L])
                Atrium[(j+L)][2].extend([j])
    return Atrium, Phases, x, y 
#Atrium, Phases, x = CreateAtrium(L,v,d)
#print(Atrium)
#print(Phases)
#print(x)