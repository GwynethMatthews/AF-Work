""" 
version 1 (one)had all neighbours in one list,
version 2 (onepointtwo) put neighbours in seperate arrays depending on if they
are up, down,
left or right
"""
import random as rnd
import cProfile
import numpy as np
import random as rnd
L=3 #size of atrium
d = 0.05 #probability of dysfunction
v = 0.2  # probability of transverse connection down
seed = 1
def CreateAtrium(L,v,d,seed):
    """Creats the atrium. Atrium[i,j] gives a site (row,column). 
    Atrium[i,j[0] gives phase. 
    Atrium[i,j][1] dysfunctionality (False = dysfunctional), 
    Atrium[i,j[2] gives neightbouring sites"""    
    Neighbours = np.ndarray(L*L, dtype = list)
    rnd.seed(seed)
    Phases = np.ndarray(L*L,dtype = float)
    Phases.fill(4)
    Functionality = np.ndarray([L*L], dtype = bool)
    index = np.indices((1, L*L))[1][0]
    Atrium  = index.reshape(L,L) # The index for that site within the long arrays
    for j in index:
        z = rnd.uniform(0,1)
        if d > z: # dysfunctional
            Functionality[j] = False
        if d <= z: # functional
            Functionality[j] = True
        if j in np.arange(0,L*L,L): # first column
            Neighbours[j] = list()
            Neighbours[j].extend([j+1])
        elif j in (np.arange(0,L*L,L)+L-1): # last column
            Neighbours[j] = list()
            Neighbours[j].extend([j-1])
        else: # body columns
            Neighbours[j] = list()
            Neighbours[j].extend([j-1,j+1])
        w = rnd.uniform(0,1)
    for j in np.arange(L*L):
        if w <= v: # transverse connections
            if j in np.arange(L*L-L,L*L):
                Neighbours[j].extend([j-(L*L-L)])
                Neighbours[j-(L*L-L)].extend([j])
            else:
                Neighbours[j].extend([j+L])
                Neighbours[(j+L)].extend([j])
    return Neighbours, Phases, Functionality, Atrium, index 
#Neighbours, Phases,Functionality, Atrium, index = CreateAtrium(L,v,d,seed)
#print(Neighbours)
#print(Phases)
#print(Functionality)
#print(Atrium)
#print(index)