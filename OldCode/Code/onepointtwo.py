""" 
version 1 (one)had all neighbours in one list,
version 2 (onepointtwo) put neighbours in seperate arrays depending on if they 
are up, down, left or right
"""
import random as rnd
import cProfile
import numpy as np
import random as rnd
#L=200 #size of atrium
#d = 0.05 #probability of dysfunction
#v = 0.2  # probability of transverse connection down
#seed = 1
def CreateAtrium(L,v,d,seed,refractory_period):
    """Creats the atrium. Atrium[i,j] gives a site (row,column). 
    Atrium[i,j[0] gives phase. 
    Atrium[i,j][1] dysfunctionality (False = dysfunctional), 
    Atrium[i,j[2] gives neightbouring sites"""    
    Neighbours_up = np.full((L*L),fill_value = None, dtype = float)
    Neighbours_down = np.full((L*L),fill_value = None, dtype = float)
    Neighbours_left = np.full((L*L),fill_value = None, dtype = float)
    Neighbours_right = np.full((L*L),fill_value = None, dtype = float)
    rnd.seed(seed)
    Phases = np.ndarray(L*L,dtype = float)
    Phases.fill(refractory_period)
    Functionality = np.ndarray([L*L], dtype = bool)
    index = np.indices((1, L*L))[1][0]
    Atrium  = index.reshape(L,L) # The index for that site within the long arrays
    w = np.random.rand(L*L)
    #print(w)
    for j in index:
        z = rnd.uniform(0,1)
        if d > z: # dysfunctional
            Functionality[j] = False
        if d <= z: # functional
            Functionality[j] = True
        if j in np.arange(0,L*L,L): # first column
            # right connection for the first column
            Neighbours_right[j] = j+1 
        elif j in (np.arange(0,L*L,L)+L-1): # last column
            # left connection for last column
            Neighbours_left[j] = j-1
        else: # body columns
            # currently all cells are connected longitudinally
            # left connection for body of cells
            Neighbours_left[j] = j-1
            # right connection for body of cells
            Neighbours_right[j] = j+1 
    for j in np.arange(L*L):
        # transverse connections (checks each one for a downward connection)
        if w[j] <= v: 
            #print(w)
            if j in np.arange(L*L-L,L*L):
                # the jth site has a downwards connection
                Neighbours_down[j] = j-(L*L-L)
                # the j-(L*L-L)th site (corresponding column in the first row
                # has an upwards connection
                Neighbours_up[j-(L*L-L)] = j
            else:
                # the jth site has a downwards connection
                Neighbours_down[j] = j+L
                # the j+Lth site has an upwards connection
                Neighbours_up[j+L] = j
    return Neighbours_up, Neighbours_down, Neighbours_right, Neighbours_left, Phases, Functionality, Atrium, index 
#Neighbours_up, Neighbours_down, Neighbours_right, Neighbours_left, Phases, Functionality, Atrium, index  = CreateAtrium(L,v,d,seed)
#print(Neighbours_up.reshape([L,L]))
#print(Neighbours_down.reshape([L,L]))
#print(Neighbours_right.reshape([L,L]))
#print(Neighbours_left.reshape([L,L]))
#print(Phases)
#print(Functionality)
#print(Atrium)
#print(index)