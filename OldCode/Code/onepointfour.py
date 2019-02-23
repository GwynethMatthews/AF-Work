
"""onepoinone uses the np.where function, 
version one used contantenate: 
excited = np.where(TempPhases == 0)[0]
# finds the indices of neighbours of the excited cells
neighbours = np.concatenate(Neighbours[excited])
# removes the neighbours that aren't resting
neighbours =  neighbours[np.where(TempPhases[neighbours] == 4)[0]]

using np.concatenate is what takes up most of the time

version two (onepointthree) seperates the Neighbour lists into right,
left, up and down, takes 300s for 10**5

"""

import onepointtwo as CA
import numpy as np
#import random as rnd
import cProfile
import pstats
#import StringIO

L = 200 # system size
d = 0.05 # dysfunctionality prob
v = 0.2  # transverse connection prob
e = 0.05 # prob of dysfunctional cell not firing
tot_time = 10**6# total number of time steps
refractory_period = 51
pace_rate = 220 # time steps between sinus beats
seed = 1
#Neighbours, Phases, Functionality, Atrium, index = CA.CreateAtrium(L,v,d,seed)
#first_col = np.arange(0, L*L, L, dtype = int)

def SinusRhythm(TempPhases, Phases,e,Functionality):
    # finds index of resting cells
    resting = np.where(TempPhases == refractory_period)[0]  
    # resting states that are in the first column
    condition1 = np.remainder(resting, L) == 0  
    # index of resting states that are in the first column
    tbe1 = np.extract(condition1, resting)
    # resting states that are in the first column and are dysfunctional
    # (does not give the index of the cell site gives the index of the 
    # shortened array of Functionality[tbe1])
    working1 = np.where(Functionality[tbe1] == False)[0] 
    # creates an array of random number between 0 and 1 based on a uniform 
    # distribution
    e_comp_val1 = np.random.rand(len(working1))
    # index of e_comp_val that are above e and correspond to excitable states
    gets_excited1 = np.where(e_comp_val1 > e)[0]
    # the index of tbe1 that corresponds to the excitable states sites index
    # working1 for dysfunctional but working cells and working2 for functional
    # cells
    working1 = working1[gets_excited1]
    working2 = np.where(Functionality[tbe1] == True)[0]
    # changes the phase of the excited states to 0
    Phases[tbe1[working1]] = 0
    Phases[tbe1[working2]] = 0

def Relaxing(TempPhases, Phases, refractory_period):
    # finds the index of cells which are not resting i.e in phases 0-3
    non_resting = np.where(TempPhases != refractory_period)[0]
    # updates the phases of the non-resting states by 1
    Phases[non_resting] += 1
    

def Excite(TempPhases, Phases, Neighbours_up, Neighbours_down, Neighbours_right, Neighbours_left, Functionality, e,refractory_period):
    # finds the index of cells that are excited
    excited = np.where(TempPhases == 0)[0]
    # finds the indices of neighbours of the excited cells
    #print('Here')
    #print(excited)
    neighbours = np.empty(0)
    have_neighbours_up = np.where(np.isnan(Neighbours_up[excited]) == False)
    neighbours_up = np.array((Neighbours_up[excited])[have_neighbours_up], dtype = int)
    #print(neighbours_up)
    have_neighbours_down = np.where(np.isnan(Neighbours_down[excited]) == False)
    neighbours_down = np.array((Neighbours_down[excited])[have_neighbours_down], dtype = int)
    #print(neighbours_down)
    have_neighbours_left = np.where(np.isnan(Neighbours_left[excited]) == False)
    neighbours_left = np.array((Neighbours_left[excited])[have_neighbours_left], dtype = int)
    #print(neighbours_left)
    have_neighbours_right = np.where(np.isnan(Neighbours_right[excited]) == False)[0]
    neighbours_right = np.array((Neighbours_right[excited])[have_neighbours_right],dtype = int)
    #print(neighbours_right)
#    condition2 = (neighbours_up in neighbours) == False
#    condition3 = (neighbours_down in neighbours) == False
#    condition4 = (neighbours_left in neighbours) == False
#    condition5 = (neighbours_right in neighbours) == False
#    neighbours.extend(np.extract(condition4,neighbours_left))
#    neighbours.extend(np.extract(condition5,neighbours_right))
#    neighbours.extend(np.extract(condition2,neighbours_up))
#    neighbours.extend(np.extract(condition3,neighbours_down))

    #neighbours = list(set(itertools.chain(neighbours_up, neighbours_down, neighbours_right, neighbours_left)))

    
    #neighbours = np.unique(np.concatenate((neighbours_up,neighbours_down,neighbours_left,neighbours_right)),axis=None)
    #c = Counter(neighbours)
    # removes the neighbours that aren't resting
    #print(neighbours)
    #neighbours = np.array(neighbours, dtype = int)
    neighbours = np.append(neighbours,neighbours_right[np.where(TempPhases[neighbours_right] == refractory_period)])
    neighbours = np.append(neighbours,neighbours_left[np.where(TempPhases[neighbours_left] == refractory_period)])
    neighbours = np.append(neighbours,neighbours_up[np.where(TempPhases[neighbours_up] == refractory_period)])
    neighbours = np.append(neighbours,neighbours_down[np.where(TempPhases[neighbours_down] == refractory_period)])
    neighbours = np.array(neighbours, dtype = int)
    # finds which neighbours are dysfunctional
    working3 = np.where(Functionality[neighbours] == False)[0]
    # creates an array of random numbers between 0 and 1 based on a uniform 
    # distribution
    e_comp_val2 = np.random.rand(len(working3))
    # index of e_comp_val that are above e and correspond to excitable states
    gets_excited2 = np.where(e_comp_val2 > e)[0]
    # the index of tbe1 that corresponds to the excitable states sites index
    # working1 for dysfunctional but working cells and working2 for functional
    # cells
    working3 = working3[gets_excited2]
    working4 = np.where(Functionality[neighbours] == True)[0]
    # changes the phase of the excited states to 0
    Phases[neighbours[working3]] = 0
    Phases[neighbours[working4]] = 0

def CMP2D(L,v,d,seed, tot_time, pace_rate, refractory_period):   
    Neighbours_up, Neighbours_down, Neighbours_right, Neighbours_left, Phases, Functionality, Atrium, index  = CA.CreateAtrium(L,v,d,seed,refractory_period)
    t = 0
    #print(Phases)
    while t < tot_time:
        #print(t)
        TempPhases1 = Phases.copy()
        #print(TempPhases1)
        if t in np.arange(0,tot_time,pace_rate):
            SinusRhythm(TempPhases1,Phases,e,Functionality)
        if np.all(TempPhases1) == False:
            Excite(TempPhases1, Phases, Neighbours_up, Neighbours_down, Neighbours_right, Neighbours_left, Functionality, e,refractory_period)
        Relaxing(TempPhases1, Phases, refractory_period)

        #print(Phases.reshape([L,L]))
        t += 1
    return Phases


pr = cProfile.Profile()
pr.enable()
CMP2D(L,v,d,seed, tot_time, pace_rate, refractory_period)
print('Atrium')
pr.disable()
#q = StringIO.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr).sort_stats(sortby)
ps.print_stats()

