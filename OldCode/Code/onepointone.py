"""onepoinone uses the np.where function, 
version one (onepointone) used contantenate: 
excited = np.where(TempPhases == 0)[0]
# finds the indices of neighbours of the excited cells
neighbours = np.concatenate(Neighbours[excited])
# removes the neighbours that aren't resting
neighbours =  neighbours[np.where(TempPhases[neighbours] == 4)[0]]

using np.concatenate is what takes up most of the time

version two (onepointthree) try seperating the Neighbour lists into right, left, up and down

"""

import one as CA
import numpy as np
import random as rnd
import cProfile
import pstats
#import StringIO

L = 200 # system size
d = 0.05 # dysfunctionality prob
v = 0.2  # transverse connection prob
e = 0.05 # prob of dysfunctional cell not firing
tot_time = 10**4 # total number of time steps
pace_rate = 220 # time steps between sinus beats
seed = 1
#Neighbours, Phases, Functionality, Atrium, index = CA.CreateAtrium(L,v,d,seed)
#first_col = np.arange(0, L*L, L, dtype = int)

def SinusRhythm(TempPhases, Phases,e,Functionality):
    # finds index of resting cells
    resting = np.where(TempPhases == 4)[0]  
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

def Relaxing(TempPhases, Phases):
    # finds the index of cells which are not resting i.e in phases 0-3
    non_resting = np.where(TempPhases != 4)[0]
    # updates the phases of the non-resting states by 1
    Phases[non_resting] += 1
    

def Excite(TempPhases, Phases, Neighbours, Functionality, e):
    # finds the index of cells that are excited
    excited = np.where(TempPhases == 0)[0]
    # finds the indices of neighbours of the excited cells
    neighbours = np.concatenate(Neighbours[excited])
    # removes the neighbours that aren't resting
    neighbours =  neighbours[np.where(TempPhases[neighbours] == 4)[0]]
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

def CMP2D(L,v,d,seed, tot_time, pace_rate):   
    Neighbours, Phases, Functionality, Atrium, index = CA.CreateAtrium(L,v,d,seed)
    t = 0
    #print(Phases)
    while t < tot_time:
        #print(t)
        TempPhases1 = Phases.copy()
        #print(TempPhases1)
        if t in np.arange(0,tot_time,pace_rate):
            SinusRhythm(TempPhases1,Phases,e,Functionality)
        if np.all(TempPhases1) == False:
            Excite(TempPhases1, Phases, Neighbours, Functionality, e)
        Relaxing(TempPhases1, Phases)

        #print(Phases)
        t += 1
    return Phases


pr = cProfile.Profile()
pr.enable()
CMP2D(L,v,d,seed, tot_time, pace_rate)
print('Atrium')
pr.disable()
#q = StringIO.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr).sort_stats(sortby)
ps.print_stats()

