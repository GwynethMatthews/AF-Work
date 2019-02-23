import TRYING as CA
import numpy as np
import random as rnd
import cProfile
import pstats
import _pickle as Pickle# import ujson 
#import StringIO
pr = cProfile.Profile()
pr.enable()
L = 200 # system size
d = 0.05 # dysfunctionality prob
v = 0.2  # transverse connection prob
e = 0.05 # prob of dysfunctional cell not firing
timeperiod = 5000 # total number of time steps
pace_rate = 2 # time steps between sinus beats
Atrium, Phases = CA.CreateAtrium(L,v,d)
#print (Atrium)
#TempAtrium = np.ndarray([L,L],dtype = list)
#Atrium = SinusBeat(Atrium)
    
def CMP2D(Atrium,Phases, e, timeperiod, pace_rate):
    t = 0    
    while t < timeperiod:
        #TempAtrium[i] = Atrium[i].copy()
        #TempAtrium = Pickle.loads(Pickle.dumps(Atrium, -1))
        TempPhases = Phases.copy()
        #print(TempPhases)
        for i in range(L):
            for j in range(L):
                if TempPhases[i,j] != 4 and TempPhases[i,j] != 0:
                        Phases[i,j] = Phases[i,j] + 1
                if TempPhases[i,j] == 0:
                    for s in Atrium[i,j][1]:
                        if TempPhases[s] == 4:
                            if Atrium[s][0] == False:  # if dysfunctional
                                if e < rnd.uniform(0,1):
                                    Phases[s] = 0
                            if Atrium[s][0] == True:
                                Phases[s] = 0
                    Phases[i,j] = 1
        if t in np.arange(0,timeperiod,pace_rate): 
            for i in range(len(Atrium)):
                if TempPhases[i,0] == 4:
                    if Atrium[i,0][0] == False:  # if dysfunctional
                        if e < rnd.uniform(0,1):
                            Phases[i,0] = 0
                    if Atrium[i,0][0] == True:
                        Phases[i,0] = 0        
        #print(Atrium)
        #print(Phases)
        #print(TempPhases)
        t +=1   
            
    return Atrium
pr = cProfile.Profile()
pr.enable()
CMP2D(Atrium,Phases, e, timeperiod, pace_rate)
print('Atrium')
pr.disable()
#q = StringIO.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr).sort_stats(sortby)
ps.print_stats()