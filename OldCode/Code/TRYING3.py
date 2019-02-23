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
timeperiod = 1000 # total number of time steps
pace_rate = 2 # time steps between sinus beats
Atrium, Phases = CA.CreateAtrium(L,v,d)
#print (Atrium)
#TempAtrium = np.ndarray([L,L],dtype = list)
#Atrium = SinusBeat(Atrium)
def Conduct(Atrium, Phases, TempPhases, e):
    ru = rnd.uniform
    for i in range(L):
        for j in range(L):
            if TempPhases[i,j] != 4 and TempPhases[i,j] != 0:
                    Phases[i,j] = Phases[i,j] + 1
            if TempPhases[i,j] == 0:
                for s in Atrium[i,j][1]:
                    if TempPhases[s] == 4:
                        if Atrium[s][0] == False:  # if dysfunctional
                            if e < ru(0,1):
                                Phases[s] = 0
                        if Atrium[s][0] == True:
                            Phases[s] = 0
                Phases[i,j] = 1    
    return Atrium, Phases

def SinusBeat(Atrium,Phases,TempPhases, e):    
    ru = rnd.uniform            
    for i in range(L):
        if TempPhases[i,0] == 4:
            if Atrium[i,0][0] == False:  # if dysfunctional
                if e < ru(0,1):
                    Phases[i,0] = 0
            if Atrium[i,0][0] == True:
                Phases[i,0] = 0       
    return Atrium, Phases        
         
def CMP2D(Atrium,Phases, e, timeperiod, pace_rate):
    t = 0    
    while t < timeperiod:
        TempPhases = Phases.copy()
        Conduct(Atrium, Phases, TempPhases, e)
        if t in np.arange(0,timeperiod,pace_rate):
            SinusBeat(Atrium,Phases,TempPhases, e)
        t +=1           
    return Atrium, Phases


pr = cProfile.Profile()
pr.enable()
CMP2D(Atrium,Phases, e, timeperiod, pace_rate)
print('Atrium')
pr.disable()
#q = StringIO.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr).sort_stats(sortby)
ps.print_stats()