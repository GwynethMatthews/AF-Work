import CreateAtrium4 as CA
import numpy as np
import random as rnd
import cProfile
import pstats

import matplotlib.pyplot as plt
import matplotlib.animation as animation
#import StringIO
pr = cProfile.Profile()
pr.enable()
L = 200 # system size
d = 0.05 # dysfunctionality prob
v = 0.2  # transverse connection prob
e = 0.05 # prob of dysfunctional cell not firing
timeperiod = 10 # total number of time steps
pace_rate = 220 # time steps between sinus beats
Atrium, Phases, x, y = CA.CreateAtrium(L,v,d)
"""FLat Array will be faster, cython"""
def Conduct(Atrium, Phases, TempPhases, e,r,c):
    ru = rnd.uniform
    if TempPhases[r,c] != 4 and TempPhases[r,c] != 0:
            Phases[r,c] = Phases[r,c] + 1
    if TempPhases[r,c] == 0:
        for n in Atrium[r,c][1]:
            if TempPhases[n] == 4:
                if Atrium[n][0] == False:  # if dysfunctional
                    if e < ru(0,1):
                        Phases[n] = 0
                if Atrium[n][0] == True:
                    Phases[n] = 0
        Phases[r,c] = 1    
    return Phases

def SinusBeat(Atrium,Phases,TempPhases, e,i):    
    ru = rnd.uniform            
    if TempPhases[i,0] == 4:
        if Atrium[i,0][0] == False:  # if dysfunctional
            if e < ru(0,1):
                Phases[i,0] = 0
        if Atrium[i,0][0] == True:
            Phases[i,0] = 0       
    return Phases        
         
def CMP2D(L,Atrium,Phases,x,y, e, timeperiod, pace_rate):
    t = 0  
    print(Phases)
    while t < timeperiod:
        TempPhases = Phases.copy()
        result1 = list(map(lambda r,c: Conduct(Atrium, Phases, TempPhases, e,r,c), x, y))
        if t in np.arange(0,timeperiod,pace_rate):
            #print(t)
            result2 = list(map(lambda i: SinusBeat(Atrium,Phases,TempPhases, e,i), np.arange(0,L)))
        t +=1
        #print(Phases)
    return Phases

"""fig, ax = plt.subplots()
model = ax.matshow(Phases)
animation = animation.FuncAnimation(fig,CMP2D(L,Atrium,Phases,x,y, e, timeperiod, pace_rate),frames = timeperiod, repeat = False,save_count = 10)
plt.show()"""
pr = cProfile.Profile()
pr.enable()
CMP2D(L,Atrium,Phases,x,y, e, timeperiod, pace_rate)
print('Atrium')
pr.disable()
#q = StringIO.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr).sort_stats(sortby)
ps.print_stats()