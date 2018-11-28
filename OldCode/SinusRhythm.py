import CreateAtrium1 as CA
import numpy as np
import random as rnd
import cProfile
import pstats
import _pickle as Pickle# import ujson 
import gc
#import StringIO
pr = cProfile.Profile()
pr.enable()
L = 200 # system size
d = 0.05 # dysfunctionality prob
v = 0.2  # transverse connection prob
e = 0.05 # prob of dysfunctional cell not firing
timeperiod = 100 # total number of time steps
pace_rate = 2 # time steps between sinus beats
Atrium = CA.CreateAtrium(L,v,d)
#print (Atrium)
TempAtrium = np.ndarray([L,L],dtype = list)
#Atrium = SinusBeat(Atrium)
    
def CMP2D(Atrium,TempAtrium, e, timeperiod, pace_rate):
    t = 0    
    while t < timeperiod:
        #TempAtrium[i] = Atrium[i].copy()
        TempAtrium = Pickle.loads(Pickle.dumps(Atrium, -1))
        #TempAtrium = json.loads(json.dumps(Atrium))
        #print(TempAtrium)
        for i in range(L):
            for j in range(L):
                if TempAtrium[i,j][0] != 4 and TempAtrium[i,j][0] != 0:
                        Atrium[i,j][0] = Atrium[i,j][0] + 1
                if TempAtrium[i,j][0] == 0:
                    for s in TempAtrium[i,j][2]:
                        if TempAtrium[s][0] == 4:
                            if TempAtrium[s][1] == False:  # if dysfunctional
                                if e < rnd.uniform(0,1):
                                    Atrium[s][0] = 0
                            if TempAtrium[s][1] == True:
                                Atrium[s][0] = 0
                    Atrium[i,j][0] = 1
        if t in np.arange(0,timeperiod,pace_rate): 
            for i in range(len(Atrium)):
                if TempAtrium[i,0][0] == 4:
                    if TempAtrium[i,0][1] == False:  # if dysfunctional
                        if e < rnd.uniform(0,1):
                            Atrium[i,0][0] = 0
                    if Atrium[i,0][1] == True:
                        Atrium[i,0][0] = 0        
        #print(Atrium)
        #print(TempAtrium)
       # print(TempAtrium1)
        t +=1   
            
    return Atrium
pr = cProfile.Profile()
pr.enable()
CMP2D(Atrium,TempAtrium, e, timeperiod, pace_rate)
print('Atrium')
pr.disable()
#q = StringIO.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr).sort_stats(sortby)
ps.print_stats()