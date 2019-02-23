import CreateAtrium3 as CA
import numpy as np
import random as rnd
import cProfile
import pstats
#import StringIO

L = 200 # system size
d = 0.05 # dysfunctionality prob
v = 0.2  # transverse connection prob
e = 0.05 # prob of dysfunctional cell not firing
time = 100 # total number of time steps
pace_rate = 5 # time steps between sinus beats
Atrium, resting, first_column = CA.CreateAtrium(L,v,d)
#print(first_column)
tbe = []
preexcited = []
excited = []
r1 = []
r2 = []
r3 = []
neighbours = []
#print(Atrium)
#def SinusRhythm(Atrium,resting):
def ChangePhase(Atrium,s,new_phase):
    Atrium[s][1] = new_phase
    return s


def Conduct(excited,Atrium):
    neighbours = []
    #neighbour_coordinate = []
    [neighbours.extend(Atrium[s][3]) for s in excited]
    #neighbours = ([Atrium[n] for n in neighbour_coordinate])
    #print(neighbours)
    return neighbours

def ToBeExcited(first_column,neighbours, t):
    tbe = neighbours
    if t in np.arange(0,time,pace_rate):
        tbe.extend(first_column)
    #print(tbe)
    return tbe

def Excite(preexcited,resting,tbe):
    preexcited = [s for s in tbe if Atrium[s][1]==4 and Atrium[s][2] == True]
    preexcited.extend([s for s in tbe if Atrium[s][1]==4 and Atrium[s][2] == False and e < rnd.uniform(0,1)]) 
    tbe = []
    return preexcited

def Relax(Atrium,preexcited, excited,r1,r2,r3,resting):
    resting.extend([ChangePhase(Atrium, s,new_phase = 4) for s in r3])
    r3 =[ChangePhase(Atrium,s,new_phase = 3) for s in r2]
    r2 =[ChangePhase(Atrium,s,new_phase = 2) for s in r1]
    r1 = [ChangePhase(Atrium,s,new_phase = 1) for s in excited]
    excited = [ChangePhase(Atrium,s,new_phase = 0) for s in preexcited]
    preexcited = []
    return preexcited,excited,r1,r2,r3,resting

def CMP(Atrium,first_row,tbe,neighbours,preexcited,excited,r1,r2,r3,resting, time, pace_rate):
    t = 0
    while t < time:
        neighbours = Conduct(excited, Atrium)       
        tbe = ToBeExcited(first_row,neighbours,t)
        t +=1
        preexcited = Excite(preexcited,resting,tbe)
        preexcited,excited,r1,r2,r3,resting = Relax(Atrium,preexcited, excited,r1,r2,r3,resting)        
pr = cProfile.Profile()
pr.enable()        
CMP(Atrium, first_column,tbe,neighbours,preexcited,excited,r1,r2,r3,resting, time, pace_rate)
print('Atrium')
pr.disable()
#q = StringIO.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr).sort_stats(sortby)
ps.print_stats()
