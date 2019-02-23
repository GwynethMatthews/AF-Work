import Flattenedcode as CA
import numpy as np
import random as rnd
import cProfile
import pstats
#import StringIO

L = 5 # system size
d = 0.05 # dysfunctionality prob
v = 0.2  # transverse connection prob
e = 0.05 # prob of dysfunctional cell not firing
time = 5 # total number of time steps
pace_rate = 5 # time steps between sinus beats
#Atrium, Phases, x, y = CA.CreateAtrium(L,v,d)
#first_col = np.arange(0,L*L,L)
tbe = []
preexcited = []
excited = []
r1 = []
r2 = []
r3 = []
neighbours = []
#print(Atrium)
#def SinusRhythm(Atrium,resting):
def ChangePhase(Phases,s,new_phase):
    Phases[s] = new_phase

def Neighbours(s,neighbour_coordinate):
    neighbour_coordinate.extend(s[2])
    return neighbour_coordinate

def Conduct(excited,Atrium):
    neighbours = []
    neighbour_coordinate = []
    neighbour_coordinate.extend(Atrium[s][2] for s in excited)
    neighbours.extend(neighbour_coordinate)
    return neighbours

def ToBeExcited(first_col,neighbours, t,time,pace_rate):
    tbe = neighbours
    if t in np.arange(0,time,pace_rate):
        tbe.extend(first_col)
    return tbe

def Excite(Atrium, Phases,preexcited,tbe,e):
    
    print(tbe)
    preexcited = [s for s in tbe if Phases[s]==4 and Atrium[s][1] == True]
    preexcited.extend([s for s in tbe if Phases[s]==4 and Atrium[s][1] == False and e < rnd.uniform(0,1)]) 
    tbe = []
    return preexcited

def Relax(Atrium,preexcited, excited,r1,r2,r3):
    r3 = r2
    r2 = r1
    r1 = excited
    excited = preexcited
    preexcited = []
    return preexcited,excited,r1,r2,r3
def Update_atrium(Phases,preexcited,excited,r1,r2,r3):
    [ChangePhase(Phases,s,0) for s in preexcited]
    [ChangePhase(Phases,s,1) for s in excited]
    [ChangePhase(Phases,s,2) for s in r1]
    [ChangePhase(Phases,s,3) for s in r2]
    [ChangePhase(Phases,s,4) for s in r3]
    
def CMP(L,v,d,e,time,pace_rate):
    t = 0
    tbe = []
    preexcited = []
    excited = []
    r1 = []
    r2 = []
    r3 = []
    neighbours = []
    Atrium, Phases, x, y = CA.CreateAtrium(L,v,d)
    first_col = np.arange(0,L*L,L)
    print(Phases)
    while t < time:
        neighbours = Conduct(excited, Atrium)     
        print(neighbours)
        tbe = ToBeExcited(first_col,neighbours,t,time,pace_rate)
        print(tbe)
        t +=1
        preexcited = Excite(Atrium, Phases, preexcited,tbe,e)
        print(preexcited)
        Update_atrium(Phases,preexcited,excited,r1,r2,r3)
        preexcited,excited,r1,r2,r3 = Relax(Atrium,preexcited, excited,r1,r2,r3)
        print(preexcited)
        print(excited)
        print(r1)
        print(r2)
        print(r3)        
        print(Phases)
        
pr = cProfile.Profile()
pr.enable()        
CMP(L,v,d,e,time,pace_rate)
print('Atrium')
pr.disable()
#q = StringIO.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr).sort_stats(sortby)
#ps.print_stats()


    