import CreateAtrium3 as CA
import numpy as np
import random as rnd
import cProfile
import pstats
import copy
#import StringIO

L = 2 # system size
d = 0.05 # dysfunctionality prob
v = 0.2  # transverse connection prob
e = 0.05 # prob of dysfunctional cell not firing
time = 1 # total number of time steps
pace_rate = 2 # time steps between sinus beats
Atrium, resting, first_column = CA.CreateAtrium(L,v,d)
tbe = []
preexcited = []
excited = []
r1 = []
r2 = []
r3 = []
nonrest = []
def ChangePhase(Atrium,s,new_phase):
    s[1] = new_phase
    Atrium[s[0]][1] = new_phase
    return s

def Neighbours(s,neighbour_coordinate):
    neighbour_coordinate = []
    neighbour_coordinate.extend(Atrium[n] for n in s[3])
    print('nc')
    print(neighbour_coordinate)
    return neighbour_coordinate

def Conduct(excited,Atrium):
    neighbours = []
    neighbour_coordinate = []
    print('e')
    print(excited)
    neighbours.extend([Neighbours(s,neighbour_coordinate) for s in excited])
    print('n')
    print(neighbours)
    return neighbours

def ToBeExcited(first_row,neighbours, t):
    tbe = []
    tbe = neighbours
    if t in np.arange(0,time,pace_rate):
        tbe.extend(first_row)
    return tbe

def Excite(preexcited,resting,tbe):
    print('tbe')
    print(tbe)
    preexcited = [y for y in tbe if y[1] == True and y[2] == True]
    preexcited.extend([x for x in tbe if x[1] == True and x[2] == False and e < rnd.uniform(0,1)]) 
    tbe = []
    return preexcited

def Relax(Atrium,preexcited, excited,r1,r2,r3,resting):
    r3 = []
    r3 = r2 
    r2 = []
    r2 = r1 
    r1 = []
    r1 = excited 
    excited = preexcited 
    preexcited = []
    return preexcited,excited,r1,r2,r3,resting

def Update_atrium(Atrium,preexcited,r3):
    [ChangePhase(Atrium,s,False) for s in preexcited]
    [ChangePhase(Atrium,s,True) for s in r3]
    
def CMP(Atrium,first_row,tbe,preexcited,excited,r1,r2,r3,resting, time, pace_rate):
    t = 0
    while t < time:
        neighbours = Conduct(excited, Atrium)   
        tbe = ToBeExcited(first_row,neighbours,t)
        t +=1
        preexcited = Excite(preexcited,resting,tbe)
        Update_atrium(Atrium,preexcited,r3)
        preexcited,excited,r1,r2,r3,resting,= Relax(Atrium,preexcited, excited,r1,r2,r3,resting)        

pr = cProfile.Profile()
pr.enable()        
CMP(Atrium, first_column,tbe,preexcited,excited,r1,r2,r3,resting, time, pace_rate)
print('Atrium')
pr.disable()
sortby = 'cumulative'
ps = pstats.Stats(pr).sort_stats(sortby)
ps.print_stats()


    