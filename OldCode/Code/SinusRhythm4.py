import CreateAtrium3 as CA
import numpy as np
import random as rnd
import cProfile
import pstats
import copy
#import StringIO

L = 200 # system size
d = 0.05 # dysfunctionality prob
v = 0.2  # transverse connection prob
e = 0.05 # prob of dysfunctional cell not firing
time = 200 # total number of time steps
pace_rate = 6 # time steps between sinus beats
Atrium, resting, first_column = CA.CreateAtrium(L,v,d)
tbe = []
preexcited = []
excited = []
r1 = []
r2 = []
r3 = []
nonrest = []
allsites = copy.deepcopy(resting)
#print(Atrium)
#def SinusRhythm(Atrium,resting):
def ChangePhase(Atrium,s,new_phase):
    s[1] = new_phase
    Atrium[s[0]][1] = new_phase
    return s

def Neighbours(s,neighbour_coordinate):
    neighbour_coordinate.extend(s[3])
    return neighbour_coordinate

def Conduct(excited,Atrium):
    neighbours = []
    neighbour_coordinate = []
    [neighbour_coordinate.extend(s[3]) for s in excited]
    neighbours = ([Atrium[n] for n in neighbour_coordinate])
    return neighbours

def ToBeExcited(first_row,neighbours, t):
    tbe = []
    tbe = neighbours
    if t in np.arange(0,time,pace_rate):
        tbe.extend(first_row)
    return tbe

def Excite(preexcited,resting,tbe):
   #print(resting)
    preexcited = [y for y in tbe if y[1] == True  and y[2] == True]
    preexcited.extend([x for x in tbe if x[1] ==True and x[2] == False and e < rnd.uniform(0,1)]) 
    tbe = []
    return preexcited

def Relax(Atrium,preexcited, excited,r1,r2,r3,resting):
    #resting = list(allsites)
#    print('resting')
#    print(resting)
#    print(r1)
#    print(r2)
#    print(excited)
    #resting1 = [resting.remove(s) for s in resting if s in r1]  
    #resting2 = [resting1.remove(p) for p in resting if p in r2] #[ChangePhase(Atrium, s,new_phase = 4) for s in r3])
    #resting3 = [resting2.remove(q) for q in resting if q in excited]
    #for r in preexcited:
    #[resting.remove(r) for r in preexcited if r in resting]
    #resting.extend(r3)
 #   print(resting)
    r3 = []
    r3 = r2 #[ChangePhase(Atrium,s,new_phase = 3) for s in r2]
    r2 = []
    r2 = r1 #[ChangePhase(Atrium,s,new_phase = 2) for s in r1]
    r1 = []
    r1 = excited #[ChangePhase(Atrium,s,new_phase = 1) for s in excited]
    excited = preexcited #[ChangePhase(Atrium,s,new_phase = 0) for s in preexcited]

    preexcited = []
    return preexcited,excited,r1,r2,r3,resting

def Update_atrium(Atrium,preexcited,r3):
    [ChangePhase(Atrium,s,False) for s in preexcited]
    [ChangePhase(Atrium,s,True) for s in r3]
    
def CMP(Atrium,first_row,tbe,preexcited,excited,r1,r2,r3,resting, time, pace_rate):
    t = 0
    while t < time:
        neighbours = Conduct(excited, Atrium)   
       # print('neighbours')
        #print(neighbours)
        tbe = ToBeExcited(first_row,neighbours,t)
#        print(t)
#        print('tbe')
#        print(tbe)
        t +=1
        preexcited = Excite(preexcited,resting,tbe)
#        print('preexcited')
        Update_atrium(Atrium,preexcited,r3)
        preexcited,excited,r1,r2,r3,resting,= Relax(Atrium,preexcited, excited,r1,r2,r3,resting)        
        
#        print('preexcited')
#        print(preexcited)
#        print('excited')
#        print(excited)
#        print('r1')
#        print(r1)
#        print('r2')
#        print(r2)
#        print('r3')
#        print(r3)
#        print('resting')
#        print(resting)
pr = cProfile.Profile()
pr.enable()        
CMP(Atrium, first_column,tbe,preexcited,excited,r1,r2,r3,resting, time, pace_rate)
print('Atrium')
pr.disable()
#q = StringIO.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr).sort_stats(sortby)
ps.print_stats()


    