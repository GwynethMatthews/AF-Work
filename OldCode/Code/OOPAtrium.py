import numpy as np
import random as rnd
import cProfile
import pstats

class Atrium:
    def __init__(self,L,v,d,e,time,pacemaker_rate):
        self.size = L*L
        self.shape = (L,L)
        #self.first_row = 
        self.excited = np.zeros([L,L])
        self.r1 = np.zeros([L,L])
        self.r2 = np.zeros([L,L])
        self.r3 = np.zeros([L,L])
        self.resting = np.ones([L,L])
        self.dysfuctional = np.empty([L,L])
        self.vertical_connections = np.empty([L,L])
        self.dys_prob = d
        self.dys_percent = e
        self.vert_prob = v
        self.t = 0
        self.pacemaker_rate = np.arange(0,time,pacemaker_rate)
        self.tot_time = time

    
        for i in range(L):
             for j in range(L):
                 if v < rnd.uniform(0,1):
                     self.vertical_connections[i,j] = False # 0 = no downwards connection
                 else:
                     self.vertical_connections[i,j] = True # 1 = downwards connection
             for j in range(L):
                 if d < rnd.uniform(0,1):
                     self.dysfuctional[i,j] = True # 1 = working
                 else:
                     self.dysfuctional[i,j] = False # 0 = dysfunctionL
       
        def Sinus_Rhythm(self, pacemaker_time, t):
            if self.t in self.pacemaker_time:
                self.excited.extend()
                
                
        
a1 = Atrium(2,0.2,0.05,0.05, 7, 70)
#print(a1.shape)
#print(a1.excited)
#print(a1.resting)
#print(a1.r1)
#print(a1.dysfuctional)
print(a1.cell_grid)
index = np.arange(a1.size, step=a1.shape[1])
print(index)
index = index[a1.cell_grid[index]]
print(index)