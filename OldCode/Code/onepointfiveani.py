
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import cProfile
import pstats
L = 200 # system size
d = 0.05 # dysfunctionality prob
v = 0.12  # transverse connection prob
e = 0.05 # prob of dysfunctional cell not firing
tot_time = int(10**6) # total number of time steps
ref_period = 51
pace_rate = np.arange(0,tot_time,220) # time steps between sinus beats
seed1 = 1
seed2 = 2
#L=200 #size of atrium
#d = 0.05 #probability of dysfunction
#v = 0.2  # probability of transverse connection down
#seed = 1
def CreateAtrium(L,v,d,seed,ref_period):
    """Creats the atrium. Each attribute stored in a seperate list"""    
    Neigh_u = np.full((L*L),fill_value = None, dtype = float)
    Neigh_d = np.full((L*L),fill_value = None, dtype = float)
    Neigh_l = np.full((L*L),fill_value = None, dtype = float)
    Neigh_r = np.full((L*L),fill_value = None, dtype = float)
    np.random.seed(seed)
    Phases = np.ndarray(L*L,dtype = float)
    Phases.fill(ref_period)
    Funct = np.ndarray([L*L], dtype = bool)
    index = np.indices((1, L*L))[1][0]
    #Atrium  = index.reshape(L,L) # The index for that site within the long arrays
    w = np.random.rand(L*L)
    #print(w)
    for j in index:
        z = np.random.uniform(0,1)
        if d > z: # dysfunctional
            Funct[j] = False
        if d <= z: # functional
            Funct[j] = True
        if j in np.arange(0,L*L,L): # first column
            # right connection for the first column
            Neigh_r[j] = j+1 
        elif j in (np.arange(0,L*L,L)+L-1): # last column
            # left connection for last column
            Neigh_l[j] = j-1
        else: # body columns
            # currently all cells are connected longitudinally
            # left connection for body of cells
            Neigh_l[j] = j-1
            # right connection for body of cells
            Neigh_r[j] = j+1 
    for j in np.arange(L*L):
        # transverse connections (checks each one for a downward connection)
        if w[j] <= v: 
            #print(w)
            if j in np.arange(L*L-L,L*L):
                # the jth site has a downwards connection
                Neigh_d[j] = j-(L*L-L)
                # the j-(L*L-L)th site (corresponding column in the first row
                # has an upwards connection
                Neigh_u[j-(L*L-L)] = j
            else:
                # the jth site has a downwards connection
                Neigh_d[j] = j+L
                # the j+Lth site has an upwards connection
                Neigh_u[j+L] = j
    return Neigh_u, Neigh_d, Neigh_r, Neigh_l, Phases, Funct

def SinusRhythm(TPhases, Phases,e,Funct):
    """Pacemaker activity"""
    # finds index of resting cells
    resting = np.where(TPhases == ref_period)[0]  
    # resting states that are in the first column
    condition1 = np.remainder(resting, L) == 0  
    # index of resting states that are in the first column
    tbe1 = np.extract(condition1, resting)
    # resting states that are in the first column and are dysfunctional
    # (does not give the index of the cell site gives the index of the 
    # shortened array of Funct[tbe1])
    working1 = np.where(Funct[tbe1] == False)[0] 
    # creates an array of random number between 0 and 1 based on a uniform 
    # distribution
    e_comp_val1 = np.random.rand(len(working1))
    # index of e_comp_val that are above e and correspond to excitable states
    gets_excited1 = np.where(e_comp_val1 > e)[0]
    # the index of tbe1 that corresponds to the excitable states sites index
    # working1 for dysfunctional but working cells and working2 for functional
    # cells
    working1 = working1[gets_excited1]
    working2 = np.where(Funct[tbe1] == True)[0]
    # changes the phase of the excited states to 0
    Phases[tbe1[working1]] = 0
    Phases[tbe1[working2]] = 0

def Relaxing(TPhases, Phases, ref_period):
    """All cells move to the next phase"""
    # finds the index of cells which are not resting i.e in phases 0 to 
    # refractor_period - 1
    non_resting = np.where(TPhases != ref_period)[0]
    # updates the phases of the non-resting states by 1
    Phases[non_resting] += 1
    

def Excite(TPhases, Phases,Neigh_u,Neigh_d,Neigh_r,Neigh_l, Funct,e,ref_period):
    """Finds neighbours and excites them"""
    # finds the index of cells that are excited
    excited = np.where(TPhases == 0)[0]
    # creates the neighbour array
    neighbours = np.empty(0)
    # finds which excited cells that have up neighbours and then makes a list 
    # of the neighborus
    have_neigh_up = np.where(np.isnan(Neigh_u[excited]) == False)
    neighbours_up = np.array((Neigh_u[excited])[have_neigh_up], dtype = int)
    # finds which excited cells that have down neighbours and then makes a list 
    # of the neighborus
    have_neigh_down = np.where(np.isnan(Neigh_d[excited]) == False)
    neighbours_down = np.array((Neigh_d[excited])[have_neigh_down], dtype = int)
    # finds which excited cells that have left neighbours and then makes a list 
    # of the neighborus
    have_neigh_left = np.where(np.isnan(Neigh_l[excited]) == False)
    neighbours_left = np.array((Neigh_l[excited])[have_neigh_left], dtype = int)
    # finds which excited cells that have right neighbours and then makes a list 
    # of the neighborus
    have_neigh_right = np.where(np.isnan(Neigh_r[excited]) == False)[0]
    neighbours_right = np.array((Neigh_r[excited])[have_neigh_right],dtype = int)
    # makes a long list of all neighbours
    neighbours = np.append(neighbours,neighbours_right[np.where(TPhases[neighbours_right] == ref_period)])
    neighbours = np.append(neighbours,neighbours_left[np.where(TPhases[neighbours_left] == ref_period)])
    neighbours = np.append(neighbours,neighbours_up[np.where(TPhases[neighbours_up] == ref_period)])
    neighbours = np.append(neighbours,neighbours_down[np.where(TPhases[neighbours_down] == ref_period)])
    #neighbours = np.array(neighbours, dtype = int)
    neighbours = np.array(np.unique(neighbours),dtype = int)
    # finds which neighbours are dysfunctional
    working3 = np.where(Funct[neighbours] == False)[0]
    # creates an array of random numbers between 0 and 1 based on a uniform 
    # distribution
    e_comp_val2 = np.random.rand(len(working3))
    # index of e_comp_val that are above e and correspond to excitable states
    gets_excited2 = np.where(e_comp_val2 > e)[0]
    # the index of tbe1 that corresponds to the excitable states sites index
    # working3 for dysfunctional but working cells and working4 for functional
    # cells
    working3 = working3[gets_excited2]
    working4 = np.where(Funct[neighbours] == True)[0]
    # changes the phase of the excitable states to 0
    Phases[neighbours[working3]] = 0
    Phases[neighbours[working4]] = 0

def CMP2D(L,v,d,seed2, tot_time, pace_rate, ref_period):   
    # creates the array
    Neigh_u, Neigh_d, Neigh_r, Neigh_l, Phases, Funct  = CreateAtrium(L,v,d,seed2,ref_period)
    t = 0
    while t < tot_time:
        #print(t)
        # runs the time step change
        CMP2d_1step(seed2, pace_rate, ref_period,t, Neigh_u, Neigh_d, Neigh_r, Neigh_l, Phases, Funct)
        t += 1
        #print(Phases.reshape(L,L))

    
def CMP2d_1step(seed, pace_rate, ref_period,t, Neigh_u, Neigh_d, Neigh_r, Neigh_l, Phases, Funct):
    TPhases1 = Phases.copy()
    #np.random.seed(seed2)
    #print(Phases.reshape(L,L))
    if t in pace_rate:
        SinusRhythm(TPhases1,Phases,e,Funct)
        #print(Phases.reshape(L,L))
    #if np.all(TPhases1) == False:
    Excite(TPhases1, Phases, Neigh_u, Neigh_d, Neigh_r, Neigh_l, Funct, e,ref_period)
    #print(Phases.reshape(L,L))
    Relaxing(TPhases1, Phases, ref_period)
    #print(Phases.reshape(L,L))
    return Phases

def update(frame_number, mat,seed, pace_rate, ref_period, Neigh_u, Neigh_d, Neigh_r, Neigh_l, Phases, Funct):
    TPhases2 = Phases.copy()
    #np.random.seed(seed)
    if frame_number in pace_rate:
        SinusRhythm(TPhases2,Phases,e,Funct)
    #if np.all(TPhases2) == False:
    Excite(TPhases2, Phases, Neigh_u, Neigh_d, Neigh_r, Neigh_l, Funct, e,ref_period)
    Relaxing(TPhases2, Phases, ref_period)
    data = Phases.reshape([L,L])
    mat.set_data(data)
    #print(data)
    return mat,

# Runs profile on the CMP2D function without animation

pr = cProfile.Profile()
pr.enable()
CMP2D(L,v,d,seed2, tot_time, pace_rate, ref_period)
print('Atrium')
pr.disable()
sortby = 'cumulative'
ps = pstats.Stats(pr).sort_stats(sortby)
ps.print_stats()
"""

# Runs animation
Neigh_u, Neigh_d, Neigh_r, Neigh_l, Phases, Funct  = CreateAtrium(L,v,d,seed2,ref_period)
fig = plt.figure(figsize = [10,10])
ax = plt.subplot()
ax.set_axis_off()
mat = ax.matshow(Phases.reshape([L,L]),cmap=plt.cm.gray_r)
mat.set_clim(0,ref_period)
#plt.colorbar(mat)
ani = animation.FuncAnimation(fig, update, frames = tot_time,fargs = (mat,seed2, pace_rate, ref_period, Neigh_u, Neigh_d, Neigh_r, Neigh_l, Phases, Funct), interval=100, repeat = None)
#plt.axis('auto')
plt.axis([0,L,0,L])
plt.show()
 
"""