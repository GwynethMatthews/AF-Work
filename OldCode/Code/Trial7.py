import numpy as np
import random as rnd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import cProfile
import pstats
from line_profiler import LineProfiler
def Atrium(L,v,d,seed1, seed2,ref_period):  
    N_up = np.full((L*L),fill_value = None, dtype = float)
    N_down = np.full((L*L),fill_value = None, dtype = float)
    N_left = np.full((L*L),fill_value = None, dtype = float)
    N_right = np.full((L*L),fill_value = None, dtype = float)
    Phases = np.ndarray(L*L,dtype = float)
    #Potentials = np.ndarray(L*L,dtype = float)
    tbe = np.ndarray(L*L)
    Phases.fill(ref_period)
    #Potentials.fill(-90)
    Dysfun = np.ndarray([L*L], dtype = bool)
    index = np.indices((1, L*L))[1][0]
    #Atrium  = index.reshape(L,L) # The index for that site within the long arrays
    np.random.seed(seed1)
    w = np.random.rand(L*L)
    np.random.seed(seed2)
    z = np.random.rand(L*L)
    for j in index:
        #z = rnd.uniform(0,1)
        if d > z[j]: # dysfunctional
            Dysfun[j] = True
        if d <= z[j]: # functional
            Dysfun[j] = False
        if j in np.arange(0,L*L,L): # first column
            # right connection for the first column
            N_right[j] = j+1 
        elif j in (np.arange(0,L*L,L)+L-1): # last column
            # left connection for last column
            N_left[j] = j-1
        else: # body columns
            # currently all cells are connected longitudinally
            # left connection for body of cells
            N_left[j] = j-1
            # right connection for body of cells
            N_right[j] = j+1 
    for j in index:
        # transverse connections (checks each one for a downward connection)
        if w[j] <= v: 
            #print(w)
            if j in np.arange(L*L-L,L*L):
                # the jth site has a downwards connection
                N_down[j] = j-(L*L-L)
                # the j-(L*L-L)th site (corresponding column in the first row
                # has an upwards connection
                N_up[j-(L*L-L)] = j
            else:
                # the jth site has a downwards connection
                N_down[j] = j+L
                # the j+Lth site has an upwards connection
                N_up[j+L] = j
    return N_up, N_down, N_right, N_left, Phases, Dysfun,index,tbe
def SinusRhythm(L,e,ref_period,TPhases, Phases,Dysfun):
    """Pacemaker activity"""
    # TAKE THIS OUT
    tbe = np.empty(L*L)
    tbe.fill(False)
    # index of the first column the pacemaker cells (sorta)
    first_col = np.arange(0,L*L,L)
    # gets rid of the non-resting cells
    first_col = first_col[TPhases[first_col]==(ref_period)]
    # picks out the dysfunctional cells
    #DOESN't CHANGE
    dysfunctional_cells = first_col[Dysfun[first_col]]
    # creates an array of random number between 0 and 1 based on a uniform 
    # distribution    
    e_comp_val1 = np.random.rand(len(dysfunctional_cells))
    # picks out dysfunctional cells that fire
    dysfunctional_cells = dysfunctional_cells[e_comp_val1 > e]
    # picks out the functional cells
    functional_cells = first_col[~Dysfun[first_col]]
    # changes the phase of the cells that get excited
    tbe[dysfunctional_cells] = True
    tbe[functional_cells] = True
    # changes the phase of the cells that get excited
    Phases[dysfunctional_cells] = 0
    Phases[functional_cells] = 0
    #Potentials[dysfunctional_cells] = 20
    #Potentials[functional_cells] = 20
    

def Relaxing(index, TPhases, Phases, ref_period):
    """All cells move to the next phase"""
    # finds the index of cells which are not resting i.e in phases 0 to 
    # refractor_period - 1
    non_resting = np.array(index[TPhases != (ref_period)],dtype = int)
    # updates the phases of the non-resting states by 1
    Phases[non_resting] = Phases[non_resting] + 1 
    #Potentials[non_resting] -= 2.2

def Excite(L,index,tbe,TPhases,Phases,e,N_up,N_down,N_right,N_left,Dysfun):
    tbe.fill(False)
    excited = index[TPhases == 0]
    neighbours_up = np.array(N_up[excited][~np.isnan(N_up[excited])],dtype = int)
    #print(neighbours_up)
    if len(neighbours_up) > 0:
        neighbours_up = neighbours_up[TPhases[neighbours_up]==ref_period]
        n_up_dys = neighbours_up[Dysfun[neighbours_up]]   
        e_comp_val2 = np.random.rand(len(n_up_dys))
        n_up_dys = n_up_dys[e_comp_val2 > e]
        n_up_fun = neighbours_up[~Dysfun[neighbours_up]]
        tbe[n_up_dys] = True
        tbe[n_up_fun] = True
    neighbours_down = np.array(N_down[excited][~np.isnan(N_down[excited])],dtype = int)
    #print(neighbours_down)
    if len(neighbours_down) > 0:
        neighbours_down = neighbours_down[TPhases[neighbours_down] == ref_period]
        n_down_dys = neighbours_down[Dysfun[neighbours_down]]   
        e_comp_val3 = np.random.rand(len(n_down_dys))
        n_down_dys = n_down_dys[e_comp_val3 > e]
        n_down_fun = neighbours_down[~Dysfun[neighbours_down]]
        tbe[n_down_dys] = True
        tbe[n_down_fun] = True
    neighbours_right = np.array(N_right[excited][~np.isnan(N_right[excited])],dtype = int)
    #print(neighbours_right)
    if len(neighbours_right) > 0:
        neighbours_right = neighbours_right[TPhases[neighbours_right]==ref_period]
        n_right_dys = neighbours_right[Dysfun[neighbours_right]]   
        e_comp_val4 = np.random.rand(len(n_right_dys))
        n_right_dys = n_right_dys[e_comp_val4 > e]
        n_right_fun = neighbours_right[~Dysfun[neighbours_right]]
        tbe[n_right_dys] = True
        tbe[n_right_fun] = True
    neighbours_left = np.array(N_left[excited][~np.isnan(N_left[excited])],dtype = int)
    #print(neighbours_left)
    if len(neighbours_left) > 0:
        neighbours_left = neighbours_left[TPhases[neighbours_left]==ref_period]
        n_left_dys = neighbours_left[Dysfun[neighbours_left]]   
        e_comp_val5 = np.random.rand(len(n_left_dys))
        n_left_dys = n_left_dys[e_comp_val5 > e]
        n_left_fun = neighbours_left[~Dysfun[neighbours_left]]
        tbe[n_left_dys] = True
        tbe[n_left_fun] = True   
    Phases[tbe == True] = 0
    #GET RID OF POTENTIALS
    #Potentials[tbe == True] = 20
    #np.place(Phases,tbe,0)
    #Phases = np.putmask(Phases,tbe==1,0)
    #Potentials = np.putmask(Potentials,tbe==1,20)
    #print(Phases.reshape(L,L))
        
    
def CMP2D(seed1,seed2,seed3,L,v,d,e,ref_period,pace_rate,tot_time):
    # creates the array
    N_up, N_down, N_right, N_left, Phases,Dysfun, index, tbe= Atrium(L,v,d,seed1, seed2,ref_period)
    t = 0
    np.random.seed(seed3)
    tbe = np.empty(L*L)
    while t < tot_time:
        #print(t)
        # runs the time step change
        CMP2d_1step(L,e,ref_period,pace_rate,tbe,Phases,Dysfun,tot_time,t,index,N_up,N_down,N_right,N_left)
        t = t + 1
        #print(Phases.reshape(L,L))    
    
def CMP2d_1step(L,e,ref_period,pace_rate,tbe,Phases,Dysfun,tot_time,t,index,N_up,N_down,N_right,N_left):
    TPhases = Phases.copy()
    #np.random.seed(seed2)
    #print(Phases.reshape(L,L))
    if np.remainder(t, pace_rate) == 0:
        SinusRhythm(L,e,ref_period,TPhases, Phases,Dysfun)
        #print(Phases.reshape(L,L))
    #if np.all(TPhases1) == False:
    Excite(L,index,tbe,TPhases,Phases,e,N_up,N_down,N_right,N_left,Dysfun)
    #print(Phases.reshape(L,L))
    Relaxing(index, TPhases, Phases, ref_period)
    #print(Phases.reshape(L,L))
    #return Phases

def update(frame_number, mat,seed3,L,e,ref_period,pace_rate,tbe,Phases,Dysfun,tot_time,index,N_up,N_down,N_right,N_left):
    TPhases = Phases.copy()
    #np.random.seed(seed3)
    if frame_number in pace_rate:
        SinusRhythm(L,e,ref_period,TPhases, Phases,Dysfun)
        #print(Phases.reshape(L,L))
    #if np.all(TPhases1) == False:
    Excite(L,index,tbe, TPhases,Phases,e,N_up,N_down,N_right,N_left,Dysfun)
    #print(Phases.reshape(L,L))
    Relaxing(index, TPhases, Phases, ref_period)
    data = Phases.reshape([L,L])
    mat.set_data(data)
    #print(data)
    return mat,

L = 200
v = 0.1
d = 0.05
e = 0.05
seed1 = 1
seed2 = 2
seed3 = 3
ref_period = 51
tot_time = 10**6
pace_rate = 220 #np.arange(0,tot_time,220)
#N_up, N_down, N_right, N_left, Phases, Dysfun, index= Atrium(L,v,d,seed1, seed2,ref_period)
"""
pr = cProfile.Profile()
pr.enable()
CMP2D(seed1,seed2,seed3,L,v,d,e,ref_period,pace_rate,tot_time)
print('Atrium')
pr.disable()
sortby = 'cumulative'
ps = pstats.Stats(pr).sort_stats(sortby)
ps.print_stats()

#CMP2D(seed1,seed2,seed3,L,v,d,e,ref_period,pace_rate,tot_time)
lp = LineProfiler()
lp_wrapper = lp(CMP2D)
lp.add_function(CMP2d_1step)
lp.add_function(SinusRhythm)
lp.add_function(Excite)
lp.add_function(Relaxing)
lp_wrapper(seed1,seed2,seed3,L,v,d,e,ref_period,pace_rate,tot_time)
lp.print_stats()
"""
L = 200
v = 0.1
d = 0.05
e = 0.05
seed1 = 1
seed2 = 2
seed3 = 3
ref_period = 51
tot_time = 10**5
pace_rate = np.arange(0,tot_time,220)
N_up, N_down, N_right, N_left, Phases, Dysfun, index,tbe= Atrium(L,v,d,seed1, seed2,ref_period)
fig = plt.figure(figsize = [15,15])
ax = plt.subplot()
ax.set_axis_off()
mat = ax.matshow(Phases.reshape([L,L]),cmap=plt.cm.gray_r)
mat.set_clim(0,ref_period)
#plt.colorbar(mat)
ani = animation.FuncAnimation(fig, update, frames = tot_time,fargs = (mat,seed3,L,e,ref_period,pace_rate,tbe,Phases,Dysfun,tot_time,index,N_up,N_down,N_right,N_left), interval=100, repeat = None)
#plt.axis('auto')
plt.axis([0,L,0,L])
plt.show()
