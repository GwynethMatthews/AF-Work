import numpy as np
import random as rnd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import cProfile
import pstats

from line_profiler import LineProfiler
def Atrium(L,v,d,seed1, seed2,rp):  
    """Creates Atrium"""
    N_up = np.full((L*L),fill_value = None, dtype = float) ### change to int
    N_down = np.full((L*L),fill_value = None, dtype = float)
    N_left = np.full((L*L),fill_value = None, dtype = float)
    N_right = np.full((L*L),fill_value = None, dtype = float)
    Phases = np.ndarray(L*L,dtype = float)
    #Potentials = np.ndarray(L*L,dtype = float)
    tbe = np.ndarray(L*L) ## get rid off
    Phases.fill(rp)
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


def CMP2D(seed1,seed2,seed3,L,v,d,e,rp,pace,time):
    """Complete model from creation of grid to end of total_time"""
    N_up, N_down, N_right, N_left, Phases,Dysfun, index, tbe= Atrium(L,v,d,seed1, seed2,rp)
    t = 0
    t_AF = 0
    np.random.seed(seed3)
    tbe = np.empty(L*L,dtype = bool)
    first_col = np.arange(0,L*L,L)
    fc_func = first_col[~Dysfun[first_col]]  
    fc_dys = first_col[Dysfun[first_col]]
    nrr = np.random.rand
    while t < time:
        TPhases = Phases.copy()
        tbe.fill(False)
        if np.remainder(t, pace) == 0:
            SinusRhythm(L,e,rp,nrr,TPhases, Phases,Dysfun,fc_func,fc_dys,tbe)
        Conduct(L,nrr,index,tbe,TPhases,Phases,e,N_up,N_down,N_right,N_left,Dysfun)
        Relaxing(index, TPhases, Phases, rp,tbe)
#        count = (Phases==0).sum()
#        if count >= 220:
#            t_AF += 1
#        if count < 220:
#            if t_AF <= 1:
#                t_AF = 0
        t = t + 1
    return t_AF

def update(frame_number, mat,seed3,L,e,rp,pace_rate,nrr,tbe,Phases,Dysfun,fc_func,fc_dys,time,index,N_up,N_down,N_right,N_left):
    """Next frame update for animation"""
    TPhases = Phases.copy()
    tbe.fill(False)
    if frame_number in pace_rate:
        SinusRhythm(L,e,rp,nrr,TPhases, Phases,Dysfun,fc_func,fc_dys,tbe)
    Conduct(L,nrr,index,tbe, TPhases,Phases,e,N_up,N_down,N_right,N_left,Dysfun)
    Relaxing(index, TPhases, Phases, rp,tbe)
    data = Phases.reshape([L,L])
    mat.set_data(data)
    return mat,

def SinusRhythm(L,e,rp,nrr,TPhases, Phases,Dysfun,fc_func,fc_dys,tbe):
    """Pacemaker activity"""
    e_comp_val1 = nrr(len(fc_dys))
    dysfunctional_cells = fc_dys[e_comp_val1 > e]
    tbe[dysfunctional_cells] = True
    tbe[fc_func] = True

def Relaxing(index, TPhases, Phases, rp,tbe):
    """All cells move to the next phase"""
    tbe = np.array(tbe, dtype= bool) 
    Phases[tbe] = 0
    non_resting = np.array(index[TPhases != (rp)],dtype = int)
    Phases[non_resting] += 1

def Conduct(L,nrr,index,tbe,TPhases,Phases,e,N_up,N_down,N_right,N_left,Dysfun):
    """Finds neighbours of excited cells and excites them"""
    excited = index[TPhases == 0] ## bring tbe from last timestep
    neighbours_up = np.array(N_up[excited][~np.isnan(N_up[excited])],dtype = int)
    if len(neighbours_up) > 0:
        neighbours_up = neighbours_up[TPhases[neighbours_up]==rp]
        n_up_dys = neighbours_up[Dysfun[neighbours_up]]   
        e_comp_val2 = nrr(len(n_up_dys))
        n_up_dys = n_up_dys[e_comp_val2 > e]
        n_up_fun = neighbours_up[~Dysfun[neighbours_up]]
        tbe[n_up_dys] = True
        tbe[n_up_fun] = True
    neighbours_down = np.array(N_down[excited][~np.isnan(N_down[excited])],dtype = int)
    if len(neighbours_down) > 0:
        neighbours_down = neighbours_down[TPhases[neighbours_down] == rp]
        n_down_dys = neighbours_down[Dysfun[neighbours_down]]   
        e_comp_val3 = nrr(len(n_down_dys))
        n_down_dys = n_down_dys[e_comp_val3 > e]
        n_down_fun = neighbours_down[~Dysfun[neighbours_down]]
        tbe[n_down_dys] = True
        tbe[n_down_fun] = True
    neighbours_right = np.array(N_right[excited][~np.isnan(N_right[excited])],dtype = int)
    #print(neighbours_right)
    if len(neighbours_right) > 0:
        neighbours_right = neighbours_right[TPhases[neighbours_right]==rp]
        n_right_dys = neighbours_right[Dysfun[neighbours_right]]   
        e_comp_val4 = nrr(len(n_right_dys))
        n_right_dys = n_right_dys[e_comp_val4 > e]
        n_right_fun = neighbours_right[~Dysfun[neighbours_right]]
        tbe[n_right_dys] = True
        tbe[n_right_fun] = True
    neighbours_left = np.array(N_left[excited][~np.isnan(N_left[excited])],dtype = int)
    #print(neighbours_left)
    if len(neighbours_left) > 0:
        neighbours_left = neighbours_left[TPhases[neighbours_left]==rp]
        n_left_dys = neighbours_left[Dysfun[neighbours_left]]   
        e_comp_val5 = nrr(len(n_left_dys))
        n_left_dys = n_left_dys[e_comp_val5 > e]
        n_left_fun = neighbours_left[~Dysfun[neighbours_left]]
        tbe[n_left_dys] = True
        tbe[n_left_fun] = True     

def RiskCurve(L,d,e,rp,pace,time):
    """Collects data for and plots the risk curve"""
    AF_risks = []
    AF_risks_std = []
    nus = [0,0.01,0.05,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.3,0.5,1]
    for j in nus:
        v_tAF = []
        for i in range(10):
            v = j
            print(v)
            seed1 = np.random.seed(i)
            seed2 = np.random.seed(i+5)
            seed3 = np.random.seed(i+10)
            t_AF = CMP2D(seed1,seed2,seed3,L,v,d,e,rp,pace,time)
            v_tAF.extend([t_AF/time])
        sample_avg = sum(v_tAF)/10
        sample_std = np.std(v_tAF)
        AF_risks.extend([sample_avg])
        AF_risks_std.extend([sample_std])
    plt.plot(nus,AF_risks)
    plt.title('Risk of AF', fontsize = 24)
    plt.xlabel('Fraction of Transverse Connections',fontsize = 20)
    plt.ylabel('Time in AF/Risk of AF',fontsize = 20)
    plt.xlim(0,1)
    plt.errorbar(nus, AF_risks, yerr=AF_risks_std)
        
        
            
    

L = 200
v = 0.19
d = 0.05
e = 0.05
seed1 = 1
seed2 = 2
seed3 = 3
rp = 50
time = 10**6
pace = 220
#RiskCurve(L,d,e,rp,pace,time)

"""
pr = cProfile.Profile()
pr.enable()
t_AF = CMP2D(seed1,seed2,seed3,L,v,d,e,rp,pace,time)
print('Atrium')
print(t_AF)
pr.disable()
sortby = 'cumulative'
ps = pstats.Stats(pr).sort_stats(sortby)
ps.print_stats()
"""
#CMP2D(seed1,seed2,seed3,L,v,d,e,rp,pace,time)
lp = LineProfiler()
lp_wrapper = lp(CMP2D)
#lp.add_function(CMP2d_1step)
lp.add_function(SinusRhythm)
lp.add_function(Conduct)
lp.add_function(Relaxing)
lp_wrapper(seed1,seed2,seed3,L,v,d,e,rp,pace,time)
lp.print_stats()
"""
L = 200
v = 0.1
d = 0.05
e = 0.05
seed1 = 1
seed2 = 2
seed3 = 3
rp = 51
time = 10**5
pace_rate = np.arange(0,time,220)
N_up, N_down, N_right, N_left, Phases, Dysfun, index,tbe= Atrium(L,v,d,seed1, seed2,rp)
tbe = np.empty(L*L)
first_col = np.arange(0,L*L,L)
fc_func = first_col[~Dysfun[first_col]]  
fc_dys = first_col[Dysfun[first_col]]
nrr = np.random.rand

fig = plt.figure(figsize = [15,15])
ax = plt.subplot()
ax.set_axis_off()
mat = ax.matshow(Phases.reshape([L,L]),cmap=plt.cm.gray_r)
mat.set_clim(0,rp)
ani = animation.FuncAnimation(fig, update, frames = time,fargs = (mat,seed3,L,e,rp,pace_rate,nrr,tbe,Phases,Dysfun,fc_func,fc_dys,time,index,N_up,N_down,N_right,N_left), interval=100, repeat = None)
plt.axis([0,L,0,L])
plt.show()
"""