import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pickle
class Atrium():
    def __init__(self,L=200,v=0.2,d=0.05,e=0.05,rp=50,tot_time = 10**6,pace_rate = 220,seed1 = 1,seed2=2,seed3=3):
        self.size = L
        self.first_col = np.arange(0,self.size*self.size,self.size)
        self.transverse_prob = v
        self.dysfunctional_prob = d
        self.nonfire_prob = e
        self.rp =rp
        self.t = 0
        self.tot_time = tot_time
        self.tot_AF = 0
        self.t_AF = 0
        self.t_SR = 0
        self.pace_rate = pace_rate
        self.pace = np.arange(0,self.tot_time,self.pace_rate)
        self.seed_dysfunc = seed1
        self.seed_connect = seed2
        self.seed_prop = seed3
        self.index = np.arange(0,L*L)
        self.n_up = np.full((L*L),fill_value = None,dtype = float)
        self.n_down = np.full((L*L),fill_value = None,dtype = float)
        self.n_right = np.full((L*L),fill_value = None,dtype = float)
        self.n_left = np.full((L*L),fill_value = None,dtype = float)
        self.phases = np.full((L*L),fill_value = self.rp)
        self.dysfunctional_cells = np.full([L*L],fill_value = False, dtype = bool)
        self.excited = np.full([L*L],fill_value = False, dtype = bool)
        self.states = [[]]*self.rp
        self.resting = np.full([L*L],fill_value = True, dtype = bool)
        self.tbe = np.full([L*L],fill_value = False, dtype = bool)
        np.random.seed(self.seed_dysfunc)
        w = np.random.rand(L*L)
        np.random.seed(self.seed_connect)
        z = np.random.rand(L*L)
        
        for j in self.index:
            if d > z[j]: # dysfunctional
                self.dysfunctional_cells[j] = False
            if d <= z[j]: # functional
                self.dysfunctional_cells[j] = True
            if j in np.arange(0,L*L,L):
                self.n_right[j] = j+1 
            elif j in (np.arange(0,L*L,L)+L-1):
                self.n_left[j] = j-1
            else:
                self.n_left[j] = j-1
                self.n_right[j] = j+1 
                
        for j in self.index:
            if w[j] <= v: 
                if j in np.arange(L*L-L,L*L):
                    self.n_down[j] = j-(L*L-L)
                    self.n_up[j-(L*L-L)] = j
                else:
                    self.n_down[j] = j+L
                    self.n_up[j+L] = j
        self.first_dys = np.array(self.first_col[~self.dysfunctional_cells[self.first_col]])
        self.first_fun = np.array(self.first_col[self.dysfunctional_cells[self.first_col]])
        
        
    def SinusRhythm(self):
        """Pacemaker activity"""
        e_comp_val1 = np.random.rand(len(self.first_dys))
        dysfunctional_cells = self.first_dys[e_comp_val1 > self.nonfire_prob]
        self.tbe[dysfunctional_cells] = True
        self.tbe[self.first_fun] = True
        
    def Relaxing(self):
        """All cells move to the next phase"""
        self.resting[self.tbe] = False
        self.resting[self.states[-1]] = True
        del self.states[-1]
        self.states.insert(0,self.index[self.tbe])
        
    def Relaxing_ani(self):
        """All cells move to the next phase"""
        self.phases[self.tbe] = 0
        self.phases[~self.resting] += 1
        self.resting[self.tbe] = False
        self.resting[self.states[-1]] = True
        del self.states[-1]
        self.states.insert(0,self.index[self.tbe])    
        
    def Conduct(self):
        """Finds neighbours of excited cells and excites them"""
        neighbours_up = self.n_up[self.states[0][~np.isnan(self.n_up[self.states[0]])],]
        neighbours_down = self.n_down[self.states[0][~np.isnan(self.n_down[self.states[0]])]]
        neighbours_right = self.n_right[self.states[0][~np.isnan(self.n_right[self.states[0]])]]
        neighbours_left = self.n_left[self.states[0][~np.isnan(self.n_left[self.states[0]])]]
        neighbours = [neighbours_up,neighbours_down,neighbours_left,neighbours_right]
        neighbours = np.array(np.concatenate(neighbours),dtype = int) 
        neighbours = neighbours[self.resting[neighbours]]
        neighbours_dys = neighbours[~self.dysfunctional_cells[neighbours]]
        e_comp_val2 = np.random.rand(len(neighbours_dys))
        neighbours_dys = neighbours_dys[e_comp_val2 > self.nonfire_prob]
        neighbours_fun = neighbours[self.dysfunctional_cells[neighbours]]
        self.tbe[neighbours_fun] = True
        self.tbe[neighbours_dys] = True
        self.tbe[self.states[0]] = False

    def CMP2D_timestep(self):
        if self.t in self.pace:
            self.SinusRhythm()
        self.Relaxing()
        self.Conduct()
            
    def CMP2D(self):
        np.random.seed(self.seed_prop)
        while self.t < self.tot_time:
            self.CMP2D_timestep()
            excited_cells = len(self.states[0])            
            if excited_cells > self.size*1.1:
                if self.t_SR == 0:
                    self.t_AF += 1
                else:
                    self.t_AF += 1
                    self.t_SR = 0
            if excited_cells <= self.size*1.1:
                if self.t_AF > 0:
                    if self.t_SR > self.pace_rate:
                        self.tot_AF += self.t_AF
                        self.t_AF = 0
                    else: 
                        self.t_AF += 1
                        self.t_SR += 1                       
            self.t += 1
        self.tot_AF += self.t_AF
    
        
###############################################################################
#RISK CURVE
#
def RiskCurve():
    """Collects data for and plots the risk curve"""
    AF_risks = []
    AF_risks_std = []
    v_tAFs = []
    #nus = [0,0.001,0.01,0.02,0.05,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.3,0.5,1]
    nus = [0,0.001,0.005,0.01,0.02,0.05,0.08,0.09,0.1,0.11,0.115,0.12,0.125,0.13,0.135,0.14,0.145,0.15,0.155,0.16,0.165,0.17,0.175,0.18,0.185,0.19,0.195,0.2,0.21,0.22,0.3,0.5,1]
    for j in nus:
        v_tAF = []
        if j < 0.12:
            repeats = 20
        if j >= 0.12 and j < 0.3:
            repeats = 50
        if j >= 0.3:
            repeats = 10            
        for i in range(repeats):
            v = j
            print(v)
            seed1 = np.random.seed(i)
            seed2 = np.random.seed(i+repeats)
            seed3 = np.random.seed(i+repeats+repeats)
            A = Atrium(v=j, seed1=seed1, seed2=seed2, seed3=seed3)
            A.CMP2D()
            v_tAF.extend([A.tot_AF/A.tot_time])
        sample_avg = sum(v_tAF)/repeats
        sample_std = np.std(v_tAF)
        v_tAFs.extend([v_tAF])
        AF_risks.extend([sample_avg])
        AF_risks_std.extend([sample_std])
    plt.plot(nus,AF_risks)
    plt.title('Risk of AF', fontsize = 24)
    plt.xlabel('Fraction of Transverse Connections',fontsize = 20)
    plt.ylabel('Time in AF/Risk of AF',fontsize = 20)
    plt.xlim(0,1)
    plt.errorbar(nus, AF_risks, yerr=AF_risks_std)
    print(v_tAFs)
    data = [nus,AF_risks,v_tAFs]
    pickle.dump(data,open( "data_risk_curve_dif_repeats.p", "wb" ) )
    
"""multi-processor"""
"""KIM FIRE PAPER"""
    

RiskCurve()
###############################################################################
#Line Profiler
#from line_profiler import LineProfiler            
#A= Atrium()        
#lp = LineProfiler()
#lp_wrapper = lp(A.CMP2D)
#lp.add_function(A.CMP2D_timestep)
#lp.add_function(A.SinusRhythm)
#lp.add_function(A.Conduct)
#lp.add_function(A.Relaxing)
#lp_wrapper()
#lp.print_stats()


###############################################################################
# Animation
#Atrium = Atrium()
#def update(frame_number, mat,Atrium):
#    """Next frame update for animation"""
#    if frame_number in Atrium.pace:
#        Atrium.SinusRhythm()
#    Atrium.Relaxing_ani()
#    Atrium.Conduct()
#    #Atrium.time = frame_number
#    #print(Atrium.time)
#    #Atrium.CMP2D_timestep()
#    data = Atrium.phases.reshape([Atrium.size,Atrium.size])
#    mat.set_data(data)
#    return mat,
#np.random.rand(Atrium.seed_prop)
#fig = plt.figure(figsize = [15,15])
#ax = plt.subplot()
#ax.set_axis_off()
#mat = ax.matshow(Atrium.phases.reshape([Atrium.size,Atrium.size]),cmap=plt.cm.gray_r)
#mat.set_clim(0,Atrium.rp)
#ani = animation.FuncAnimation(fig, update, frames = Atrium.tot_time,fargs = (mat,Atrium), interval=100, repeat = None)
#plt.axis([0,Atrium.size,0,Atrium.size])
#plt.show()