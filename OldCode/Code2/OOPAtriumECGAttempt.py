import numpy as np
import matplotlib.gridspec as mgs
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pickle
class Atrium():
    """Creates the myocardium's structure"""
    def __init__(self,L=200,v=0.17,d=0.05,e=0.05,rp=50,tot_time = 10**4
                 ,pace_rate = 220,seed1 = 1,seed2=2,seed3=4):
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
        self.y = np.indices((self.size,self.size))[0]
        self.x = np.indices((self.size,self.size))[1]
        self.x_rel_probe = self.y - (self.size/2) + 0.5
        self.y_rel_probe = self.x - (self.size/2) + 0.5
        self.n_up = np.full((L*L),fill_value = None,dtype = float)
        self.n_down = np.full((L*L),fill_value = None,dtype = float)
        self.n_right = np.full((L*L),fill_value = None,dtype = float)
        self.n_left = np.full((L*L),fill_value = None,dtype = float)
        self.phases = np.full((L*L),fill_value = self.rp)
        self.time_for_ECG = np.arange(-500,1)
        self.potentials = np.zeros(len(self.time_for_ECG))
        self.V = np.full((L*L),fill_value = -90.0)
        self.dysfunctional_cells = np.full([L*L],fill_value = False
                                           , dtype = bool)
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
        self.first_dys = np.array(self.first_col
                                  [~self.dysfunctional_cells[self.first_col]])
        self.first_fun = np.array(self.first_col
                                  [self.dysfunctional_cells[self.first_col]])
        
        
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
        self.V[self.tbe] = 20.
        self.V[~self.resting] -= 2.2
        self.resting[self.tbe] = False
        self.resting[self.states[-1]] = True
        #self.V[self.states[-1]] = -90.
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
        self.Relaxing_ani()
        self.Conduct()
            
    def CMP2D(self):
        np.random.seed(self.seed_prop)
        num_ex_cells = []
        ecg_values = []
        in_AF = []
        time = np.arange(0,self.tot_time)
        while self.t < self.tot_time:
            self.CMP2D_timestep()
            excited_cells = len(self.states[0])
            num_ex_cells.extend([excited_cells])
            ECG_value = self.ECG([((self.size/2)+1.5),((self.size/2)+1.5)])
            ecg_values.extend([ECG_value])
            #print(excited_cells)
            if excited_cells > self.size*1.1:
                if self.t_SR == 0:
                    self.t_AF += 1
                    in_AF.extend([True])
                else:
                    self.t_AF += 1
                    self.t_SR = 0
                    in_AF.extend([True])
            if excited_cells <= self.size*1.1:
                if self.t_AF > 0:
                    if self.t_SR > self.pace_rate:
                        in_AF.extend([False])
                        self.tot_AF += self.t_AF
                        self.t_AF = 0
                    else: 
                        in_AF.extend([True])
                        self.t_AF += 1
                        self.t_SR += 1  
                else:
                    in_AF.extend([False])
            self.t += 1
        data = [num_ex_cells,ecg_values,in_AF,time]    
        pickle.dump(data,open( "data_page67_0.17.p", "wb" ) )
        #self.tot_AF += self.t_AF
    
    
    def ECG(self,LoP):
        volt = self.V.reshape(self.size,self.size)
        numerator = (((self.x[1:,1:]-LoP[0])*(volt[1:,1:]-volt[1:,:-1])) - 
                     ((self.y[1:,1:]-LoP[1])*(volt[1:,1:]-volt[:-1,1:])))
        denominator = (((self.x[1:,1:]-LoP[0])**2)+
                       ((self.y[1:,1:]-LoP[1])**2))**(3/2)
        values = numerator/denominator
        ECG_value1 = sum(values.flatten())
        return ECG_value1
    
    
A = Atrium()
A.CMP2D()    
    
###############################################################################
#RISK CURVE

#def RiskCurve(repeats):
#    """Collects data for and plots the risk curve"""
#    AF_risks = []
#    AF_risks_std = []
#    v_tAFs = []
#    #nus = [0,0.001,0.01,0.02,0.05,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.3,0.5,1]
#    #nus_K = [0.02,0.04,0.06,0.08,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18
#           #,0.19,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.3,1]
#    #nus = [0,0.005,0.01,0.02,0.05,0.08,0.09,0.1,0.11,0.12,0.13,
#           #0.15,0.17,0.18,0.19,0.2,0.21,0.3,0.5,1]
#    nus = [0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.3]       
#    for j in nus:
#        v_tAF = []
#        for i in range(repeats):
#            v = j
#            print(v)
#            seed1 = np.random.seed(i)
#            seed2 = np.random.seed(i+repeats)
#            seed3 = np.random.seed(i+repeats+repeats)
#            A = Atrium(v=j, seed1=seed1, seed2=seed2, seed3=seed3)
#            A.CMP2D()
#            v_tAF.extend([A.tot_AF/A.tot_time])
#        sample_avg = sum(v_tAF)/repeats
#        sample_std = np.std(v_tAF)
#        v_tAFs.extend([v_tAF])
#        AF_risks.extend([sample_avg])
#        AF_risks_std.extend([sample_std])
#    plt.plot(nus,AF_risks)
#    plt.title('Risk of AF', fontsize = 24)
#    plt.xlabel('Fraction of Transverse Connections',fontsize = 20)
#    plt.ylabel('Time in AF/Risk of AF',fontsize = 20)
#    plt.xlim(0,1)
#    plt.errorbar(nus, AF_risks, yerr=AF_risks_std)
#    data = [nus,AF_risks,v_tAFs]
#    pickle.dump(data,open( "data_risk_curve_extra_data.p", "wb" ) )

#RiskCurve(50)
    
    
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

#def update2(frame_number, mat,A,denominator):
#    """Next frame update for animation with ECG"""
#    if frame_number in A.pace:
#        A.SinusRhythm()
#    A.Relaxing_ani()
#    A.Conduct()
#    data = A.phases.reshape([A.size,A.size])
#    mat.set_data(data)
#    ECG_value = A.ECG(denominator)
#    A.potentials = np.roll(A.potentials,-1)
#    A.potentials[-1] = ECG_value
#    ax2.clear()
#    ax2.plot(A.time_for_ECG,A.potentials)
#    ax2.set_xlim(-500,0) 
#    #ax2.set_ylim(-400,400)
#    ax2.set_xlabel('Time')
#    ax2.set_ylabel('Potential')
#    return mat,
#
def update1(frame_number, mat,A):
    """Next frame update for animation without ECG"""
    if frame_number in A.pace:
        A.SinusRhythm()
    A.Relaxing_ani()
    A.Conduct()
    data = A.phases.reshape([A.size,A.size])
    mat.set_data(data)
    return mat,
###############################################################################
# Without ECG
#A = Atrium()
#
#np.random.rand(A.seed_prop)
#fig1 = plt.figure(figsize = [15,15])
#ax = fig1.subplots(1,1)
#mat1 = ax.matshow(A.phases.reshape([A.size,A.size]),cmap=plt.cm.gray_r)
#mat1.set_clim(0,A.rp)
#ani1 = animation.FuncAnimation(fig1, update1, frames = A.tot_time
#                               ,fargs = (mat1,A), interval=100, repeat = None)
#plt.axis([0,A.size,0,A.size])

###############################################################################
# With ECG
#A = Atrium()
#LoP = [100.5,100.5]
#fig2 = plt.figure()
#mgs.GridSpec(4,5)
#ax1 = plt.subplot2grid((4,5), (0,1), colspan=3, rowspan=3)
#ax2 = plt.subplot2grid((4,5), (3,0), colspan=5)
##ax1.set_axis_off()
#mat2 = ax1.matshow(A.phases.reshape([A.size,A.size]),cmap=plt.cm.gray_r)
#mat2.set_clim(0,A.rp)
#ani2 = animation.FuncAnimation(fig2, update2, frames = A.tot_time
#                    ,fargs = (mat2,A,LoP), interval=100, repeat = None)

###############################################################################


#plt.show()
