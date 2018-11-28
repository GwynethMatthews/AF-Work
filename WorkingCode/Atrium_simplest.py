"""

"""
import numpy as np
import pickle

class Atrium():
    """Creates the myocardium's structure.
    """
    def __init__(self,L=200,v=0.2,d=0.05,e=0.05,rp=50,tot_time = 10**6
                 ,pace_rate = 220,seed1 = 1,seed2=2,seed3=3):
        self.size = L 
        self.first_col = np.arange(0,self.size*self.size,self.size)
        self.transverse_prob = v
        self.dysfunctional_prob = d
        self.nonfire_prob = e
        self.rp =rp
        self.t = 0
        self.tot_time = tot_time
        self.tot_AF = 0 # overall
        self.t_AF = 0 # in this episode
        self.t_SR = 0 # in this period of SR
        self.pace_rate = pace_rate # number of timesteps between sinus beats
        self.pace = np.arange(0,self.tot_time,self.pace_rate) # timesteps of SR beat
        self.seed_dysfunc = seed1
        self.seed_connect = seed2
        self.seed_prop = seed3
        self.index = np.arange(0,L*L) # cell positions in each array
        self.y = np.indices((self.size,self.size))[0] # y coordinate for cells
        self.x = np.indices((self.size,self.size))[1] # x coordinate for cells
        self.n_up = np.full((L*L),fill_value = None,dtype = float) # neighbours above
        self.n_down = np.full((L*L),fill_value = None,dtype = float) # neighbours below
        self.n_right = np.full((L*L),fill_value = None,dtype = float) # neighbours to right
        self.n_left = np.full((L*L),fill_value = None,dtype = float) #neighbours to left
        self.phases = np.full((L*L),fill_value = self.rp) # state cell is in (0 = excited, rp = resting)
        self.time_for_ECG = np.arange(-500,1) # time for ECG plot
        self.potentials = np.zeros(len(self.time_for_ECG)) # ECG plot values
        self.V = np.full((L*L),fill_value = -90.0) # voltage depending on state of cell given in Supplementary Material
        self.dysfunctional_cells = np.full([L*L],fill_value = False
                                           , dtype = bool)
        self.states = [[]]*self.rp # list of lists containing cells in each state except resting
        self.resting = np.full([L*L],fill_value = True, dtype = bool) # can they be excited
        self.tbe = np.full([L*L],fill_value = False, dtype = bool) # cells to be excited in next timestep
        
        #setting connections and dysfubctional cells
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
                    
        # functional and dysfunctional cells for first column (speeds up Sinus Rhythm)
        self.first_dys = np.array(self.first_col
                                  [~self.dysfunctional_cells[self.first_col]])
        self.first_fun = np.array(self.first_col
                                  [self.dysfunctional_cells[self.first_col]])
        
        
    def SinusRhythm(self):
        """Pacemaker activity. Sets first column of cells to true in tbe if excited"""
        e_comp_val1 = np.random.rand(len(self.first_dys))
        dysfunctional_cells = self.first_dys[e_comp_val1 > self.nonfire_prob]
        self.tbe[dysfunctional_cells] = True
        self.tbe[self.first_fun] = True
        
    def Relaxing(self):
        """All cells move to the next phase. tbe cells get excited, states 
        move down the refractory phase until resting"""
        self.resting[self.tbe] = False
        self.resting[self.states[-1]] = True
        del self.states[-1]
        self.states.insert(0,self.index[self.tbe])
        
    def Relaxing_ani(self):
        """All cells move to the next phase. tbe cells get excited, states 
        move down the refractory phase until resting. Includes change to phases
        and voltages"""
        self.phases[self.tbe] = 0
        self.phases[~self.resting] += 1
        self.V[self.tbe] = 20.
        self.V[~self.resting] -= 2.2
        self.resting[self.tbe] = False
        self.resting[self.states[-1]] = True
        del self.states[-1]
        self.states.insert(0,self.index[self.tbe])    
        
    def Conduct(self):
        """Finds neighbours of excited cells and sets their tbe to True"""
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
        # Try changing this to 
        self.tbe[neighbours_fun] = True
        self.tbe[neighbours_dys] = True
        self.tbe[self.states[0]] = False

    def CMP2D_timestep(self):
        """A single timestep"""
        if np.remainder(self.t,self.pace_rate)==0:
            self.SinusRhythm()
        self.Relaxing()
        self.Conduct()
            
    def CMP2D(self):
        """The basic model. CM2D_timestep() runs for tot_time."""
        np.random.seed(self.seed_prop)
        while self.t < self.tot_time:
            self.CMP2D_timestep()
            self.t += 1
    
    def ECG(self,LoP):
        """Calculates the ECG value for a single timestep"""
        volt = self.V.reshape(self.size,self.size)
        numerator = (((self.x[1:,1:]-LoP[0])*(volt[1:,1:]-volt[1:,:-1])) - 
                     ((self.y[1:,1:]-LoP[1])*(volt[1:,1:]-volt[:-1,1:])))
        denominator = (((self.x[1:,1:]-LoP[0])**2)+
                       ((self.y[1:,1:]-LoP[1])**2))**(3/2)
        values = numerator/denominator
        ECG_value1 = sum(values.flatten())
        return ECG_value1
    
    def CMP2D_page67(self):
        """Runs CMP2D() and collects data for ECG and number of excited cells. 
        Used to replicate page67 of Kishan's thesis"""
        np.random.seed(self.seed_prop)
        num_ex_cells = []
        ecg_values = []
        time = np.arange(0,self.tot_time)
        while self.t < self.tot_time:
            self.CMP2D_timestep()
            excited_cells = len(self.states[0])
            num_ex_cells.extend([excited_cells])
            ECG_value = self.ECG([((self.size/2)+0.5),((self.size/2)+0.5)])
            ecg_values.extend([ECG_value])                      
            self.t += 1
        data = [num_ex_cells,ecg_values,time]    
        pickle.dump(data,open( "data_page67.p", "wb" ) )

    def CMP2D_time_AF(self):
        """Runs CMP2D and collects data on the ammount of time in AF"""
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
