import numpy as np
import pickle

class Atrium:
    """Creates the myocardium's structure.
    """
    def __init__(self, hexagonal = False, L = 200, v_para = 0.6, v_tran_1 = 0.6,
                 v_tran_2 = 0.6, d = 0.05, e = 0.05, rp = 50, tot_time = 10**6, 
                 pace_rate = 220, seed1 = 10, seed2 = 20, seed3 = 30, seed4 = 40):
        
        # System Parameters
        self.hexagonal = hexagonal
        self.size = L 
        
        self.parallel_prob = v_para
      
        if self.hexagonal == False:
            self.transverse_prob = v_tran_1
            
        if self.hexagonal == True:
            self.transverse_prob_l = v_tran_1
            self.transverse_prob_r = v_tran_2
        self.tot_time = tot_time
        self.dysfunctional_prob = d
        self.nonfire_prob = e
        self.rp = rp        
        self.pace_rate = pace_rate
        self.pace = np.arange(0,self.tot_time,self.pace_rate)
        
        # System cell positions
        self.first_col = np.arange(0,self.size*self.size,self.size)
        self.index = np.arange(0,L*L) # cell positions in each array
        self.position = self.index.reshape(self.size,self.size)
        self.y = np.indices((self.size,self.size))[0] # y coordinate for cells
        self.x = np.indices((self.size,self.size))[1] # x coordinate for cells 
        self.dysfunctional_cells = np.empty([L*L], dtype = bool) 
        self.dysfunctional_cells.fill(False)
        
        #System seeds
        self.seed_dysfunc = seed1
        self.seed_connect_tran = seed2
        self.seed_connect_para = seed3
        self.seed_prop = seed4

        # Measurable Variables
        self.t = 0
        self.tot_AF = 0 # overall
        self.t_AF = 0 # in this episode
        self.t_SR = 0 # in this period of SR

        # Neighbours
        if self.hexagonal == False:
            self.neighbours = np.empty(((L*L)*4),dtype = float)
            self.neighbours.fill(None)
            self.start_n_right = self.size*self.size
            self.start_n_down = self.size*self.size*2
            self.start_n_left = self.size*self.size*3
            
        if self.hexagonal == True:
            self.neighbours = np.empty(((L*L)*6),dtype = float)
            self.neighbours.fill(None)
            self.start_n_up_left = 0
            self.start_n_up_right = self.size*self.size
            self.start_n_right = self.size*self.size*2
            self.start_n_down_right = self.size*self.size*3
            self.start_n_down_left = self.size*self.size*4
            self.start_n_left = self.size*self.size*5
        

        # For ECG
        self.time_for_ECG = np.arange(-500,1) # time for ECG plot
        self.potentials = np.zeros(len(self.time_for_ECG)) # ECG plot values
        
        # State of system
        self.phases = np.empty((L*L)) # state cell is in (0 = excited, rp = resting)
        self.phases.fill(self.rp)
        self.V = np.empty((L*L)) # voltage depending on state of cell given in Supplementary Material
        self.V.fill(-90.0)
        self.states = [[]]*self.rp # list of lists containing cells in each state except resting
        self.resting = np.empty([L*L],dtype = bool) # can they be excited
        self.resting.fill(True)
        self.tbe = np.empty([L*L], dtype = bool) # cells to be excited in next timestep
        self.tbe.fill(False)
        # Setting connections and dysfunctional cells
        if self.hexagonal == False:
            np.random.seed(self.seed_dysfunc)
            num_rand_dysfunc = np.random.rand(L*L)
            np.random.seed(self.seed_connect_tran)
            num_rand_tran = np.random.rand(L*L)
            np.random.seed(self.seed_connect_para)
            num_rand_para = np.random.rand(L*L)
    
            for j in self.index:
                if self.dysfunctional_prob > num_rand_dysfunc[j]: # dysfunctional
                    self.dysfunctional_cells[j] = False
                if self.dysfunctional_prob <= num_rand_dysfunc[j]: # functional
                    self.dysfunctional_cells[j] = True
    
            for j in self.index:
                if num_rand_para[j] <= self.parallel_prob:
                    if j in np.arange(0,self.size*self.size,self.size):
                        self.neighbours[j+self.start_n_right] = int(j+1)
                        self.neighbours[j+1+self.start_n_left] = int(j)
                        
                    elif j in (np.arange(0,self.size*self.size,self.size)+L-1):
                        self.neighbours[j+self.start_n_right] = None
                        
                    else:
                        self.neighbours[j+self.start_n_right] = int(j+1)             
                        self.neighbours[j+1+self.start_n_left] = int(j)
                        
            for j in self.index:
                if num_rand_tran[j] <= self.transverse_prob: 
                    if j in np.arange(self.size*self.size-self.size,self.size*self.size):
                        self.neighbours[j+self.start_n_down] = j-(self.size*self.size-self.size)
                        self.neighbours[j-(self.size*self.size-self.size)] = j
                    
                    else:
                        self.neighbours[j+self.start_n_down] = j+self.size
                        self.neighbours[j+self.size] = j
        
        
        if self.hexagonal == True:
            np.random.seed(self.seed_dysfunc)
            num_rand_dysfunc = np.random.rand(L*L)
            np.random.seed(self.seed_connect_tran)
            num_rand_tran1 = np.random.rand(L*L)
            num_rand_tran2 = np.random.rand(self.size*self.size)
            np.random.seed(self.seed_connect_para)
            num_rand_para = np.random.rand(L*L)
            
            for j in self.index:
                if d > num_rand_dysfunc[j]: # dysfunctional
                    self.dysfunctional_cells[j] = False
                if d <= num_rand_dysfunc[j]: # functional
                    self.dysfunctional_cells[j] = True
                
                if num_rand_para[j] <= self.parallel_prob:
                    if j in (np.arange(0,self.size*self.size,self.size)+self.size-1):
                        self.neighbours[j+self.start_n_right] = None
                    else:
                        self.neighbours[j+self.start_n_right] = int(j+1)
                        self.neighbours[j+1+self.start_n_left] = int(j)
            
            for j in self.index:
                #even
                if j in self.position[np.arange(0,self.size,2)]:
                    if num_rand_tran1[j] <= self.transverse_prob_l:
                        if j not in self.first_col:
                            self.neighbours[j+self.start_n_down_left] = j+L-1
                            self.neighbours[j+L-1+ self.start_n_up_right] = j
                            
                    if num_rand_tran2 [j] <= self.transverse_prob_r:
                        self.neighbours[j+self.start_n_down_right] = j+L
                        self.neighbours[j+L+self.start_n_up_left] = j
                        
                #odd
                if j in self.position[np.arange(1,self.size,2)]:
                    if j in np.arange(L*L-L,L*L):
                        if num_rand_tran1[j] <= self.transverse_prob_l:
                            self.neighbours[j+self.start_n_down_left] = j-(L*L-L)
                            self.neighbours[j-(L*L-L)+self.start_n_up_right] = j

                    else:
                        if num_rand_tran1[j] <= self.transverse_prob_l:
                            self.neighbours[j+self.start_n_down_left] = j+L
                            self.neighbours[j+L+self.start_n_up_right] = j
                            
                        if num_rand_tran2 [j] <= self.transverse_prob_r:
                            if j not in (np.arange(0,self.size*self.size,self.size)+self.size-1):
                                self.neighbours[j+self.start_n_down_right] = j+L+1
                                self.neighbours[j+L+1+self.start_n_up_left] = j
                                
        # Functional and dysfunctional cells for first column (speeds up Sinus Rhythm)
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
        if self.hexagonal == False:
            neighbours = self.neighbours[np.array([self.states[0],self.states[0]+self.start_n_down,self.states[0]+self.start_n_left,self.states[0]+self.start_n_right])]
            neighbours = np.array(neighbours[~np.isnan(neighbours)],dtype = int) 

        if self.hexagonal == True:
            neighbours = self.neighbours[np.array([self.states[0]+self.start_n_up_left,self.states[0]+self.start_n_up_right,self.states[0]+self.start_n_right,self.states[0]+self.start_n_down_right,self.states[0]+self.start_n_down_left,self.states[0]+self.start_n_left])]
            neighbours = np.array(neighbours[~np.isnan(neighbours)],dtype = int) 
        
        neighbours = neighbours[self.resting[neighbours]]
        neighbours_dys = neighbours[~self.dysfunctional_cells[neighbours]]
        
        e_comp_val2 = np.random.rand(len(neighbours_dys))
        neighbours_dys = neighbours_dys[e_comp_val2 > self.nonfire_prob]
        neighbours_fun = neighbours[self.dysfunctional_cells[neighbours]]

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
        
        
def RiskCurve(repeats,nus_trans,nus_para,seed1,seed2,seed3,seed4):
    """Collects data for and plots the risk curve. Saves data."""
    AF_risks = []
    AF_risks_std = []
    v_tAFs = []
    nus_para = nus_para
    nus = nus_trans 
    for j in nus:
        v_tAF = []
        for i in range(repeats):
            A = Atrium(hexagonal = True, v_para=nus_para,v_tran_1= j, seed1=seed1, seed2=seed2, seed3=seed3,seed4 = seed4)
            A.CMP2D_time_AF()
            v_tAF.extend([A.tot_AF/A.tot_time])
            seed1 +=20
            seed2 +=20
            seed3 +=20
            seed4 +=20
        sample_avg = sum(v_tAF)/repeats
        sample_std = np.std(v_tAF)
        v_tAFs.extend([v_tAF])
        AF_risks.extend([sample_avg])
        AF_risks_std.extend([sample_std])
    data = [nus_para,nus,AF_risks,v_tAFs]
    return data


def Risk_Plane(repeats, nus_para):
    data_full = []
    seed1 = 1
    seed2 = 2
    seed3 = 3
    seed4 = 4
    for k in nus_para:
        seed1 += 200
        seed2 += 200
        seed3 += 200
        seed4 += 200
        nus_trans = [0,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1]
        data = RiskCurve(repeats,nus_trans,k,seed1,seed2,seed3,seed4)
        data_full.extend([data])
    pickle.dump(data_full,open( "risk_plane.p", "wb" ) )
    
repeats = 1
nus_para = [0,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1]
Risk_Plane(repeats,nus_para)