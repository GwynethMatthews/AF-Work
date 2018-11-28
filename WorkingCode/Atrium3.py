import numpy as np
import pickle
"""Atrium3 is the excites if the sum of all the 1/N values is greater than the threshold"""
class Atrium():
    """Creates the myocardium's structure.
    """
    def __init__(self, hexagonal = False, L = 200, v_para = 1, v_tran_1 = 1,
                 v_tran_2 = 0.6, threshold = 0.1, rp = 50, tot_time = 10**6, 
                 pace_rate = 220, seed1 = 10, seed2 = 20, seed3 = 30):
        
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
        self.rp = rp        
        self.pace_rate = pace_rate
        self.pace = np.arange(0,self.tot_time,self.pace_rate)
        self.threshold = threshold
        
        # System cell positions
        self.first_col = np.arange(0,self.size*self.size,self.size)
        self.last_col = np.arange(0,self.size*self.size,self.size)+L-1
        self.index = np.arange(0,L*L) # cell positions in each array
        self.position = self.index.reshape(self.size,self.size)
        self.y = np.indices((self.size,self.size))[0] # y coordinate for cells
        self.x = np.indices((self.size,self.size))[1] # x coordinate for cells 
        
        #System seeds
        self.seed_connect_tran = seed1
        self.seed_connect_para = seed2
        self.seed_prop = seed3

        # Measurable Variables
        self.excitations = np.zeros(self.size*self.size)
        self.t = 0
        self.tot_AF = 0 # overall
        self.t_AF = 0 # in this episode
        self.t_SR = 0 # in this period of SR

        # For ECG
        self.time_for_ECG = np.arange(-500,1) # time for ECG plot
        self.potentials = np.zeros(len(self.time_for_ECG)) # ECG plot values
        
        # State of system
        self.phases = np.full((L*L),fill_value = self.rp) # state cell is in (0 = excited, rp = resting)
        self.V = np.full((L*L),fill_value = -90.0) # voltage depending on state of cell given in Supplementary Material
        self.states = [[]]*self.rp # list of lists containing cells in each state except resting
        self.resting = np.full([L*L],fill_value = True, dtype = bool) # can they be excited
        self.tbe = np.full([L*L],fill_value = False, dtype = bool) # cells to be excited in next timestep
        
        # Setting connections and dysfunctional cells
        if self.hexagonal == False:
            
            # setting up the neighbours
            self.neighbours = np.full(((L*L)*4),fill_value = None,dtype = float)
            self.start_n_right = self.size*self.size
            self.start_n_down = self.size*self.size*2
            self.start_n_left = self.size*self.size*3
            
            # setting up the seeds
            np.random.seed(self.seed_connect_tran)
            num_rand_tran = np.random.rand(L*L)
            np.random.seed(self.seed_connect_para)
            num_rand_para = np.random.rand(L*L)
    
            # PROBABLY A QUICKER WAY OF DOING THIS SINCE DON'T NEED ANYTHING 
            # EXCEPT THE self.list_of_neighbours BUT THAT USES self.neighbours
            for j in self.index:
                if num_rand_para[j] <= self.parallel_prob:
                    if j in self.first_col:
                        self.neighbours[j+self.start_n_right] = int(j+1)
                        self.neighbours[j+1+self.start_n_left] = int(j)
                    
                    elif j in self.last_col:
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
            
            # list of all neighbours for Sq model (includes None's)    
            self.list_of_neighbours =  [[self.neighbours[i],
                             self.neighbours[i+self.start_n_right],
                             self.neighbours[i+self.start_n_left],
                             self.neighbours[i+self.start_n_down]] for i in self.index]    
        
        
        if self.hexagonal == True:
            
            # setting up the neighbours
            self.neighbours = np.full(((L*L)*6),fill_value = None,dtype = float)
            self.start_n_up_left = 0
            self.start_n_up_right = self.size*self.size
            self.start_n_right = self.size*self.size*2
            self.start_n_down_right = self.size*self.size*3
            self.start_n_down_left = self.size*self.size*4
            self.start_n_left = self.size*self.size*5
            
            # setting up the seeds
            np.random.seed(self.seed_connect_tran)
            num_rand_tran1 = np.random.rand(L*L)
            num_rand_tran2 = np.random.rand(self.size*self.size)
            np.random.seed(self.seed_connect_para)
            num_rand_para = np.random.rand(L*L)
            
            # PROBABLY A QUICKER WAY OF DOING THIS SINCE DON'T NEED ANYTHING 
            # EXCEPT THE self.list_of_neighbours BUT THAT USES self.neighbours
            for j in self.index:

                if num_rand_para[j] <= self.parallel_prob:
                    if j in self.last_col:
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
                            if j not in self.last_col:
                                self.neighbours[j+self.start_n_down_right] = j+L+1
                                self.neighbours[j+L+1+self.start_n_up_left] = j
            
            # list of all neighbours for Hex model (includes None's)                     
            self.list_of_neighbours =  [[self.neighbours[i+self.start_n_right],
                             self.neighbours[i+self.start_n_down_right],
                             self.neighbours[i+self.start_n_up_left],
                             self.neighbours[i+self.start_n_up_right],
                             self.neighbours[i+self.start_n_down_left],
                             self.neighbours[i+self.start_n_left]] for i in self.index]  
        
        
        #list of lists of neighbours for each cell (length L*L, removes all None's)
        self.neighbour_list = np.array([[x for x in 
                                         self.list_of_neighbours[i] if str(x) != 'nan'] for i in self.index])

        
    def SinusRhythm(self):
        """Pacemaker activity. Assumes all first column cells get excited """
        self.tbe[self.first_col] = True

        
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
        """Finds neighbours of excited cells and sets their tbe to True if over the threshold"""
        
        # neighbours_list (different from self.neighbour_list) is a list of 
        # lists of all the resting cells of the excited cells
        neighbours_list = [[x for x in self.neighbour_list[i] if self.resting[int(x)]==True] for i in self.states[0]]
        
        # resting neighbours is the number of resting cells connected to each 
        # excited cell (length is the same as that of self.states[0])
        resting_neighbours = list(map(len,neighbours_list))
        
        # inward current is the sum of 1/N given to each neighbour cell
        inward_current = np.zeros(self.size*self.size)
        for i in range(len(neighbours_list)):
            inward_current[np.array(neighbours_list[i],dtype = int)] += float(1)/np.array(resting_neighbours[i])
            
        # resting neighbours get excited if they are over the threshold
        get_excited = np.where(inward_current >= self.threshold)[0]
        
        self.tbe[get_excited] = True
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
            self.excitations[self.states[0]] += 1
            
            # is in SR
            if max(self.excitations[self.first_col]) == max(self.excitations):
                if self.t_AF != 0 :
                    self.tot_AF += self.t_AF
                    self.t_AF = 0
                    
            # is in AF        
            if max(self.excitations[self.first_col]) < max(self.excitations):
                self.t_AF += 1
                     
            self.t += 1
        self.tot_AF += self.t_AF
        print(self.tot_AF)
 
    
    def CMP2D_timestep_perc(self):
        """A single timestep"""
        self.SinusRhythm()
        self.Relaxing()
        self.Conduct()
        while len(self.states[0]) != 0:
            self.Relaxing()
            self.Conduct()
            self.t += 1
            for i in self.states[0]:
                if i % self.size == self.size-1:
                    return 1

        return 0