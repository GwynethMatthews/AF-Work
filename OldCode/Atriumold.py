import numpy as np
import scipy.stats as stats
""" Standard Conduction Model"""


class Atrium():
    """Creates the myocardium's structure.
    Model 1 = standard CMP model, 
    Model 2 = source model,
    Model 3 = source-sink model
    """
    def __init__(self, hexagonal = False, model = 1, L = 200, v_para = 1,
                 v_tran_1 = 0.1, v_tran_2 = 0.6, d = 0.05, threshold_cells = 1,
                 threshold = 0.25, e = 0.05, rp = 50, tot_time = 10**6, 
                 pace_rate = 220, s1 = 1, s2 = 2, s3 = 3, s4 = 4):
        
        # System Parameters
        self.hexagonal = hexagonal
        self.model = model
        self.L = L 
        
        if self.hexagonal == False:
            self.transverse_prob = v_tran_1
            
        if self.hexagonal == True:
            self.transverse_prob_l = v_tran_1
            self.transverse_prob_r = v_tran_2
            
        self.parallel_prob = v_para    
        self.tot_time = tot_time
        self.rp = rp        
        self.pace_rate = pace_rate
        self.pace = np.arange(0, tot_time, pace_rate)
        
        if self.model == 1:
            self.dysfunctional_prob = d
            self.nonfire_prob = e
            self.dysfunctional_cells = np.full([L*L],fill_value = False, dtype = bool) 
            
        if self.model == 2:
            self.nonfire_prob = e
            self.thresh_cells = threshold_cells
        
        if self.model == 3:
            self.threshold = threshold
            
        
        # System cell positions
        self.first_col = np.arange(0, L * L, L)
        self.last_col = np.arange(0, L * L, L) + L - 1
        
        self.index = np.arange(0, L*L) # cell positions in each array
        self.position = self.index.reshape(L, L)
        
        self.y = np.indices((L, L))[0] # y coordinate for cells
        self.x = np.indices((L, L))[1] # x coordinate for cells
        
        
        #System seeds
        self.seed_dysfunc = s1
        self.seed_connect_tran = s2
        self.seed_connect_para = s3
        self.seed_prop = s4

        # Measurable Variables
        self.excitations = np.zeros(L*L)
        self.AF = False
        self.sources  = []
        self.t = 0
        self.tot_AF = 0 # overall
        self.t_AF = 0 # in this episode
        self.t_SR = 0 # in this period of SR

        # Neighbours
        a = np.full((L*L),fill_value = None,dtype = float)
        if self.hexagonal == False: # up, right, down, left 
            self.neighbours = np.array([a]*4)
            
        if self.hexagonal == True: # up_left, up_right, right, down_right, down_left, left
            self.neighbours = np.array([a]*6)      

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
            np.random.seed(self.seed_dysfunc)
            num_rand_dysfunc = np.random.rand(L*L)
            
            np.random.seed(self.seed_connect_tran)
            num_rand_tran = np.random.rand(L*L)
            
            np.random.seed(self.seed_connect_para)
            num_rand_para = np.random.rand(L*L)
   
            if self.model == 1:
                for j in self.index:
                    
                    if self.dysfunctional_prob > num_rand_dysfunc[j]: # dysfunctional
                        self.dysfunctional_cells[j] = False
                        
                    if self.dysfunctional_prob <= num_rand_dysfunc[j]: # functional
                        self.dysfunctional_cells[j] = True

            for j in self.index:
                
                if num_rand_para[j] <= self.parallel_prob:
                    
                    if j in np.arange(0, L * L, L):
                        self.neighbours[1][j] = int(j+1)
                        self.neighbours[3][j + 1] = int(j)
                    

                    elif j in (np.arange(0, L * L, L) + L - 1):
                        self.neighbours[1][j] = None
                        
                        
                    else:
                        self.neighbours[1][j] = int(j+1)             
                        self.neighbours[3][j + 1] = int(j)
                        
            for j in self.index:
                
                if num_rand_tran[j] <= self.transverse_prob: 
                    
                    if j in np.arange((L * L) - L, L * L):
                        self.neighbours[2][j] = j - ((L * L) - L)
                        self.neighbours[0][j - ((L * L) - L)] = j
                    
                    
                    else:
                        self.neighbours[2][j] = j + L
                        self.neighbours[0][j + L] = j
            
                
            self.list_of_neighbours =  [[self.neighbours[0][i],
                 self.neighbours[1][i],
                 self.neighbours[2][i],
                 self.neighbours[3][i]] for i in self.index] 
                   
        
        if self.hexagonal == True:
            np.random.seed(self.seed_dysfunc)
            num_rand_dysfunc = np.random.rand(L*L)
            
            np.random.seed(self.seed_connect_tran)
            num_rand_tran1 = np.random.rand(L*L)
            num_rand_tran2 = np.random.rand(L*L)
            
            
            np.random.seed(self.seed_connect_para)
            num_rand_para = np.random.rand(L*L)
            
            if self.model == 1:
                
                for j in self.index:
                
                    if d > num_rand_dysfunc[j]: # dysfunctional
                        self.dysfunctional_cells[j] = False
                        
                    if d <= num_rand_dysfunc[j]: # functional
                        self.dysfunctional_cells[j] = True
                        
            for j in self.index:
                
                if num_rand_para[j] <= self.parallel_prob:
                    
                    if j in self.last_col:
                        self.neighbours[2][j] = None
                        
                    else:
                        self.neighbours[2][j] = int(j + 1)
                        self.neighbours[5][j + 1] = int(j)
            
            for j in self.index:
                #even
                if j in self.position[np.arange(0, L, 2)]:
                    
                    if num_rand_tran1[j] <= self.transverse_prob_l:
                        
                        if j not in self.first_col:
                            self.neighbours[4][j] = j + L - 1
                            self.neighbours[1][j + L - 1] = j
                            
                    if num_rand_tran2 [j] <= self.transverse_prob_r:
                        self.neighbours[3][j] = j + L
                        self.neighbours[0][j + L] = j
                        
                #odd
                if j in self.position[np.arange(1, L, 2)]:
                    
                    if j in np.arange(L * L - L, L * L):
                        
                        if num_rand_tran1[j] <= self.transverse_prob_l:
                            self.neighbours[4][j] = j - ((L * L) - L)
                            self.neighbours[1][j - ((L * L) - L)] = j

                    else:
                        
                        if num_rand_tran1[j] <= self.transverse_prob_l:
                            self.neighbours[4][j] = j+L
                            self.neighbours[1][j + L] = j
                            
                        if num_rand_tran2 [j] <= self.transverse_prob_r:
                            if j not in self.last_col:
                                self.neighbours[3][j] = j + L + 1
                                self.neighbours[0][j + L + 1] = j
                                
                
            self.list_of_neighbours =  [[self.neighbours[0][i],
                 self.neighbours[1][i],
                 self.neighbours[2][i],
                 self.neighbours[3][i],
                 self.neighbours[4][i],
                 self.neighbours[5][i]] for i in self.index]  
        
            
        self.neighbour_list = np.array([[x for x in 
                                 self.list_of_neighbours[i] if str(x) 
                                 != 'nan'] for i in self.index])
                    
        
        if self.model == 1:
            
            self.first_dys = np.array(self.first_col
                                      [~self.dysfunctional_cells[self.first_col]])
            self.first_fun = np.array(self.first_col
                                      [self.dysfunctional_cells[self.first_col]])
        
        
    def SinusRhythm(self):
        
        if self.model == 1:
        
            e_comp_val1 = np.random.rand(len(self.first_dys))
            dysfunctional_cells = self.first_dys[e_comp_val1 > self.nonfire_prob]
    
            self.tbe[dysfunctional_cells] = True
            self.tbe[self.first_fun] = True
            
        else:
            self.tbe[self.first_col] = True
        
        
    def Relaxing(self):

        self.resting[self.tbe] = False
        self.resting[self.states[-1]] = True
        
        del self.states[-1]
        self.states.insert(0, self.index[self.tbe])
        
    def Relaxing_ani(self):

        self.phases[self.tbe] = 0
        self.phases[~self.resting] += 1
        
        self.V[self.tbe] = 20.
        self.V[~self.resting] -= 2.2
        
        self.resting[self.tbe] = False
        self.resting[self.states[-1]] = True
        
        del self.states[-1]
        self.states.insert(0, self.index[self.tbe])   
        
    def Conduct(self):
        
        x = self.states[0]
        
        if self.model == 3:
            
            neighbours_list = [[y for y in self.neighbour_list[i] if self.resting[int(y)] == True] for i in x]

            resting_neighbours = list(map(len,neighbours_list))
            inward_current = np.zeros(self.L * self.L)
            
            for i in range(len(neighbours_list)):
                
                if resting_neighbours[i] != 0:
                    inward_current[np.array(neighbours_list[i],dtype = int)] += float(1) / np.array(resting_neighbours[i])
            
            #get_excited = np.where(inward_current >= self.threshold)[0]
            
            threshold = np.random.rand(self.L * self.L)* (0.5) 
            get_excited = self.index[inward_current >= threshold]
            
            self.tbe[get_excited] = True
        
        else:
            
            if self.hexagonal == False:
                
                neighbours = np.array([self.neighbours[0][x],
                                       self.neighbours[1][x],
                                       self.neighbours[2][x],
                                       self.neighbours[3][x]])
    
            if self.hexagonal == True:
                
                neighbours = np.array([self.neighbours[0][x],
                                       self.neighbours[1][x],
                                       self.neighbours[2][x],
                                       self.neighbours[3][x],
                                       self.neighbours[4][x],
                                       self.neighbours[5][x]])

            
            neighbours = np.array(neighbours[~np.isnan(neighbours)], dtype = int) 
            neighbours = neighbours[self.resting[neighbours]]
            
            if self.model == 1:
                
                neighbours_fun = neighbours[self.dysfunctional_cells[neighbours]]
                neighbours_dys = neighbours[~self.dysfunctional_cells[neighbours]]
            
            if self.model == 2:
                
                excitees, counts = np.unique(neighbours,return_counts = True)
                neighbours_fun = excitees[counts >= self.thresh_cells]
                neighbours_dys = excitees[counts < self.thresh_cells]
                
            e_comp_val2 = np.random.rand(len(neighbours_dys))
            neighbours_dys = neighbours_dys[e_comp_val2 > self.nonfire_prob]
    
            self.tbe[neighbours_fun] = True
            self.tbe[neighbours_dys] = True
        
        self.tbe[x] = False

    def TimeInAF(self):
        self.excitations[self.states[0]] += 1
        #if self.t in [396,397,398,399,400,401]:
        #    print(self.t)
        #    print(self.phases[23280])
        #    print(self.resting[23280])
        pacemaker_mode = int(np.average(self.excitations[self.first_col]) + 0.5)
       #print(pacemaker_mode)
       # print(self.excitations[9004])
        if self.AF == False:

            possible_AF_locs = np.where(self.excitations > pacemaker_mode)[0]
        
            neighbours_list = [[y for y in self.neighbour_list[i] if self.tbe[int(y)] == True] for i in possible_AF_locs]
            tbe_neighbours = np.array(list(map(len,neighbours_list)))
   
            if len(np.where(tbe_neighbours > 1)[0]) > 0:
  
                self.AF = True
                self.sources.extend(possible_AF_locs[np.where(tbe_neighbours)])
                #print(self.sources)

                #print(self.t)
                
        if self.AF == True:
            if max(self.excitations[np.array(self.sources)]) <= pacemaker_mode:
                self.AF = False
                #print(self.AF)
                #print(self.t)
        
    def CMP2D_timestep(self):

        if np.remainder(self.t, self.pace_rate) == 0:
            self.SinusRhythm()
            
        self.Relaxing()
        self.Conduct()
        self.TimeInAF()
        self.t += 1
        
    def CMP2D_timestep_ani(self):

        if np.remainder(self.t, self.pace_rate) == 0:
            self.SinusRhythm()
            
        self.Relaxing_ani()
        self.Conduct()
        self.TimeInAF()
        self.t += 1
            
    def CMP2D(self):

        np.random.seed(self.seed_prop)
        
        while self.t < self.tot_time:
            
            self.CMP2D_timestep()
           
#A = Atrium()