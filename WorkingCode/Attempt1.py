import numpy as np
#import scipy.stats as stats
""" Standard Conduction Model"""


class Atrium():
    """Creates the myocardium's structure.
    hexagonal: False = Square Lattice, True = Hexagonal Lattice
    model: 1 = standard CMP model, 2 = source-sink model
    L: Lattice length
    v_para: nu parallel 
    v_tran_1: nu transfers in square lattice or nu of main diagonal in hexagonal
    v_tran_2: nu of minor diagonal in hexagonal
    d: fraction of dysfunctional cells in standard model
    e: fraction of time that dysfunctional cells do not fire
    threshold: value equal to and above which cells always excite in the source-sink model
    p: the probability with which cells under the threshold excite in the source-sink model
    rp: refractory period
    tot_time: total time of simulation
    pace_rate: rate of pacemaker rhythm
    s1: seed for dysfunctional cell selection
    s2: seed for transverse connection selection
    s3: seed for parallel connection selection
    s4: seed for cell firing selection (e or p)
    """
    def __init__(self, hexagonal = False, model = 1, L = 200, v_para = 1,
                 v_tran_1 = 0.1, v_tran_2 = 0.6, d = 0.05, e = 0.05, 
                 threshold = 0.25, p = 0.05, rp = 50, tot_time = 10**6,
                 pace_rate = 220, s1 = 1, s2 = 2, s3 = 3, s4 = 4):
        
        # System Parameters
        self.hexagonal = hexagonal ## hexagonal lattice or square
        self.model = model ## standard of source sink
        self.L = L ## length of lattice
        
        if self.hexagonal == False:
            self.transverse_prob = v_tran_1 ## up/down
            
        if self.hexagonal == True:
            self.transverse_prob_r = v_tran_1 ##  down right to up left
            self.transverse_prob_l = v_tran_2 ##  down left to up right
            
        self.parallel_prob = v_para    ## left/right
        self.tot_time = tot_time ## total number of steps
        self.rp = rp        ## refractory period 
        self.pace_rate = pace_rate ## timesteps between pacemaker beats
        self.threshold = threshold ## above which cells always excite in source-sink model
        self.p = p ## probability of cells bellow threshold exciting
        
        self.index = np.arange(0, L*L) ## cell positions in each array
        self.position = self.index.reshape(L, L)
        
        self.first_col = np.arange(0, L * L, L) ## cell index of left most column
        self.not_first_col = self.index[self.index % L != 0] ## cell index for all cells except the left most column
        self.last_col = np.arange(0, L * L, L) + L - 1  ## cell index of right most column
        self.last_row = np.arange((L * L) - L, L * L) ## cell index of last row
        
        # System seeds
        self.seed_dysfunc = s1 ## seed for dyfunctional cells in standard model
        self.seed_connections = s2 ## seed for cell connections
        self.seed_prop = s4 ## seed for dynamical probabilities e/p
        
        # Neighbours        
        if self.hexagonal == False:
            ## up, right, down, left 
            self.pos_neighbours = np.array([np.full((L*L),
                    fill_value = None,dtype = float)]*4) ## all neighbours (i.e. if all nu == 1)
            self.neighbours = np.array([np.zeros(L*L,
                    dtype = bool)]*4) ## whether they actually have the neighbours
            
            
            np.random.seed(self.seed_connections) ## setting random seed for connections
        
            num_rand_tran = np.random.rand(L*L) ## if < nu then cell will have a downward connection
            num_rand_para = np.random.rand(L*L) ## if < nu then cell will have a connection to the right
   
            for j in self.index:
                
                if j in self.last_col:
                    ## last column cannot have a connection to the right
                    self.pos_neighbours[1][j] = None
                    self.neighbours[1][j] = False
                    
                else:
                    ## all other cells can have a connection to the right and 
                    ## the cell to their right will have a connection to the left
                    self.pos_neighbours[1][j] = int(j+1)             
                    self.pos_neighbours[3][j + 1] = int(j)
                    
                    ## setting connections for this lattice
                    if num_rand_para[j] <= self.parallel_prob:
                        self.neighbours[1][j] = True
                        self.neighbours[3][j + 1] = True
                        
                
                #if num_rand_tran[j] <= self.transverse_prob: 
                    
                if j in self.last_row:
                    ## last row connects to first row
                    self.pos_neighbours[2][j] = j - ((L * L) - L)
                    self.pos_neighbours[0][j - ((L * L) - L)] = j
                    
                    ## setting connections for this lattice
                    if num_rand_tran[j] <= self.transverse_prob:
                        self.neighbours[2][j] = True
                        self.neighbours[0][j - ((L * L) - L)] = True
                
                else:
                    ## all other cells connect to the cell below
                    self.pos_neighbours[2][j] = j + L
                    self.pos_neighbours[0][j + L] = j
                    
                    ## setting connections for this lattice
                    if num_rand_tran[j] <= self.transverse_prob:
                        self.neighbours[2][j] = True
                        self.neighbours[0][j + L] = True
                        
            self.neighbour_list = np.array([[self.pos_neighbours[0][x][self.neighbours[0][x]],
                                         self.pos_neighbours[1][x][self.neighbours[1][x]],
                                         self.pos_neighbours[2][x][self.neighbours[2][x]],
                                         self.pos_neighbours[3][x][self.neighbours[3][x]]] for x in self.index])
                        
        if self.hexagonal == True: # up_left, up_right, right, down_right, down_left, left
            self.pos_neighbours = np.array([np.full((L*L),
                    fill_value = None, dtype = float)]*6) ## all neighbours (i.e. if all nu == 1) 
            self.neighbours = np.array([np.zeros(L*L,
                    dtype = bool)]*6) ## whether they actually have the neighbours
    
    
            np.random.seed(self.seed_connections) ## setting seed for connections
            
            num_rand_tran1 = np.random.rand(L*L) ## if < nu then cell will have a downward right connection
            num_rand_tran2 = np.random.rand(L*L) ## if < nu then cell will have a downward left connection 
            num_rand_para = np.random.rand(L*L) ## if < nu then cell will have a connection to the right
            
            for j in self.index:
                    
                if j in self.last_col:
                    ## last column can't have connections to the right 
                    ## regardless of odd or even
                    self.pos_neighbours[2][j] = None
                    self.neighbours[2][j] = False
                        
                else:
                    ## all other cells can have a connection to the right and 
                    ## the cell to their right will have a connection to the left
                    self.pos_neighbours[2][j] = int(j + 1)
                    self.pos_neighbours[5][j + 1] = int(j)
                    
                    ## setting connections for this lattice
                    if num_rand_para[j] <= self.parallel_prob:
                        self.neighbours[2][j] = True
                        self.neighbours[5][j + 1] = True
                        
            
                #even
                if j in self.position[np.arange(0, L, 2)]:
                        
                    if j not in self.first_col:
                        ## sets down right connections for all cells except the first column 
                        self.pos_neighbours[4][j] = j + L - 1
                        self.pos_neighbours[1][j + L - 1] = j
                        
                        ## connection goes down and to the left of cell j 
                        if num_rand_tran2[j] <= self.transverse_prob_l:
                            self.neighbours[4][j] = True
                            self.neighbours[1][j + L - 1] = True
                    
                    ## all cells in the even rows can have down right connections
                    ## connection goes down and to the right of cell j 
                    self.pos_neighbours[3][j] = j + L
                    self.pos_neighbours[0][j + L] = j
                    
                    if num_rand_tran1[j] <= self.transverse_prob_r:
                        self.neighbours[3][j] = True
                        self.neighbours[0][j + L] = True
                #odd
                else:
                #if j in self.position[np.arange(1, L, 2)]:
                    
                    if j in self.last_row:
                        ## last row is odd and connects to top row
                        self.pos_neighbours[4][j] = j - ((L * L) - L)
                        self.pos_neighbours[1][j - ((L * L) - L)] = j
                        
                        if num_rand_tran2[j] <= self.transverse_prob_l:
                            self.neighbours[4][j] = True
                            self.neighbours[1][j - ((L * L) - L)] = True
    
                    else: ## not in last row
                        
                        self.pos_neighbours[4][j] = j + L
                        self.pos_neighbours[1][j + L] = j
                        
                        if num_rand_tran2[j] <= self.transverse_prob_l:
                            self.neighbours[4][j] = True
                            self.neighbours[1][j + L] = True
                            
                        
                    if j not in self.last_col:
                        ## last column cannot have a down right connection
                        
                        if j in self.last_row:
                            ## last row connects to first row
                            self.pos_neighbours[3][j] = j - ((L * L) - L) + 1
                            self.pos_neighbours[0][j - ((L * L) - L) + 1] = j
                            
                            if num_rand_tran1 [j] <= self.transverse_prob_r:
                                self.neighbours[3][j] = True
                                self.neighbours[0][j - ((L * L) - L) + 1] = True
                       
                        else:## not in last row
                            
                            self.pos_neighbours[3][j] = j + L + 1
                            self.pos_neighbours[0][j + L + 1] = j
                            
                            if num_rand_tran1 [j] <= self.transverse_prob_r:
                                
                                self.neighbours[3][j] = True
                                self.neighbours[0][j + L + 1] = True
                   
            self.neighbour_list = np.array([[self.pos_neighbours[0][x][self.neighbours[0][x]],
                                         self.pos_neighbours[1][x][self.neighbours[1][x]],
                                         self.pos_neighbours[2][x][self.neighbours[2][x]],
                                         self.pos_neighbours[3][x][self.neighbours[3][x]],
                                         self.pos_neighbours[4][x][self.neighbours[4][x]],
                                         self.pos_neighbours[5][x][self.neighbours[5][x]]] for x in self.index])
    
    
        self.phases = np.full((L*L),fill_value = self.rp) # state cell is in (0 = excited, rp = resting)       
        self.states = [[]]*self.rp # list of lists containing cells in each state except resting
        self.resting = np.full([L*L],fill_value = True, dtype = bool) # can they be excited
        self.tbe = np.full([L*L],fill_value = False, dtype = bool) # cells to be excited in next timestep
        
        self.t = 0
        
        np.random.seed(self.seed_prop)
    
        
    def SinusRhythm(self):
        
        self.tbe[self.first_col] = True
        
    def Relaxing(self):

        self.resting[self.tbe] = False
        self.resting[self.states[-1]] = True
        
        del self.states[-1]
        self.states.insert(0, self.index[self.tbe])
        
    def Relaxing_ani(self):

        self.phases[self.tbe] = 0
        self.phases[~self.resting] += 1
        
        # Needed for ECG
#        self.V[self.tbe] = 20.
#        self.V[~self.resting] -= 2.2
        
        self.resting[self.tbe] = False
        self.resting[self.states[-1]] = True
        
        del self.states[-1]
        self.states.insert(0, self.index[self.tbe])  
        
        
    def Conduct(self):
        
        x = self.states[0]
            
        neighbours_list = [i[self.resting[i]] for i in self.neighbour_list[x]]

        inward_current = np.zeros(self.L * self.L)
        
        for i in range(len(neighbours_list)):
            
            if len(neighbours_list[i]) != 0:
                inward_current[neighbours_list[i]] += float(1) / len(neighbours_list[i])
        
        receive_current = self.index[inward_current > 0]
       
        get_excited = receive_current[inward_current[receive_current] >= self.threshold]
        
        possible_excited = receive_current[inward_current[receive_current] < self.threshold]
        
        e_comp_val3 = np.random.rand(len(possible_excited))
        possible_excited = possible_excited[e_comp_val3 <= self.p]

        self.tbe[possible_excited] = True
        self.tbe[get_excited] = True
        
        self.tbe[x] = False
        
        
        
    def TimeInAF(self):

        x = self.states[0]
        #not_first_col = self.not_first_col
        if len(x) > 0: 
            
            #x = self.states[0]
            
            self.excitation_rate[x] = self.t - self.last_excitation[x]
            self.last_excitation[x] = self.t
        
            a = sum(self.excitation_rate[x])/len(self.excitation_rate[x])#np.mean(self.excitation_rate[x])

            #if self.t > self.pace_rate:
            if a < self.pace_rate * 0.9:
                self.AF = True
                self.t_AF += 1
                #print(self.AF)
                #print(self.t)
            
            else:
                self.AF = False
                #print(self.AF)
    
    def CMP2D_timestep2(self):

        if np.remainder(self.t, self.pace_rate) == 0:
            self.SinusRhythm()
            
        self.Relaxing()
        self.Conduct()
        self.TimeInAF()
        self.t += 1
    
    def CMP2D_timestep_ani2(self):

        if np.remainder(self.t, self.pace_rate) == 0:
            self.SinusRhythm()
            
        self.Relaxing_ani()
        self.Conduct()

        self.TimeInAF()
        self.t += 1
        
    def CMP2D(self):

        
        while self.t < self.tot_time:
            self.CMP2D_timestep2()
            
        