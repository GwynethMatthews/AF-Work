
import numpy as np

class Atrium():
    """Creates the myocardium's structure.
    hexagonal: False = Square Lattice, True = Hexagonal Lattice
    model: 1 = standard CMP model, 2 = source-sink model
    L: Lattice length
    v_para: nu parallel 
    v_tran_1: nu transfers in square lattice or nu of main diagonal in hexagonal
    v_tran_2: nu of minor diagonal in hexagonal
    threshold: value equal to and above which cells always excite in the source-sink model
    p: the probability with which cells under the threshold excite in the source-sink model
    rp: refractory period
    tot_time: total time of simulation
    pace_rate: rate of pacemaker rhythm
    s1: seed for transverse connection selection
    s2: seed for parallel connection selection
    s3: seed for cell firing selection (e or p)
    """
    def __init__(self, hexagonal = False, L = 200, v_para = 1,
                 v_tran_1 = 0.1, v_tran_2 = 0.6, 
                 threshold = 0.25, p = 0.05, rp = 50, tot_time = 10**6,
                 pace_rate = 220, s1 = 1, s2 = 2, s3 = 3):
        
        # System Parameters
        self.hexagonal = hexagonal
        self.L = L
        
        self.transverse_prob_l = v_tran_1
        self.transverse_prob_r = v_tran_2
        self.parallel_prob = v_para  
        
    
        self.tot_time = tot_time
        self.rp = rp        
        self.pace_rate = pace_rate
        self.pace = np.arange(0, tot_time, pace_rate)
        
        self.threshold = threshold
        self.p = p
            
        
        # System cell positions
        self.index = np.arange(0, L*L) # cell positions in each array
        self.position = self.index.reshape(L, L)

        self.first_col = np.arange(0, L * L, L)
        self.not_first_col = self.index[self.index % L != 0]
        self.last_col = np.arange(0, L * L, L) + L - 1
        
        #System seeds
        self.seed_connect_tran = s1
        self.seed_connect_para = s2
        self.seed_prop = s3

        # Measurable Variables
        self.excitation_rate = np.zeros(L*L, dtype = int)
        self.last_excitation = np.full((L*L), fill_value = -self.pace_rate)
        self.AF = False
        self.fail_safe = False
        self.stop = False
        self.time_extinguished = 0
        #self.sources  = []
        self.t = 0
        self.tot_AF = 0 # overall
        self.t_AF = 0 # in this episode
        #self.t_SR = 0 # in this period of SR

        # Neighbours
        a = np.full((L*L),fill_value = None,dtype = float)
        self.neighbours = np.array([a]*6)      

        # State of system
        
        self.states = [[]]*self.rp # list of lists containing cells in each state except resting
        self.resting = np.full([L*L],fill_value = True, dtype = bool) # can they be excited
        self.tbe = np.full([L*L],fill_value = False, dtype = bool) # cells to be excited in next timestep
        
        # Setting connections and dysfunctional cells
         

        np.random.seed(self.seed_connect_tran)
        num_rand_tran1 = np.random.rand(L*L)
        num_rand_tran2 = np.random.rand(L*L)
        
        
        np.random.seed(self.seed_connect_para)
        num_rand_para = np.random.rand(L*L)
        
                    
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
                    self.neighbours[0][j + self.L] = j
                    
            #odd
            else:
            #if j in self.position[np.arange(1, L, 2)]:
                
                if j in np.arange(self.L * self.L - self.L, self.L * self.L): ### last row
                    
                    if num_rand_tran1[j] <= self.transverse_prob_l:
                        self.neighbours[4][j] = j - ((self.L * self.L) - self.L)
                        self.neighbours[1][j - ((self.L * self.L) - self.L)] = j

                else:
                    
                    if num_rand_tran1[j] <= self.transverse_prob_l:
                        self.neighbours[4][j] = j+L
                        self.neighbours[1][j + self.L] = j
                        
                if num_rand_tran2 [j] <= self.transverse_prob_r:
                    
                    if j not in self.last_col:
                        
                        if j in np.arange(self.L * self.L - self.L, self.L * self.L):
                            self.neighbours[3][j] = j - ((self.L*self.L) - self.L) + 1
                            self.neighbours[0][j - ((self.L*self.L) - self.L) + 1] = j
                            
                        else:    
                            self.neighbours[3][j] = j + self.L + 1
                            self.neighbours[0][j + self.L + 1] = j
                            
            
        self.list_of_neighbours =  [[self.neighbours[0][i],
             self.neighbours[1][i],
             self.neighbours[2][i],
             self.neighbours[3][i],
             self.neighbours[4][i],
             self.neighbours[5][i]] for i in self.index]  
        
            
        self.neighbour_list = np.array([[x for x in 
                                 self.list_of_neighbours[i] if str(x) 
                                 != 'nan'] for i in self.index])
            
    def SinusRhythm(self):
        
            self.tbe[self.first_col] = True
           
    def Relaxing(self):

        self.resting[self.tbe] = False
        self.resting[self.states[-1]] = True
        
        del self.states[-1]
        self.states.insert(0, self.index[self.tbe])

    
    def Conduct(self):
        
        x = self.states[0]
            
        neighbours_list = [[y for y in self.neighbour_list[i] if self.resting[int(y)] == True] for i in x]

        resting_neighbours = list(map(len,neighbours_list))
        inward_current = np.zeros(self.L * self.L)
        
        for i in range(len(neighbours_list)):
            
            if resting_neighbours[i] != 0:
                inward_current[np.array(neighbours_list[i],dtype = int)] += float(1) / np.array(resting_neighbours[i])
        
        receive_current = self.index[inward_current > 0]
<<<<<<< HEAD
        print(inward_current)
=======
        
>>>>>>> 9ef767596c9c4cc3e3587084cba1fb3187197814
        get_excited = receive_current[inward_current[receive_current] >= self.threshold]
        
        possible_excited = receive_current[inward_current[receive_current] < self.threshold]
        e_comp_val3 = np.random.rand(len(possible_excited))
        possible_excited = possible_excited[e_comp_val3 <= self.p]

        self.tbe[possible_excited] = True
        self.tbe[get_excited] = True
        
        self.tbe[x] = False
        
    


    def TimeInAF(self):

        if len(self.states[0]) > 0: 
            
            x = self.states[0]
            
            self.excitation_rate[x] = np.abs(self.last_excitation[x] - self.t)
            self.last_excitation[x] = self.t
        
            a = np.mean(self.excitation_rate[x])

            if a < self.pace_rate * 0.9:
                self.AF = True
                self.tot_AF += 1
            
            else:
                self.AF = False



# parameters array of [v,t,p] of length 60, i th ejob variable
def OnePacemakerBeat(parameters, seeds, itr):
    data_full = []
    
    for l in range(len(parameters[itr])):
        repeat_data = []
        print(l)
        for m in range(2):
            print(m)
            A = Atrium(hexagonal = True, L = 200, v_para = parameters[itr][l][0],
                         v_tran_1 = parameters[itr][l][0], v_tran_2 = parameters[itr][l][0],
                         threshold = parameters[itr][l][1], p = parameters[itr][l][2], rp = 50, tot_time = 10**6,
                         pace_rate = 220, s1 = seeds[itr][l][m][0], s2 = seeds[itr][l][m][1], s3 = seeds[itr][l][m][2])
    
            A.SinusRhythm()
            A.Relaxing()
            A.Conduct()
            while A.stop == False:
                
                if len(A.states[0]) != 0:
                    
                    if A.t != 5000:
                        A.Relaxing()
                        A.Conduct()
                        A.t += 1     
                        
                        A.TimeInAF()
                        
                        if A.AF == True:
                            A.t_AF = A.t
                            A.stop = True
                            #A.tot_AF += 1
                            
            
                    else:
                        
                        A.fail_safe = True
                        A.time_extinguished = A.t
                        A.stop = True
                    
                
                else: 
                    A.time_extinguished = A.t
                    A.stop = True
        # nu, threshold, p, s1, s2, s3, whether it extinguishes at 5000, whether it
        # eneters AF, when it eneters AF (0 if not AF), time wave is extinguished 
        # (5000 if fail safe, 0 if enteres AF)
        
            data = np.array([parameters[itr][l][0], parameters[itr][l][1], parameters[itr][l][2], 
                         A.seed_connect_tran, A.seed_connect_para, A.seed_prop, 
                         A.fail_safe, A.AF, A.t_AF, A.time_extinguished], dtype = float)
            repeat_data.extend([data])
            
        data_full.extend([repeat_data])
        
    data_full = np.array(data_full, dtype = float)
    np.save('test_data',data_full)
    
    
parameters = []
for i in np.linspace(0.01, 1, 10, endpoint = True): # nu values
    for j in np.array([0.25,0.3,0.4,0.45,0.5,0.75]): # threshold values
        for k in np.linspace(0, 0.99, 10, endpoint = True): # p values
            parameters.extend([[i,j,k]])
parameters = np.array(parameters).reshape((12,50,3))
s = np.random.randint(0,2**31,(1000,60,2,3))
