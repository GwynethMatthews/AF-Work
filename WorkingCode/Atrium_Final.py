import numpy as np

class Atrium:
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
    boundary: True for pipe, False fro sq slab
    pacemaker_line: True for first column, False for quarter circle corner
    radius: size of square SA node if pacemaker_line == False
    charge_conservation: if True then an excited cell with inward_current > threshold gives inward_current/N to each resting neighbour, if False gives 1/N
    t_under_on: whether or not the inward_current is kept if the cell doesn't reach the threshold
    t_under: the time the current is kept for
    """

    def __init__(self, hexagonal=False, Lx=200, Ly=200, rp=50, tot_time=10**6, nu_para=0.6, nu_trans=0.6,
                 pace_rate=220, p_nonfire=0.25, seed_connections=1, seed_prop=4, boundary=False,
                 pacemaker_line = True, radius = 3, charge_conservation = True, t_under_on = False, t_under = 3):

        global inward_current
        # System Parameters
        self.hexagonal = hexagonal
        self.boundary = boundary
        self.pacemaker_line = pacemaker_line
        self.radius = radius
        
        self.Lx = Lx
        self.Ly = Ly
        
        self.charge_conservation = charge_conservation

        self.t_under_on = t_under_on
        self.t_under = t_under

        self.nu_para = nu_para
        self.nu_trans = nu_trans
            
        self.tot_time = tot_time
        self.rp = rp        
        self.pace_rate = pace_rate
        self.pace = np.arange(0, tot_time, pace_rate)
        
        self.p_nonfire = p_nonfire

        # System cell positions
        self.index = np.arange(0, Lx * Ly)   # cell positions in each array
        self.position = self.index.reshape(Ly, Lx)     # Reshape works as (num rows, num columns)
        self.inward_current = np.zeros_like(self.index)

        self.current = np.ones_like(self.index, dtype = float)


        self.first_col = np.arange(0, Lx * Ly, Lx)
        self.not_first_col = self.index[self.index % Lx != 0]
        self.last_col = np.arange(0, Lx * Ly, Lx) + Lx - 1
        self.last_row = np.arange((self.Lx * self.Ly) - self.Lx, self.Lx * self.Ly)
        
        # System seeds
        self.seed_connections = seed_connections
        self.seed_prop = seed_prop

        # Measurable Variables for AF time
        # All dummy variables that get overridden in function
        self.excitation_rate = None
        self.last_excitation = None
        self.number_of_excitations = None
        self.inward_current = np.zeros_like(self.index, dtype = float)
        self.zero_current = np.full([self.Lx*self.Ly], fill_value=True, dtype=bool) 
        
        self.AF = None
        self.sources = None
        self.t = None
        self.tot_AF = None  # overall
        self.t_AF = None  # in this episode
        self.t_SR = None  # in this period of SR
        self.set_AF_measuring_vars()

        # State of system
        self.phases = np.full((Lx * Ly), fill_value=self.rp)         # state cell is in (0 = excited, rp = resting)
        self.V = np.full((Lx * Ly), fill_value=-90.0)                # voltage depending on state of cell
        
        # Added by Gwyn for Trial and Error data collection
        self.resting_cells = np.zeros((501))
        self.time_for_graphs = np.arange(-500,1,1)
        self.resting_cells_over_time = []
        self.receive_current = 0
        self.fail_safe = False
        self.time_extinguished = 0
        self.stop = False
        self.propagated = False
        self.propagation_time = 0
        self.time_increase_rp = self.pace_rate + self.rp
        
        self.states = [[]] * self.rp              # list of lists containing cells in each state except resting

        self.resting = np.full([self.Lx * Ly], fill_value=True, dtype=bool)         # can they be excited
        self.to_be_excited = np.full([self.Lx * Ly], fill_value=False, dtype=bool)        # cells to be excited next timestep

        self.t_under_states = [[]] * self.t_under


        self.neighbours = None    # Dummy that gets overwritten in create_neighbours
        self.list_of_neighbours = None  # Dummy that gets overwritten in create_neighbours
        
        self.create_neighbours()

        self.neighbour_list = np.array([np.array([x for x in 
                                 self.list_of_neighbours[i] if str(x) 
                                 != 'nan'],dtype = int) for i in self.index])
        self.pacemaker_cells = None
        self.create_pacemaker_cells()

    def set_AF_measuring_vars(self):

        self.excitation_rate = np.zeros(self.Lx * self.Ly, dtype=int)
        self.last_excitation = np.full((self.Lx * self.Ly), fill_value=-self.pace_rate)
        self.number_of_excitations = np.zeros(self.Lx * self.Ly, dtype=int)


        self.AF = False
        self.sources = []
        self.t = 0
        self.tot_AF = 0  # overall
        self.t_AF = 0  # in this episode
        self.t_SR = 0  # in this period of SR

    def create_neighbours(self):
        # Setting connections
        a = np.full((self.Lx * self.Ly), fill_value=None, dtype=float)
        
        if self.hexagonal:  # if self.hexagonal == True:
            neighbours = np.array([a] * 6)  # up_left, up_right, right, down_right, down_left, left

            np.random.seed(self.seed_connections)
            rand_nums = np.random.rand(3, self.Lx * self.Ly)

            num_rand_tran1 = rand_nums[0]
            num_rand_tran2 = rand_nums[1]
            num_rand_para = rand_nums[2]

            for j in self.index:

                if num_rand_para[j] <= self.nu_para:

                    if j not in self.last_col:
                        neighbours[2][j] = int(j + 1)
                        neighbours[5][j + 1] = int(j)

                # even
                if j in self.position[np.arange(0, self.Ly, 2)]:

                    if num_rand_tran1[j] <= self.nu_trans:
                        if j not in self.first_col:
                            neighbours[4][j] = j + self.Lx - 1
                            neighbours[1][j + self.Lx - 1] = j

                    if num_rand_tran2[j] <= self.nu_trans:
                        neighbours[3][j] = j + self.Lx
                        neighbours[0][j + self.Lx] = j

                # odd
                else:
                #if j in self.position[np.arange(1, Ly, 2)]:
                    if num_rand_tran1[j] <= self.nu_trans:
                        
                        if j in self.last_row:
                            if self.boundary == True: ### last row
                                neighbours[4][j] = j - ((self.Lx * self.Ly) - self.Lx)
                                neighbours[1][j - ((self.Lx * self.Ly) - self.Lx)] = j
                            
                        else:
                            neighbours[4][j] = j + self.Lx
                            neighbours[1][j + self.Lx] = j
                            

                    if num_rand_tran2[j] <= self.nu_trans:
                        if j not in self.last_col:
                            if j in self.last_row:
                                if self.boundary == True:
                                    neighbours[3][j] = j - ((self.Lx * self.Ly) - self.Lx) + 1
                                    neighbours[0][j - ((self.Lx * self.Ly) - self.Lx) + 1] = j
                              
                            else:    
                                neighbours[3][j] = j + self.Lx + 1
                                neighbours[0][j + self.Lx + 1] = j
                            
            self.list_of_neighbours = [[neighbours[0][i],
                                        neighbours[1][i],
                                        neighbours[2][i],
                                        neighbours[3][i],
                                        neighbours[4][i],
                                        neighbours[5][i]] for i in self.index]
            
        else:    # Square lattice
            neighbours = np.array([a] * 4)

            np.random.seed(self.seed_connections)
            num_rand_tran = np.random.rand(self.Lx * self.Ly)
            num_rand_para = np.random.rand(self.Lx * self.Ly)

            for j in self.index:

                if num_rand_para[j] <= self.nu_para:

                    if j not in self.last_col:
                        neighbours[1][j] = int(j + 1)
                        neighbours[3][j + 1] = int(j)


                if num_rand_tran[j] <= self.nu_trans:

                    if j in self.last_row:
                        if self.boundary == True:
                            neighbours[2][j] = j - ((self.Lx * self.Ly) - self.Lx)
                            neighbours[0][j - ((self.Lx * self.Ly) - self.Lx)] = j

                    else:
                        neighbours[2][j] = j + self.Lx
                        neighbours[0][j + self.Lx] = j

            self.list_of_neighbours = [[neighbours[0][i],
                                        neighbours[1][i],
                                        neighbours[2][i],
                                        neighbours[3][i]] for i in self.index]

        self.neighbours = neighbours

        self.neighbour_list = np.array([np.array([x for x in 
                                 self.list_of_neighbours[i] if str(x) 
                                 != 'nan'],dtype = int) for i in self.index])
            
    def create_pacemaker_cells(self):
        if self.pacemaker_line == True:
            self.pacemaker_cells = self.first_col
            
        else:
            self.pacemaker_cells = np.concatenate([np.arange(self.radius) + (self.Lx*i) for i in range(self.radius)])    ### Might be Ly???

    def change_resting_cells(self):
        self.resting[self.to_be_excited] = False    # Sets recently excited cells to False (not resting)
        self.resting[self.states[-1]] = True     # Sets previously refractory cells to resting

    def cycle_through_states(self):
        del self.states[-1]    # Pushes last cells out of states (not resting)
        self.states.insert(0, self.index[self.to_be_excited])     # Inserts new excited cells

    def relaxing(self):
        self.change_resting_cells()
        self.cycle_through_states()

    def excitation_tracker(self, excited_cells):
        self.excitation_rate[excited_cells] = self.t - self.last_excitation[excited_cells]
        self.last_excitation[excited_cells] = self.t
        self.number_of_excitations[excited_cells] += 1 

    def AF_checker(self, excited_cells):
        # print('here')
        #b = self.excitation_rate[excited_cells[self.number_of_excitations[excited_cells] > 1]]
        b = self.excitation_rate[excited_cells]
        if len(b) > 0:
            a = np.mean(b)


            if a < self.pace_rate - 1: # and beat <= self.pace_rate * 1.5:
                self.AF = True
                #print(self.AF)
                self.t_AF += 1
                #print(self.t_AF)
                ##print(self.AF)
                #print(self.t)
    
            else:
                self.AF = False

    def excite_cells(self, excited_cells):
        self.to_be_excited[excited_cells] = True

    def unexcite_cells(self, excited_cells):
        self.to_be_excited[excited_cells] = False

    def time_in_AF(self, excited_cells):
        if excited_cells.size > 0:      # Checks if non empty
            
            self.excitation_tracker(excited_cells)    # Tracks the excitation rates and last excitation of cells
            #self.AF_checker(excited_cells)          # Checks and updates the AF time

    def sinus_rhythm(self):
        return None   # Dummy that gets overriden by inheriting classes

    def conduct(self):
        return None   # Dummy that gets overriden by inheriting classes

    def cmp_no_sinus(self):
        self.relaxing()    # Changes self.states, cycles cells through the refractory states

        excited_cells = self.states[0]
        if len(excited_cells) > 0:    # NOTE: added this so that the new method of finding inwards_current doesn't break the animation
            self.conduct(excited_cells)

        self.unexcite_cells(excited_cells)  # Were excited, now refractory

        self.time_in_AF(excited_cells)
        self.t += 1

        
    def cmp_timestep(self):
        self.sinus_rhythm() # Now checks if correct time inside this function
        self.cmp_no_sinus()


    def cmp_animation(self):
        self.phases[~self.resting] += 1        
        self.phases[self.to_be_excited] = 0  # Needed for animation

        self.cmp_no_sinus()

    def cmp_full(self):
        np.random.seed(self.seed_prop)   # Sets seed for all dysfunctional firings etc.
   
        while self.t < self.tot_time:
            self.cmp_timestep()
            
        self.tot_AF += self.t_AF
    
    def change_connections(self, new_nu_para, new_nu_trans):
        self.nu_trans = new_nu_trans
        self.nu_para = new_nu_para
        self.create_neighbours()
    
    def change_rp(self, increment):
        new_rp = self.rp + increment
        
        if increment >= 0:
            self.states.extend([[]]*(new_rp-self.rp))
            
            self.phases[self.phases == self.rp] = new_rp
            self.phases[self.resting] = new_rp
            
        elif increment < 0:
            self.resting[np.concatenate(self.states[increment:])] = True
            
            self.phases[np.concatenate(self.states[increment:])] = new_rp
            self.phases[self.phases == self.rp] = new_rp
            
            del self.states[increment:]
            
        self.rp = new_rp

    def resting_cells_over_time_collect(self):
        self.resting_cells_over_time.extend([len(self.resting[self.resting == True])])
        
    def find_propagation_time(self):
        if self.propagated == False:
            if sum(self.number_of_excitations[self.last_col]) > 0:
                self.propagated = True
                self.propagation_time = self.t
                
    def pacing_with_change_of_rp_then_no_sinus(self, time_between_pace_and_change_of_rp, increment, num_paces):
        """ time_between_pace_and_change_of_rp is the time between the pace and t_c where all cells excited at t > t_c
        will have the new rp (e.g. if = 0 then all cells that excite after a new pace will have the new refractory period)
        increment is the change in rp (if set to 0 then rp doesn't change, normal pacing)"""

        np.random.seed(self.seed_prop)   # Sets seed for all dysfunctional firings etc.

        while self.t < (num_paces * self.pace_rate):
            self.sinus_rhythm() # self.sinus_rhythm checks whether t % pace_rate == 0 
            
            if self.t == (self.time_increase_rp + time_between_pace_and_change_of_rp):
                # changes the time for the next increase in rp
                self.time_increase_rp += increment + self.pace_rate
                self.change_rp(increment)  
                
            self.cmp_no_sinus()
            
        while self.t < self.tot_time:
            self.cmp_no_sinus()
            
        self.tot_AF += self.t_AF
        
    def pacing_with_change_of_rp_ani(self, time_between_pace_and_change_of_rp, increment):
#         """ time_between_pace_and_change_of_rp is the time between the pace and t_c where all cells excited at t > t_c
#        will have the new rp (e.g. if = 0 then all cells that excite after a new pace will have the new refractory period)
#        increment is the change in rp (if set to 0 then rp doesn't change, normal pacing)"""
        
        self.sinus_rhythm() # self.sinus_rhythm checks whether t % pace_rate == 0 
        
        if self.t == (self.time_increase_rp + time_between_pace_and_change_of_rp):
            # changes the time for the next increase in rp
            self.time_increase_rp += increment + self.pace_rate
            self.change_rp(increment)  

        self.cmp_animation()    # Doesn't have a sinus rhythm
        #self.cmp_no_sinus()
            

class DysfuncModel(Atrium):
    
    def __init__(self, seed_dysfunc=1, dysfunctional_prob=0.05, hexagonal=False, rp=50, tot_time=10**6, nu_para=0.6, nu_trans=0.6,
                 pace_rate=220, p_nonfire=0.05, seed_connections=1, seed_prop=4, boundary=True, pacemaker_line=True, radius = 3):
        super(DysfuncModel, self).__init__(hexagonal, rp, tot_time, nu_para, nu_trans, pace_rate, p_nonfire, seed_connections, seed_prop, boundary, pacemaker_line, radius)     # Calls Atrium init function
        
        self.seed_dysfunc = seed_dysfunc
        
        self.dysfunctional_prob = dysfunctional_prob
        self.dysfunctional_cells = np.full([self.Lx * self.Ly], fill_value=False, dtype=bool)

        self.set_dysfunctional_cells()

        functional_first_col_positions = self.dysfunctional_cells[self.first_col]

        self.first_dys = np.array(self.first_col[~functional_first_col_positions])
        self.first_fun = np.array(self.first_col[functional_first_col_positions])

    def set_dysfunctional_cells(self):
        np.random.seed(self.seed_dysfunc)
        num_rand_dysfunc = np.random.rand(self.Lx * self.Ly)

        for j in self.index:
            if self.dysfunctional_prob <= num_rand_dysfunc[j]:  # functional
                self.dysfunctional_cells[j] = True

    def sinus_rhythm(self):
        if self.t % self.pace_rate == 0:
            dysfunc_fire_rand_nums = np.random.rand(len(self.first_dys))
            successful_dysfunctional_cells = self.first_dys[dysfunc_fire_rand_nums > self.p_nonfire]

            self.excite_cells(successful_dysfunctional_cells)
            self.excite_cells(self.first_fun)

    def resting_neighbours(self, excited_cells):
        if self.hexagonal:
            neighbours = np.array([self.neighbours[0][excited_cells], self.neighbours[1][excited_cells],
                                   self.neighbours[2][excited_cells], self.neighbours[3][excited_cells],
                                   self.neighbours[4][excited_cells], self.neighbours[5][excited_cells]])

        if not self.hexagonal:  # square lattice, rarely used so ignore check
            neighbours = np.array([self.neighbours[0][excited_cells], self.neighbours[1][excited_cells],
                                   self.neighbours[2][excited_cells], self.neighbours[3][excited_cells]])

        neighbours = np.array(neighbours[~np.isnan(neighbours)], dtype=int)
        neighbours = neighbours[self.resting[neighbours]]

        return neighbours

    def get_dysfunctional_neighbour_indices(self, neighbours):

        return self.dysfunctional_cells[neighbours]

    def conduct(self, excited_cells):
        neighbours = self.resting_neighbours(excited_cells)

        dysfunc_neighbour_indices = self.get_dysfunctional_neighbour_indices(neighbours)

        func_neighbours = neighbours[dysfunc_neighbour_indices]
        dysfunc_neighbours = neighbours[~dysfunc_neighbour_indices]

        dysfunc_fire_rand_nums = np.random.rand(len(dysfunc_neighbours))
        dysfunc_neighbours = dysfunc_neighbours[dysfunc_fire_rand_nums > self.p_nonfire]

        self.excite_cells(func_neighbours)
        self.excite_cells(dysfunc_neighbours)


class SourceSinkModel(Atrium):
    
    def __init__(self, threshold=1, hexagonal=True, Lx=70, Ly=100, rp=70, tot_time=5*10**4, nu_para=1, nu_trans=1,
                 pace_rate=220, p_nonfire=0.75, seed_connections=1, seed_prop=4, boundary=True, pacemaker_line=True, radius=3, charge_conservation = False, t_under_on = False, t_under = 3):

        super(SourceSinkModel, self).__init__(hexagonal, Lx, Ly, rp, tot_time, nu_para, nu_trans, pace_rate, p_nonfire, seed_connections, seed_prop, boundary, pacemaker_line, radius, charge_conservation, t_under_on, t_under)       # Calls Atrium init function

        self.threshold = threshold

    def sinus_rhythm(self):
        if self.t % self.pace_rate == 0:
            self.to_be_excited[self.pacemaker_cells] = True
    
    def ectopic_beat(self, location_of_cells):
        self.to_be_excited[location_of_cells] = True


    def get_inward_current(self, neighbours_list,resting_neighbours,excited_cells):

        if self.charge_conservation == True:
            j = 0
            
            for i in neighbours_list:
                if len(i) != 0:
                    self.inward_current[i] += float(self.current[excited_cells[j]])/ len(i)

                j += 1
        else:
            for i in neighbours_list:
                if len(i) != 0:
                    self.inward_current[i] += float(1)/ len(i)
                    
        #return inward_current

    def cycle_through_t_under_states(self,keep_current):
        self.zero_current[self.t_under_states[-1]] = True
        self.zero_current[self.states[0]] = True
        
        del self.t_under_states[-1]
        
        self.t_under_states.insert(0, keep_current) 
        
        self.zero_current[self.t_under_states[0]] = False 
        
    def cells_miss_threshold_p_constant(self, receive_current):
        """Returns the cells which are excited even though they receive less than the threshold. 
        Probability of excitation is constant"""
        possible_excited = receive_current[self.inward_current[receive_current] < self.threshold]

        miss_threshold_fire_rand_nums = np.random.rand(len(possible_excited))

        possible_excited = possible_excited[miss_threshold_fire_rand_nums > self.p_nonfire]  # Fire if over p_nonfire
        
        if self.t_under_on == True: 
            keep_current = possible_excited[miss_threshold_fire_rand_nums <= self.p_nonfire]
            
            self.cycle_through_t_under_states(keep_current)
        
        return possible_excited

    def cells_miss_threshold_as_a_function(self, receive_current):
        """Returns the cells which are excited even though they receive less than the threshold. 
        Probability of excitation is 1 - (p_nonfire to the power of the current received)"""
        possible_excited = receive_current[self.inward_current[receive_current] < self.threshold]

        miss_threshold_fire_rand_nums = np.random.rand(len(possible_excited))
        get_excited = possible_excited[miss_threshold_fire_rand_nums > self.p_nonfire**(self.inward_current[possible_excited])]
        
        if self.t_under_on == True: 
            keep_current = possible_excited[miss_threshold_fire_rand_nums <= self.p_nonfire**(self.inward_current[possible_excited])]
        
            self.cycle_through_t_under_states(keep_current)
        
        return get_excited            

    def find_resting_neighbours(self, excited_cells):
        neighbours_list = [j[self.resting[j]] for j in self.neighbour_list[excited_cells]]

        resting_neighbours = list(map(len, neighbours_list))

        return neighbours_list, resting_neighbours
    
    def reset_cells(self):
        if self.t_under_on == True:
            self.inward_current[self.zero_current] = 0
        
        if self.t_under_on == False:
            self.inward_current = np.zeros(self.Lx*self.Ly)
        


    def conduct(self, excited_cells):
        neighbours_list, resting_neighbours = self.find_resting_neighbours(excited_cells)
        
        self.reset_cells()
        
        self.get_inward_current(neighbours_list,resting_neighbours,excited_cells)  # amount of current received
        
        receive_current = self.index[self.inward_current > 0]  # Indices which receive any current from neighbours
        
        hit_thresh_so_excite = receive_current[self.inward_current[receive_current] >= self.threshold]
        miss_thresh_but_still_excite = self.cells_miss_threshold_as_a_function(receive_current)
        
        if self.charge_conservation == True:
            self.current = np.ones_like(self.index, dtype = float)
            self.current[hit_thresh_so_excite] = self.inward_current[hit_thresh_so_excite]
            #self.current[miss_thresh_but_still_excite] = self.inward_current[miss_thresh_but_still_excite]
        
        self.to_be_excited[miss_thresh_but_still_excite] = True
        self.to_be_excited[hit_thresh_so_excite] = True
        

