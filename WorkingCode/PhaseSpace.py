
import numpy as np
import sys
import math
import Atrium_Final as AF

#job_number = int(sys.argv[1])


#input_param = np.load('parameters.npy')
#input_seeds = np.load('seeds.npy')

#print(np.shape(input_param))

def OnePacemakerBeat(parameters, seeds, itr):
    data_full = []
    
    for l in range(len(parameters[itr])):   # New parameter set (nu, p, tau)
        repeat_data = []
        
        for i in range(500):  ### Number of repeats   

            nu = parameters[itr][l][0]

            tau = int(parameters[itr][l][1])
            p = parameters[itr][l][2]
            t_under_on = False # parameters[itr][l][3]
            pace = tau + 2
            
            avg_resting = None
            std_resting = None
            med_resting = None
            max_resting = None
            min_resting = None
            position_of_min = None
            position_of_max = None
            Fraction_of_resting_cells_last_timestep = None
            
            A = AF.SourceSinkModel(hexagonal=True, threshold=1, p_nonfire=p, pace_rate=pace,
                       Lx=70, Ly=100, tot_time=15000, nu_para=nu, nu_trans=nu, rp=tau,
                       seed_connections=seeds[itr][l][i][0], seed_prop=seeds[itr][l][i][1], 
                       charge_conservation = False, t_under = 1, t_under_on = t_under_on)

            AF_start = 40 * pace
            
            resting_cells_at_last_beat = 7000 # i.e. terminates before last beat (not possible)
            
            np.random.seed(A.seed_prop)

            A.cmp_timestep()   ### With one sinus beat       
            
            while A.stop == False:
                if A.t < A.tot_time: 
                    
                    if A.t < 31 * A.pace_rate:
                        # pacing
                        A.pacing_with_change_of_rp(time_between_pace_and_change_of_rp = 0,
                                 increment = -1)

                        A.find_propagation_time() # time till first cell in last column excites for the first time
                        
                        if A.propagated == True:
                            AF_start = int((31 * pace) + (5 * A.propagation_time))
                    
                    else:
                        if A.t == 31 * A.pace_rate: # fraction of resting cells at the last beat
                            resting_cells_at_last_beat = float(len(A.resting[A.resting == True])) / (A.Lx*A.Ly)
                        
                        if len(A.states[0]) != 0:   # continues to propagate
                            
                            A.cmp_no_sinus()
                            
                            if A.t > AF_start:
                                A.t_AF += 1
                                A.resting_cells_over_time_collect()

                        else:
                            # terminates
                            A.time_extinguished = A.t
                            A.stop = True
                            
                            Fraction_of_resting_cells_last_timestep = len(A.resting[A.resting == True]) / (A.Lx * A.Ly)
                            
                            
                else:
                    # reaches tot_time
                    A.fail_safe = True
                    A.time_extinguished = A.t
                    A.stop = True

                    Fraction_of_resting_cells_last_timestep = len(A.resting[A.resting == True]) / (A.Lx * A.Ly)
                    
            
            if A.t_AF > 0:
                A.AF = True
            
            if len(A.resting_cells_over_time) > 200:     ### Want to ignore last 200 values
            
                resting_cells_minus_slice = np.array(A.resting_cells_over_time[:-200])
            
                avg_resting = np.mean(resting_cells_minus_slice) / (A.Lx * A.Ly)    
                std_resting = np.std(resting_cells_minus_slice) / (A.Lx * A.Ly)
                med_resting = np.median(resting_cells_minus_slice) / (A.Lx * A.Ly)
                max_resting = max(resting_cells_minus_slice) / (A.Lx * A.Ly)
                min_resting = min(resting_cells_minus_slice) / (A.Lx * A.Ly)
                
                position_of_min = np.where(resting_cells_minus_slice == min(resting_cells_minus_slice))[0][0] + AF_start # time of min
                position_of_max = np.where(resting_cells_minus_slice == max(resting_cells_minus_slice))[0][0] + AF_start # time of max
                
                
            # nu
            # tau
            # p
            ### whether the charge under threshold is conserved
            ### pace_rate
            ### nu_para/nu_trans
            # repeat number
            # seed_connection
            # seed_propagation
            # A.fail_safe = whether it extinguishes at tot_time
            # A.AF = whether it eneters AF
            # A.t_AF = how long it was in AF for
            # A.time_extinguished = time wave is extinguished == tot_time if doesn't terminate
            # time of AF starting
            # average number of resting cells
            # std of number of resting cells
            # median of number of resting cells
            # min of number of resting cells
            # max of number of resting cells
            # fraction of resting cells at last beat
            # fraction of resting cells at last time step either self.states[0] == 0 or self.t == tot_time
            # position of min
            # position of max
            
            data = np.array([parameters[itr][l][0]*100, parameters[itr][l][1], parameters[itr][l][2]*100, 
                             i, A.seed_connections, A.seed_prop,
                             A.fail_safe, A.AF, A.t_AF, A.time_extinguished, AF_start, avg_resting, std_resting,
                             med_resting, min_resting, max_resting, resting_cells_at_last_beat,
                             Fraction_of_resting_cells_last_timestep, 
                             position_of_min, position_of_max])
    
            print(data)

    
            repeat_data.extend([data])

        data_full.extend([repeat_data])

    data_full = np.array(data_full)
    np.save('onesr_data_'+str(itr),data_full)


#OnePacemakerBeat(parameters=input_param, seeds=input_seeds, itr=job_number)

parameters = []
# nu, tau, p, whether the charge under threshold is conserved, amount to add to tau to get pace_rate 

for j in np.arange(50, 102, 2): # tau values
    for i in np.linspace(0.4, 1, 61, endpoint = True): # nu values
    #for i in [0.4,0.5,0.6,0.7]:#np.linspace(0.35, 1, 4, endpoint = True):
        for k in [0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,0.25, 0.26,0.27,0.28,0.29,0.3]: # p values
        #for k in [0.2,0.05,0.09]:
            parameters.extend([[i,j,k]])
            
parameters = np.array(parameters).reshape((806,61,3))
#parameters = np.array(parameters).reshape((26,12,3))
s = np.random.randint(0, 2**31, (806, 61, 500, 2),dtype='int')

#OnePacemakerBeat(parameters, s, 19)
#np.save('parameters', parameters)
#np.save('seeds', s)

