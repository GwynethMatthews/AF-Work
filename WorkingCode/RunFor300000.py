#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 16:26:38 2019

@author: Jack
"""

import numpy as np
import Atrium_Final as AF
import sys

job_number = int(sys.argv[1])
#
#
#full_data = []
#
#for i in range(832):
#    full_data.append(np.load('PhaseData/onesr_data_' + str(i) + '.npy'))
#    
#print(full_data[300][2][0])
    
#print(full_data[344][6][0][7])    ### 7th value is if it hits fail safe

#data_num = 0
#
#run_3k = []
#
#while data_num < 2000:
#    
#    p_rand = np.random.randint(1, 10) # Vals from 0.01 to 0.09
#    job_num = np.random.randint(100, 416) * 2    # Keeps only low p values
#    repeat_number = np.random.randint(0, 100)
#    
#    data = full_data[job_num][p_rand][repeat_number]
#    
#    if data[6]:  
#        data_num += 1
#        
#        valid_data = [round(data[0])/100, round(data[1]), round(data[2])/100, round(data[4]), round(data[5])]
#        
#        run_3k.append(valid_data)
#    
#np.save('run_3k.npy', run_3k)         

def decreasing_p_large_time(run_3k, itr):
    
    two_params_data = []
    
    for k in range(2):       # 2 Param sets per job

        params = run_3k[(2 * itr) + k]          # Itr k does k = 2 * itr and k + 1
        
        nu = params[0]
        tau = int(params[1])
        p = params[2]
        seed1 = int(params[3])
        seed2 = int(params[4])
        
        pace = tau + 2
        
        full_data = []
               
        for i in range(2):       # i == 0 is decreasing p, i == 1 is constant p
            print(i)
        
            avg_resting = None
            std_resting = None
            med_resting = None
            max_resting = None
            min_resting = None
            position_of_min = None
            position_of_max = None
            Fraction_of_resting_cells_last_timestep = None
            
            regular_resting_list = []
                    
            A = AF.SourceSinkModel(hexagonal=True, threshold=1, p_nonfire=p, pace_rate=pace,
                                   Lx=70, Ly=100, tot_time=500000, nu_para=nu, nu_trans=nu, rp=tau,
                                   seed_connections=seed1, seed_prop=seed2, 
                                   charge_conservation = False, t_under = 1, t_under_on = False)
        
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
                            resting_cells_at_last_beat = int(len(A.resting[A.resting == True]))
                            print('here')
     
                        if len(A.states[0]) != 0:   # continues to propagate
                            
                            if i == 0 and A.t % 2000 == 0 and A.t > 15000 and A.p_nonfire > 0.0006:
                                A.p_nonfire -= 0.0003
                                
                            if A.t % 2000 == 0:
                                regular_resting_cells = int(len(A.resting[A.resting == True]))
                                regular_resting_list.append(regular_resting_cells)
                        
                            A.cmp_no_sinus()
                            
                            if A.t > AF_start:
                                A.t_AF += 1
                                A.resting_cells_over_time_collect()
        
                        else:
                            # terminates
                            A.time_extinguished = A.t
                            A.stop = True
                            
                            Fraction_of_resting_cells_last_timestep = int(len(A.resting[A.resting == True]))
        
                            
                else:
                    # reaches tot_time
                    A.fail_safe = True
                    A.time_extinguished = A.t
                    A.stop = True
        
                    Fraction_of_resting_cells_last_timestep = int(len(A.resting[A.resting == True]))
                    
            
            if A.t_AF > 0:
                A.AF = True
            
            if len(A.resting_cells_over_time) > 200:     ### Want to ignore last 200 values
            
                resting_cells_minus_slice = np.array(A.resting_cells_over_time[:-200])
            
                avg_resting = np.mean(resting_cells_minus_slice)    
                std_resting = np.std(resting_cells_minus_slice)
                med_resting = np.median(resting_cells_minus_slice)
                max_resting = max(resting_cells_minus_slice)
                min_resting = min(resting_cells_minus_slice)
                
                position_of_min_slice = np.where(resting_cells_minus_slice == min(resting_cells_minus_slice))[0][0] + AF_start # time of min
                position_of_max_slice = np.where(resting_cells_minus_slice == max(resting_cells_minus_slice))[0][0] + AF_start # time of max
                position_of_min = np.where(np.array(A.resting_cells_over_time) == min(A.resting_cells_over_time))[0][-1] + AF_start # time of min
                position_of_max = np.where(np.array(A.resting_cells_over_time) == max(A.resting_cells_over_time))[0][-1] + AF_start
                
            num_readings = int(np.floor((A.time_extinguished - AF_start)/5000.)) 
            time_splits = np.array_split(A.resting_cells_over_time[-(num_readings * 5000):], num_readings)
            mov_sum = np.array(list(map(sum, time_splits)), dtype = int)
                # nu
                # tau
                # p
                # whether the charge under threshold is conserved
                # pace_rate
                # nu_para/nu_trans
        #     repeat number
        #     seed_connection
        #     seed_propagation
        #     A.fail_safe = whether it extinguishes at tot_time
        #     A.AF = whether it eneters AF
        #     A.t_AF = how long it was in AF for
        #     A.time_extinguished = time wave is extinguished == tot_time if doesn't terminate
        #     time of AF starting
        #     average number of resting cells
        #     std of number of resting cells
        #     median of number of resting cells
        #     min of number of resting cells
        #     max of number of resting cells
        #     fraction of resting cells at last beat
        #     fraction of resting cells at last time step either self.states[0] == 0 or self.t == tot_time
        #     position of min
        #     position of max
            
            data = np.array([nu*100, tau, p*100,    #[0] # 0, 1, 2
                             A.seed_connections, A.seed_prop,       # 3, 4
                             A.fail_safe, A.AF, A.t_AF, A.time_extinguished, AF_start,   # 5, 6, 7, 8, 9
                             avg_resting, std_resting, med_resting, min_resting, max_resting,   # 10, 11, 12, 13, 14
                             resting_cells_at_last_beat, Fraction_of_resting_cells_last_timestep,   # 15, 16
                             position_of_min_slice, position_of_max_slice, position_of_min, position_of_max], dtype = int)    # 17, 18, 19, 20
        
            data = np.array([data, regular_resting_list, mov_sum])        # [1] resting cells at p changes, [2] sum every 5000
        
            full_data.extend([data])
        
        two_params_data.extend([full_data])
       
    two_params_data = np.array(two_params_data)
    
    np.save('long_run_decrease_p_' + str(itr) + '.npy', two_params_data)
 
#    
run_3k = np.load('run_3k.npy')
decreasing_p_large_time(run_3k, job_number)

#
#for i in range(50):
#    print(i)
#    decreasing_p_large_time(run_3k, i)