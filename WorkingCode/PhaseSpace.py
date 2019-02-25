
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
    
    for l in range(len(parameters[itr])):
        repeat_data = []
        print(l)
        
        for i in range(2):  ### Number of repeats   

            #print(parameters[itr][l])
            nu_para = parameters[itr][l][0]
            nu_trans = nu_para/parameters[itr][l][5] # either a third or a quarter of nu_para

            tau = int(parameters[itr][l][1])
            p = parameters[itr][l][2]
            t_under_on = parameters[itr][l][3]
            pace = tau + parameters[itr][l][4]
            
            A = AF.SourceSinkModel(hexagonal=True, threshold=1, p_nonfire=p, pace_rate=pace,
                       Lx=100, Ly=100, tot_time=15000, nu_para=nu_para, nu_trans=nu_trans, rp=tau,
                       seed_connections=seeds[itr][l][i][0], seed_prop=seeds[itr][l][i][1], 
                       charge_conservation = False, t_under = 1, t_under_on = t_under_on)

            AF_start = 20 * pace
            Number_of_resting_cells = A.Ly*A.Lx + 5

            np.random.seed(A.seed_prop)

            A.cmp_timestep()   ### With one sinus beat       
            
            while A.stop == False:
                if A.t < A.tot_time:
                    if A.t < 10 * pace:
                        A.cmp_timestep()

                        A.find_propagation_time()
                        if A.propagated == True:
                            AF_start = int((10 * pace) + (2 * A.propagation_time))
                    
                    else:
                        if len(A.states[0]) != 0:                        
                            A.cmp_no_sinus()
                            
                            if A.t > AF_start:
                                A.t_AF += 1
                                if A.t == AF_start + 1:
                                    print('here')
                                #print('AF')
                                #Fraction_of_resting_cells = len(A.resting[A.resting == True])/float(A.Lx*A.Ly)
                        else:
                            A.time_extinguished = A.t
                            A.stop = True
                            print(len(A.resting[A.resting == True]))
                            Number_of_resting_cells = len(A.resting[A.resting == True])
                            
                            
                else:
                    A.fail_safe = True
                    A.time_extinguished = A.t
                    A.stop = True
                    print(len(A.resting[A.resting == True]))
                    Number_of_resting_cells = len(A.resting[A.resting == True])
                    
            #A.t_AF = int(A.time_extinguished - AF_start)
            
            if A.t_AF > 0:
                A.AF = True
                    
            # nu
            # tau
            # p
            # whether the charge under threshold is conserved
            # pace_rate
            # nu_para/nu_trans
            # repeat number
            # seed_connection
            # seed_propagation
            # A.fail_safe = whether it extinguishes at tot_time
            # A.AF = whether it eneters AF
            # A.t_AF = how long it was in AF for
            # A.time_extinguished = time wave is extinguished == tot_time if doesn't terminate
            # time of AF starting
            # Fraction of cells at time of termination 

            data = np.array([parameters[itr][l][0]*100, parameters[itr][l][1], parameters[itr][l][2]*100, 
                             parameters[itr][l][3], A.pace_rate, parameters[itr][l][5], i, A.seed_connections, A.seed_prop,
                             A.fail_safe, A.AF, A.t_AF, A.time_extinguished, AF_start, Number_of_resting_cells], dtype = 'int')

    
            repeat_data.extend([data])

        data_full.extend([repeat_data])

    data_full = np.array(data_full, dtype = 'int')
    np.save('onesr_data_'+str(itr),data_full)


#OnePacemakerBeat(parameters=input_param, seeds=input_seeds, itr=job_number)

parameters = []
# nu, tau, p, whether the charge under threshold is conserved, amount to add to tau to get pace_rate 
for m in [True, False]:
    for j in np.array([50, 70, 90, 110, 130], dtype=int): # tau values
        for n in [2,10,220-j]:
            for i in np.linspace(0.35, 1, 27, endpoint = True): # nu values
            #for i in np.linspace(0.35, 1, 4, endpoint = True):
                for k in [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.1,0.15,0.16,0.17,0.18,0.19,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1]: # p values
                #for k in [0,0.5,0.9]:
                    for l in [3,4]:
                        parameters.extend([[i,j,k,m,n,l]])
            
parameters = np.array(parameters).reshape((972,40,6))
#parameters = np.array(parameters).reshape((72,10,6))
s = np.random.randint(0, 2**31, (972, 40, 2, 2),dtype='int')
#OnePacemakerBeat(parameters, s, 25)
#np.save('parameters', parameters)
#np.save('seeds', s)

