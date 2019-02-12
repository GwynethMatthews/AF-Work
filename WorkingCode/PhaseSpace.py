
import numpy as np
import sys
import math
import Atrium_Final as AF
job_number = int(sys.argv[1])
#job_number = 280

input_param = np.load('parameters.npy')
input_seeds = np.load('seeds.npy')

#print(np.shape(input_param))

def OnePacemakerBeat(parameters, seeds, itr):
    data_full = []
    
    for l in range(len(parameters[itr])):
        repeat_data = []
        print(l)
        
        for i in range(2):  ### Number of repeats   

            nu = parameters[itr][l][0]
            tau = int(parameters[itr][l][1])
            p = parameters[itr][l][2]
            t_under_on = parameters[itr][l][3]
            pace = tau + parameters[itr][l][4]
            
            A = AF.SourceSinkModel(hexagonal=True, threshold=1, p_nonfire=p, pace_rate=pace,
                       Lx=100, Ly=100, tot_time=10000, nu_para=nu, nu_trans=nu, rp=tau,
                       seed_connections=seeds[itr][l][i][0], seed_prop=seeds[itr][l][i][1], 
                       charge_conservation = False, t_under = 1, t_under_on = t_under_on)

            AF_start = int((10 * pace) + (2.5 * A.Lx))

            np.random.seed(A.seed_prop)

            A.cmp_timestep()   ### With one sinus beat       
            
            while A.stop == False:
                if A.t < A.tot_time:
                    if A.t < 10 * pace:
                        A.cmp_timestep()
                    
                    else:
                        if len(A.states[0]) != 0:                        
                            A.cmp_no_sinus()
                            
                            if A.t > AF_start:
                                A.t_AF += 1
                                #print('AF')
                        else:
                            A.time_extinguished = A.t
                            A.stop = True
                            
                else:
                    A.fail_safe = True
                    A.time_extinguished = A.t
                    A.stop = True
                    
            #A.t_AF = int(A.time_extinguished - AF_start)
            
            if A.t_AF > 0:
                A.AF = True
                    
            # nu
            # tau
            # p
            # whether the charge under threshold is conserved
            # pace_rate
            # repeat number
            # seed_connection
            # seed_propagation
            # A.fail_safe = whether it extinguishes at tot_time
            # A.AF = whether it eneters AF
            # A.t_AF = how long it was in AF for
            # A.time_extinguished = time wave is extinguished == tot_time if doesn't terminate

            data = np.array([np.ceil(parameters[itr][l][0]*1000), parameters[itr][l][1], np.ceil(parameters[itr][l][2]*100), 
                             parameters[itr][l][3], A.pace_rate, i, A.seed_connections, A.seed_prop,
                             A.fail_safe, A.AF, A.t_AF, A.time_extinguished], dtype = 'int')
    
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
                for k in np.linspace(0, 0.95, 20, endpoint = True): # p values
                    parameters.extend([[i,j,k,m,n]])
            
parameters = np.array(parameters).reshape((900,18,5))
s = np.random.randint(0, 2**31, (900, 18, 2, 2),dtype='int')

np.save('parameters', parameters)
np.save('seeds', s)