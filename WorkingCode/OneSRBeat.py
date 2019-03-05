
import numpy as np
import sys
import Atrium_Final as AF
job_number = int(sys.argv[1])

input_param = np.load('parameters.npy')
input_seeds = np.load('seeds.npy')

def OnePacemakerBeat(parameters=input_param, seeds=input_seeds, itr=job_number):
    data_full = []
    
    nu = 0.55

    for l in range(len(parameters[itr])):
        repeat_data = []
	#print(l)

        for i in range(100):  ### Number of repeats   
            
            nu = parameters[itr][l][0]
            tau = int(parameters[itr][l][1])
            p = parameters[itr][l][2]

            t_under_on = parameters[itr][l][3]
            pace = tau + parameters[itr][l][4]
            
            A = AF.SourceSinkModel(hexagonal=True, threshold=1, p_nonfire=p, pace_rate=pace,
                       Lx=100, Ly=100, tot_time=10000, nu_para=nu, nu_trans=nu, rp=tau,
                       seed_connections=seeds[itr][l][i][0], seed_prop=seeds[itr][l][i][1], 
                       boundary=True, pacemaker_line=True,radius = 3, charge_conservation = False,
                       t_under = 1, t_under_on = t_under_on)

            AF_start = int((10 * pace) + (2.5 * A.L))# A.tot_time # what is this for?

            np.random.seed(A.seed_prop)

            A.cmp_timestep()   ### With one sinus beat       
            
            while A.stop == False:
                if len(A.states[0]) != 0:
                    
                    if A.t < A.tot_time:
                        #if A.t < 10 * pace:
                        #    A.cmp_timestep()
                    
                        #else:
                            #A.AF = True # Got rid of because it makes the condition below redundant 
                            #AF_start = A.t # This will update every loop so have just set it to int((10 * pace) + (2.5 * A.L))
                        
                        A.cmp_no_sinus() 
                        #A.cmp_no_timestep() # Assumed this was meant to be A.cmp_no_sinus()

                    else:
                        A.fail_safe = True
                        A.time_extinguished = A.t
                        A.stop = True

                else:
                    A.time_extinguished = A.t
                    A.stop = True
                    
                
            # nu, tau, p,whether the charge under threshold is conserved, pace_rate, repeat number, s1, s2,
            # A.fail_safe = whether it extinguishes at before tot_time, A.AF = whether it
            # eneters AF, A.t_AF = how long it was in AF for, 
            # A.time_extinguished = time wave is extinguished == tot_time if doesn't terminate

            data = np.array([parameters[itr][l][0]*100, parameters[itr][l][1]*100, parameters[itr][l][2]*100, 
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
            for i in np.linspace(0, 1, 10, endpoint = True): # nu values
                for k in np.linspace(0, 1, 10, endpoint = True): # p values
                    parameters.extend([[i,j,k,m,n]])
            
parameters = np.array(parameters).reshape((1000,3,5))
s = np.random.randint(0, 2**31, (1000, 90, 100, 2),dtype='int')

np.save('parameters', parameters)
np.save('seeds', s)

