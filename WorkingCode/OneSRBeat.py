
import numpy as np
import sys
import Atrium_Final as AF
job_number = int(sys.argv[1])

# parameters array of [v,t,p] of length 60, i th ejob variable
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
            pace = tau + 1
            
            A = AF.SourceSinkModel(hexagonal=True, threshold=1, p_nonfire=p, pace_rate=pace,
                       L=100, tot_time=10000, nu_para=nu, nu_trans=nu, rp=tau,
                       seed_connections=seeds[itr][l][i][0], seed_prop=seeds[itr][l][i][1], boundary=True, pacemaker_line=True)

            AF_start = A.tot_time

            np.random.seed(A.seed_prop)

            A.cmp_timestep()   ### With one sinus beat       
            
            while A.stop == False:
                if len(A.states[0]) != 0:

                    if A.t < 10 * pace:
                        A.cmp_timestep()
                    
                    elif A.t < A.tot_time:
                        if A.t > int((10 * pace) + (2.5 * A.L)):
                            A.AF = True
                            AF_start = A.t
                        
                        A.cmp_no_timestep()

                    else:
                        A.fail_safe = True
                        A.time_extinguished = A.t
                        A.stop = True

                else:
                    A.time_extinguished = A.t
                    A.stop = True
                    
            A.t_AF = int(A.time_extinguished - AF_start)
            
            if A.t_AF > 0:
                A.AF = True
                    
            # nu, tau, p, s1, s2, whether it extinguishes at 5000, whether it
            # eneters AF, when it eneters AF (0 if not AF), time wave is extinguished
            # (5000 if fail safe, 0 if enteres AF)

            data = np.array([parameters[itr][l][0]*100, parameters[itr][l][1]*100, parameters[itr][l][2]*100,
                             i, A.seed_connections, A.seed_prop,
                             A.fail_safe, A.AF, A.t_AF, A.time_extinguished], dtype = 'int')
    
            repeat_data.extend([data])

        data_full.extend([repeat_data])

    data_full = np.array(data_full, dtype = 'int')
    np.save('onesr_data_'+str(itr),data_full)


OnePacemakerBeat(parameters=input_param, seeds=input_seeds, itr=job_number)

parameters = []
for i in np.linspace(0.4, 1, 100, endpoint = True): # nu values
    for j in np.array([40, 50, 60, 70, 80, 90, 100, 110, 120], dtype=int): # tau values
        for k in np.linspace(0, 0.3, 100, endpoint = True): # p values
            parameters.extend([[i,j,k]])
            
parameters = np.array(parameters).reshape((1000,90,3))
s = np.random.randint(0, 2**31, (1000, 90, 100, 2),dtype='int')

np.save('parameters', parameters)
np.save('seeds', s)