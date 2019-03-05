import numpy as np 
import Atrium_Final as AF 

nu = np.linspace(0.4, 1, 31)

sc = np.random.randint(0,2**31,(1000))
sp = np.random.randint(0,2**31,(1000))
data = []
for i in range(len(nu)):
    print(nu[i])
    A = AF.SourceSinkModel(hexagonal=True, threshold=1, p_nonfire=0.1, Lx=70,Ly = 100, 
                       tot_time=2*10**4, nu_para=nu[i], rp = 70, pace_rate = 72,
                       nu_trans=nu[i], seed_connections=sc[i], seed_prop=sp[i])
    current_data = []
    A.cmp_timestep()   ### With one sinus beat       
    AF_start = 30000        
    while A.stop == False:            
        if A.t < A.tot_time:
            
            if A.t < 31 * A.pace_rate:
                # pacing
                A.pacing_with_change_of_rp(time_between_pace_and_change_of_rp = 0,
                         increment = -1)

                A.find_propagation_time() # time till first cell in last column excites for the first time
                
                if A.propagated == True:
                    AF_start = int((31 * A.pace_rate) + (5 * A.propagation_time))
            
            else:
 
                if len(A.states[0]) != 0:   # continues to propagate
                    
                    A.cmp_no_sinus()
                    
                    if A.t > AF_start:
                        A.t_AF += 1
                        current_data.extend([sum(A.inward_current)/len(A.inward_current[A.inward_current != 0])])

                else:
                    # terminates
                    A.stop = True
    
        else:
            # reaches tot_time
            A.stop = True
    data.extend([current_data])