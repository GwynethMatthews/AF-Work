import Atrium_Final as AC
import numpy as np
values_to_check = []
sc = np.random.randint(0,2**31,(1000))
sp = np.random.randint(0,2**31,(1000))

for i in range(30):
    print(i)
    A = AC.SourceSinkModel(hexagonal=True, threshold=1, p_nonfire=0.17, Lx=70,Ly = 100, 
                       tot_time=5*10**5, nu_para=0.8, rp = 70, pace_rate = 72,
                       nu_trans=0.8, seed_connections=sc[i], seed_prop=sp[i])
    
    np.random.seed(A.seed_prop)   # Sets seed for all dysfunctional firings etc.
    number_of_episodes = 0
    AF_stop = False
    while AF_stop == False:# < A.tot_time:
        if A.t < A.tot_time:
                    
            if A.t < 31 * A.pace_rate:
                # pacing
                A.pacing_with_change_of_rp(time_between_pace_and_change_of_rp = 0,
                         increment = -1)

                A.find_propagation_time() # time till first cell in last column excites for the first time
                
            
            else:

                if len(A.states[0]) != 0:   # continues to propagate
                    
                    if i == 0 and A.t % 2000 == 0 and A.t > 15000 and A.p_nonfire > 0.0006:
                        A.p_nonfire -= 0.0003
                        

                    A.cmp_no_sinus()
                
                elif len(A.states[0] == 0) and A.t > 20000:
                    values_to_check.extend([A.nu_para, A.threshold, A.p_nonfire, A.rp, A.seed_connections, A.seed_prop, A.t])

                    AF_stop = True
                
                else:
                    AF_stop = True
        else:
            AF_stop = True
                    
            
#            if A.AF == True:
#                r = 1
#            if A.AF == False:
#                r = 0
#            A.cmp_timestep()
#            
#            if A.AF == True:
#                number_of_episodes +=1
#                
#            if A.AF == False and number_of_episodes > 0:
#                values_to_check.extend([A.nu_para, A.threshold, A.p_nonfire, A.rp, A.seed_connections, A.seed_prop, A.t])
#                print([A.nu_para, A.threshold, A.p_nonfire, A.rp, A.seed_connections, A.seed_prop, A.t])
#            
#
#   
#            