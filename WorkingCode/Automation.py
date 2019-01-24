import Atrium_new_Jack as AC
import numpy as np
values_to_check = []
sc = np.random.randint(0,2**31,(1000))
sp = np.random.randint(0,2**31,(1000))

for i in range(30):
    print(i)
    A = AC.SourceSinkModel(hexagonal=True, threshold=1, p_nonfire=0.01, L=100, 
                       tot_time=3*10**5, nu_para=0.4, rp = 30,
                       nu_trans=0.4, seed_connections=1626465539, seed_prop=1652879485)
    np.random.seed(A.seed_prop)   # Sets seed for all dysfunctional firings etc.
    number_of_episodes = 0
    while A.t < A.tot_time:
        
        if A.AF == True:
            r = 1
        if A.AF == False:
            r = 0
        A.cmp_timestep()
        
        if A.AF == True:
            number_of_episodes +=1
            #print(A.t)
        if A.AF == False and number_of_episodes > 1 and r ==1:
            values_to_check.extend([A.nu_para, A.threshold, A.p_nonfire, A.rp, A.seed_connections, A.seed_prop, A.t])
            print([A.nu_para, A.threshold, A.p_nonfire, A.rp, A.seed_connections, A.seed_prop, A.t])
        if A.t % 150 ==0 and A.nu_para < 1: #,22000,42000,62000,82000]:
            A.change_connections(A.nu_para + 0.01, A.nu_trans + 0.01)
            