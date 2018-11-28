import Atrium as AC
import numpy as np
import pickle

def Percolation(repeats,nu_para,nu_trans):
    seed1 = 1
    seed2 = 2
    seed3 = 3
    seed4 = 4
    #nu_perc_prob = []
    #nu_perc_values = []
    data_full = []
    for i in nu_para:
        #perc_sum = 0
        nu_para_perc_prob = []
        nu_para_perc_values = []
        print(i)
        for j in nu_trans:
            perc_values = []
            perc_sum = 0
            for k in range(repeats):
                A = AC.Atrium(hexagonal = False, v_para=i,v_tran_1= j,d = 0, seed1=seed1, seed2=seed2, seed3=seed3,seed4 = seed4)
                seed1 += 10
                seed2 += 10
                seed3 += 10
                seed4 += 10
                
                percolation = A.CMP2D_timestep_perc()
                perc_sum += percolation
                perc_values.extend([percolation])
                
            perc_prob = float(perc_sum)/repeats
            nu_para_perc_prob.extend([perc_prob])
            nu_para_perc_values.extend([perc_values])
            #nu_trans_data = [j,perc_prob,perc_values]
            
        #data = [i,nu_trans,nu_para_perc_prob,nu_para_perc_values]    
        data_full.extend([i,nu_trans.tolist(),nu_para_perc_prob,nu_para_perc_values])
        #nu_perc_prob.extend(nu_para_perc_prob)
        #nu_perc_values.extend(nu_para_perc_values)
      
    pickle.dump(data_full,open( "percolation2.p", "wb" ) )


def Percolation2(repeats,nus):
    seed1 = 1
    seed2 = 2
    seed3 = 3
    seed4 = 4
    nu_perc_prob = []
    nu_perc_values = []
    data_full = []
    for i in nus:
        print(i)
        perc_values = []
        perc_sum = 0
        for k in range(repeats):
            A = AC.Atrium(hexagonal = False, v_para=i,v_tran_1= i,d = 0, seed1=seed1, seed2=seed2, seed3=seed3,seed4 = seed4)
            seed1 += 10
            seed2 += 10
            seed3 += 10
            seed4 += 10
            
            percolation = A.CMP2D_timestep_perc()
            perc_sum += percolation
            perc_values.extend([percolation])
        perc_prob = float(perc_sum)/repeats
        nu_perc_prob.extend([perc_prob])
        nu_perc_values.extend([perc_values])
        #nu_trans_data = [j,perc_prob,perc_values]
            
        #data = [i,nu_trans,nu_para_perc_prob,nu_para_perc_values]    
    data_full = [nus,nu_perc_prob,nu_perc_values]
        #nu_perc_prob.extend(nu_para_perc_prob)
        #nu_perc_values.extend(nu_para_perc_values)
      
    
    pickle.dump(data_full,open( "percolation_same_nu_hex_extra.p", "wb" ) )
    
nu_para = np.linspace(0.46, 0.52, 13, endpoint = True)

Percolation2(50, nu_para)

#pickle.dump(results,open( "percolation_0_0.45.p", "wb" ) )            
            
            
            
            