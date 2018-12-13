import Atrium_new as AC
import pickle
import numpy as np

def CMP2D_page67(A):
    """Runs CMP2D() and collects data for ECG and number of excited cells. 
    Used to replicate page67 of Kishan's thesis"""
    np.random.seed(A.seed_prop)
    num_ex_cells = []
    excitation_rate = []
    #ecg_values = []
    time = np.arange(0,A.tot_time)
    
    while A.t < A.tot_time:
        A.CMP2D_timestep2()
        excited_cells = len(A.states[0])
        num_ex_cells.extend([excited_cells])
        excitation_rate.extend([np.mean(A.excitation_rate[A.states[0]])])
        #ECG_value = Atrium.ECG([((Atrium.size/2)+0.5),((Atrium.size/2)+0.5)])
        #ecg_values.extend([ECG_value])                      

    data = [num_ex_cells,excitation_rate,time]    
    pickle.dump(data,open( "excitation_over_time.p", "wb" ) )

A = AC.Atrium(hexagonal = True, model = 2, L = 100, v_para = 0.49,
                 v_tran_1 = 0.49, v_tran_2 = 0.49, d = 0.05,
                 threshold = 0.5, p = 0.75, e = 0.05, rp = 30, tot_time = 10000, 
                 pace_rate = 220, s1 = 100, s2 = 4760, s3 = 3306, s4 = 476)