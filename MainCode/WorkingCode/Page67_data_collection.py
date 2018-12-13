import Atrium as AC
import pickle
import numpy as np

def CMP2D_page67(Atrium):
    """Runs CMP2D() and collects data for ECG and number of excited cells. 
    Used to replicate page67 of Kishan's thesis"""
    np.random.seed(Atrium.seed_prop)
    num_ex_cells = []
    ecg_values = []
    time = np.arange(0,Atrium.tot_time)
    
    while Atrium.t < Atrium.tot_time:
        Atrium.CMP2D_timestep()
        excited_cells = len(Atrium.states[0])
        num_ex_cells.extend([excited_cells])
        ECG_value = Atrium.ECG([((Atrium.size/2)+0.5),((Atrium.size/2)+0.5)])
        ecg_values.extend([ECG_value])                      

    data = [num_ex_cells,ecg_values,time]    
    pickle.dump(data,open( "data_page67.p", "wb" ) )

Atrium = AC.Atrium()