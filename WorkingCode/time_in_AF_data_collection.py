import Atrium as AC
import numpy as np


def CMP2D_time_AF(Atrium):
    """Runs CMP2D and collects data on the ammount of time in AF"""
    np.random.seed(Atrium.seed_prop)
    while Atrium.t < Atrium.tot_time:
        Atrium.CMP2D_timestep()
        Atrium.excitations[Atrium.states[0]] += 1
        
        # is in SR
        if Atrium.t_AF != 0 and (Atrium.excitations[Atrium.first_col]).max() == (Atrium.excitations).max():
            Atrium.tot_AF += Atrium.t_AF
            Atrium.t_AF = 0
                
        # is in AF        
        if (Atrium.excitations[Atrium.first_col]).max() < (Atrium.excitations).max():
            Atrium.t_AF += 1
                 
    Atrium.tot_AF += Atrium.t_AF
    print(Atrium.tot_AF)
    
Atrium = AC.Atrium()