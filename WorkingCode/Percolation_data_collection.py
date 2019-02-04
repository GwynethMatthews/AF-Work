import Atrium_Final as AC
import numpy as np

def CMP2D_timestep_perc(Atrium):
    """A single timestep"""
    Atrium.ectopic_beat([4950,4951,5049,5050,5051,5150,5151])
    Atrium.cmp_no_sinus()
    #c = 0
    while len(Atrium.states[0]) != 0:
        #if len(Atrium.receive_current):
        #c += Atrium.receive_current # len(Atrium.states[0])/float(len(Atrium.receive_current))
        #print(c)
        #print(Atrium.states[0])
        Atrium.cmp_no_sinus()
        Atrium.t += 1
        #c += len(Atrium.states[0])
        a = sum(Atrium.number_of_excitations[Atrium.first_col])
        b = sum(Atrium.number_of_excitations[Atrium.last_col])
        
        if a > 0:

            if b > 0:
                print(Atrium.receive_current )
                return float(Atrium.receive_current )/Atrium.t

    return 0

def Percolation_all_data(nu, seeds):
    data_all = [] 
    for i in nu:
       # nu_data = []
        print(i)
        perc_time = 0
        for k in range(50):
            print(k)
            A = AC.SourceSinkModel(hexagonal=True, threshold=1, p_nonfire=0.05, pace_rate= 220,
                       L=200, tot_time= 10000, nu_para=i, nu_trans=i, rp = 50,
                       seed_connections=seeds[int(i)][int(k)][0], seed_prop=seeds[int(i)][int(k)][1], boundary = True, pacemaker_line = True, radius = 3)
            np.random.seed(A.seed_prop)
            perc_time += CMP2D_timestep_perc(A)
        nu_perc_time = perc_time/50
       
        data_all.extend([[i, nu_perc_time]])
        
    data_full = np.array(data_all)
    print(data_full)
    #np.save('p1_fraction_of_neighbours_excited_out_of_all_resting_neighbours', data_full)
    #np.save('p1_average_number_of_excited_cells', data_full)
   #np.save('p0.1_average_time_to_percolation', data_full)
    np.save('p0.05_average_inward_current_over_excited_cells', data_full)

nu = np.linspace(0, 1, 101, endpoint = True)

s = np.random.randint(0, 2**31, (101,50,2),dtype='int')

Percolation_all_data(nu, s)