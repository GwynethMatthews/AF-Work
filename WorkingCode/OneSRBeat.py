
import numpy as np
import sys
import Atrium_Final as AF
job_number = int(sys.argv[1])

# parameters array of [v,t,p] of length 60, i th ejob variable
input_param = np.load('parameters.npy')
input_seeds = np.load('seeds.npy')

def OnePacemakerBeat(parameters=input_param, seeds=input_seeds, itr=job_number):
    data_full = []

    for l in range(len(parameters[itr])):


        repeat_data = []
	#print(l)
        for i in range(100):  ### Number of repeats

    	    

            A = AF.SourceSinkModel(hexagonal = True, L = 200, v_para = parameters[itr][l][0],
                         v_tran_1 = parameters[itr][l][0], v_tran_2 = parameters[itr][l][0],
                         threshold = parameters[itr][l][1], p = parameters[itr][l][2], rp = 50, tot_time = 10**6,
                         pace_rate = 220, s1 = seeds[itr][l][i][0], s2 = seeds[itr][l][i][1], s3 = seeds[itr][l][i][2])

            np.random.seed(A.seed_prop)

            A.SinusRhythm()
            A.Relaxing()
            A.Conduct()
            while A.stop == False:

                if len(A.states[0]) != 0:

                    if A.t != 5000:
                        A.Relaxing()
                        A.Conduct()
                        A.t += 1

                        A.TimeInAF()

                        if A.AF == True:
                            A.t_AF = A.t
                            A.stop = True


                    else:

                        A.fail_safe = True
                        A.time_extinguished = A.t
                        A.stop = True


                else:
                    A.time_extinguished = A.t
                    A.stop = True
            # nu, threshold, p, s1, s2, s3, whether it extinguishes at 5000, whether it
            # eneters AF, when it eneters AF (0 if not AF), time wave is extinguished
            # (5000 if fail safe, 0 if enteres AF)

            data = np.array([parameters[itr][l][0]*100, parameters[itr][l][1]*100, parameters[itr][l][2]*100,
                             i, A.seed_connect_tran, A.seed_connect_para, A.seed_prop,
                             A.fail_safe, A.AF, A.t_AF, A.time_extinguished], dtype = 'int')
            repeat_data.extend([data])

        data_full.extend([repeat_data])

    data_full = np.array(data_full, dtype = 'int')
    np.save('onesr_data_'+str(itr),data_full)


OnePacemakerBeat(parameters=input_param, seeds=input_seeds, itr=job_number)

#parameters = []
#for i in np.linspace(0.01, 1, 100, endpoint = True): # nu values
#    for j in np.array([0.25,0.3,0.4,0.45,0.5,0.75]): # threshold values
#        for k in np.linspace(0.01, 1, 100, endpoint = True): # p values
#            parameters.extend([[i,j,k]])
#parameters = np.array(parameters).reshape((1000,60,3))
#s = np.random.randint(0, 2**31, (1000, 60, 100, 3),dtype='int')
#
#np.save('parameters',parameters)
#np.save('seeds',s)