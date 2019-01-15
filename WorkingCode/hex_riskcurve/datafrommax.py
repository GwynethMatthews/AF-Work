import numpy as np 
data = np.concatenate([np.load('onesr_data_0.npy'),np.load('onesr_data_1.npy'), np.load('onesr_data_2.npy')])

results2 = []
for i in data:
    results1 = []
    for j in i:
        #print(j)
        if j[1] == 25:
            results1.extend([j])
    if results1 != []:
        results2.extend([results1])
print(results2[0][0])