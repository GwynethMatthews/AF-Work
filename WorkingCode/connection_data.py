import pickle
import numpy as np
import matplotlib.pyplot as plt

results = pickle.load(open('connection_data_set1_t0.5_p0.5.p','rb'))

for i in results:
    current_values = i[1][0][0]
    number_of_connections = i[1][0][1]
    for j in i[1][1:]:
        number_of_connections = np.add(number_of_connections, j[1])
    number_of_connections = number_of_connections/250000    
   # print(i[0])
   # print(current_values)
   # print(number_of_connections)
    #print(sum(number_of_connections))
    #if i[0] == 0.46:
    plt.bar(current_values[:10], number_of_connections[:10], width=0.01)
   # plt.scatter(current_values[:10], number_of_connections[:10])