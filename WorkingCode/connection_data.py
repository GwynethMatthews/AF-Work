import pickle
import numpy as np
import matplotlib.pyplot as plt

results = pickle.load(open('connection_data_set1_t0.5_p0.5.p','rb'))
average_number = 0
for i in results:
    current_values = i[1][0][0]
    number_of_excited_cells = i[1][0][1]
    for j in i[1][1:]:
        number_of_excited_cells = np.add(number_of_excited_cells, j[1])
    number_of_excited_cells = number_of_excited_cells/250000    
   # print(i[0])
#    print(current_values)
#    print(number_of_excited_cells)
#    print(sum(number_of_excited_cells[15:]))
    
    average_number += sum(number_of_excited_cells[:15])/float(sum(number_of_excited_cells))
    #if i[0] == 0.46:
    plt.pie(number_of_excited_cells[:10])
   # plt.scatter(current_values[:10], number_of_connections[:10])
print(average_number/len(results))