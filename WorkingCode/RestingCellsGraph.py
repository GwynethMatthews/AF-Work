import numpy as np
import matplotlib.pyplot as plt

min_AF_time = 500
bins = 50

data = []

for i in np.arange(832):
    results = np.load('onesr_data_%i.npy' %i)
    data.extend(results)

AF_probs = []
avg_resting_cells = []   

for i in data:
    for j in i:
        if j[7] == 1 and j[8]> min_AF_time:

            if j[6] == 0:
                AF_probs.extend([1])
               
            if j[6] == 1:
                AF_probs.extend([0])
         
            avg_resting_cells.extend([j[11]])

zipped= list(zip(avg_resting_cells,AF_probs))
zipped.sort()

a = np.array_split(zipped,bins)

x = []
y = []
for i in a:
    unzip1 = []
    unzip2 = []
    for j in i:
        unzip1.extend([j[0]])
        unzip2.extend([j[1]])
    
    x.extend([sum(unzip1)/len(unzip1)])
    y.extend([sum(unzip2)/len(unzip2)])
    
    
plt.scatter(x,y,marker = 'x',s = 25)
plt.xlabel('Fraction of Resting Cells', fontsize = 24)
plt.ylabel('Probability of Termination', fontsize = 24)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)