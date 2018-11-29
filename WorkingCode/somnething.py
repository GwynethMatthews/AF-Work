import Atrium_simplest as AC
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

cmap = mcolors.LinearSegmentedColormap.from_list('cmap',[(1,1,1), (0,0,0), (1,0,0), (0,1,0), (0,0,1)])

L = 200
A = AC.Atrium(L=200,v=0.1,d=0.05,e=0.05,rp=50,tot_time = 10**6
                 ,pace_rate = 220,seed1 = 1,seed2=2,seed3=3)

ups = np.isnan(A.n_up)
downs = np.isnan(A.n_down)
print(A.n_up)
print(ups)
print(A.n_down)
print(downs)

not_connected_up = np.where(ups == True)[0]
not_connected_down = np.where(downs == True)[0]
not_connected_both = []
not_connectu = []
not_connectd = []
for i in not_connected_up:
    if i in not_connected_down:
        not_connected_both.extend([i])
    else:
        not_connectu.extend([i])
for i in not_connected_down:
    if i not in not_connected_up:
        not_connectd.extend([i])
    
        
#print(not_connected)
connections = np.full((L*L), fill_value = 0)
connections[not_connected_both] = 1
connections[not_connectd] = 2
connections[not_connectu] = 3
connections[A.dysfunctional_cells] = 4

fig1 = plt.figure(figsize = [15,15])
ax = fig1.subplots(1, 1)
mat1 = ax.matshow(connections.reshape([A.size, A.size]), cmap=cmap)
mat1.set_clim(0, 4)
#white not connected up or down
# black is connected in at lest one direction
