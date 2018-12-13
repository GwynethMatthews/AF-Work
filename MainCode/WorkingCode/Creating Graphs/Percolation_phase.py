import pickle
import matplotlib.pyplot as plt
import numpy as np
results = pickle.load(open("percolation_phase_sq.p", "rb"))
#levels = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
#levels = [0,0.01,0.5,0.99,1]

fig, ax1 = plt.subplots()
phase_space = ax1.pcolormesh(results[1], results[2], results[3])
plt.colorbar(phase_space, ax = ax1)
ax1.set_title("Percolation Phase Space for the Square CMP Model",fontsize = 20)
ax1.set_xlabel('v_para',fontsize = 20)
ax1.set_ylabel('v_tran',fontsize = 20)  


results1 = pickle.load(open("percolation_graph_nu_hex.p", "rb")) 
fig, ax2 = plt.subplots()
ax2.errorbar(results1[1], results1[2],yerr = results1[3],fmt='x',c = 'b') 
ax2.set_title("Percolation Curve for Hexagonal CMP Model",fontsize = 20)
ax2.set_xlabel('v',fontsize = 20)
ax2.set_ylabel('Probability of Percolation',fontsize = 20)  

results2 = pickle.load(open("percolation_graph_nu_sq.p", "rb"))
fig, ax3 = plt.subplots()
ax3.errorbar(results2[1], results2[2],yerr = results2[3],fmt='x',c = 'r') 
ax3.set_title("Percolation Curve for Square CMP Model",fontsize = 20)
ax3.set_xlabel('v',fontsize = 20)
ax3.set_ylabel('Probability of Percolation',fontsize = 20)

fig, ax4 = plt.subplots()
ax4.errorbar(results1[1], results1[2],yerr = results1[3],fmt='x',c = 'b',label = 'Hexagonal') 
ax4.errorbar(results2[1], results2[2],yerr = results2[3],fmt='x',c = 'r',label = 'Square') 
ax4.set_title("Comparison of Percolation Curve for Square and Hexagonal CMP Models",fontsize = 20)
ax4.set_xlabel('v',fontsize = 20)
ax4.set_ylabel('Probability of Percolation',fontsize = 20)
ax4.legend(fontsize = 20)
y1 = [0, 1]
x1 = [2*np.sin(np.pi/18), 2*np.sin(np.pi/18)]
y2 = [0, 1]
x2 = [0.5, 0.5]
ax4.plot(x1,y1)
ax4.plot(x2,y2)
plt.show()