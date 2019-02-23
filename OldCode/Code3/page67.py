"""Reproduces page 67 of Kishan's thesis"""
import pickle
import numpy as np
import matplotlib.pyplot as plt

results = pickle.load(open("data_page67_0.17.p","rb"))

results_cells = np.array(results[0],dtype = float)
ECG_values = np.array(results[1],dtype = float)

AF = np.array(results[2])
SR = ~AF
below1 = np.ma.masked_where(AF, results_cells)
above1 = np.ma.masked_where(SR, results_cells)
below2 = np.ma.masked_where(AF, ECG_values)
above2 = np.ma.masked_where(SR, ECG_values)
result_made_up_above = np.ones(len(results[3]))
result_made_up_below = np.zeros(len(results[3]))
below3 = np.ma.masked_where(AF, result_made_up_above)
above3 = np.ma.masked_where(SR, result_made_up_above)
threshold = np.full(len(results[3]),fill_value= 220)
f1, (ax1, ax2,ax3) = plt.subplots(3, 1, sharey= 'row')
ax1.plot(results[3],below1,'g',results[3],above1,'r')
ax1.plot(results[3],threshold, '--')
ax2.plot(results[3],below2,'g',results[3],above2,'r')
ax1.set_ylim(0,900)
ax1.set_xlim(0,10000)
ax3.set_xlabel('Timesteps', fontsize=22)
ax2.set_ylabel('Electrogram', fontsize=22)
ax2.set_xlim(0,10000)
ax3.set_xlim(0,10000)
ax1.set_ylabel('Number of \n excited cells', fontsize=22, multialignment='center')
ax3.plot(results[3],below3,'g')
ax3.plot(results[3],above3,'r')
ax3.set_ylim(0,1)
ax3.fill_between(results[3],result_made_up_above,result_made_up_below,AF,color ='r')
ax3.fill_between(results[3],result_made_up_above,result_made_up_below,SR,color ='g')
ax3.text(750, 0.4, "SR",
         horizontalalignment='center', fontsize=40)
ax3.text(2500, 0.4, "AF",
         horizontalalignment='center', fontsize=40)
ax3.text(4000, 0.4, "SR",
         horizontalalignment='center', fontsize=40)
ax3.text(6000, 0.4, "AF",
         horizontalalignment='center', fontsize=40)
ax3.text(8850, 0.4, "SR",
         horizontalalignment='center', fontsize=40)
ax3.set_yticks([])#ax3.set_xaxis('off')
ax1.set_xticklabels([])
ax2.set_xticklabels([])
ax3.tick_params(axis='x', labelsize=20)
ax2.tick_params(axis='y', labelsize=20)
ax1.tick_params(axis='y', labelsize=20)