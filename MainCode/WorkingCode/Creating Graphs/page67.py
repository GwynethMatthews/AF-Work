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
f1, (ax5, ax6,ax7) = plt.subplots(3, 1, sharey= 'row')
ax5.plot(results[3],below1,'g',results[3],above1,'r')
ax5.plot(results[3],threshold, '--')
ax6.plot(results[3],below2,'g',results[3],above2,'r')
ax5.set_ylim(0,900)
ax5.set_xlim(0,10000)
ax7.set_xlabel('Timesteps', fontsize=22)
ax6.set_ylabel('Electrogram', fontsize=22)
ax6.set_xlim(0,10000)
ax7.set_xlim(0,10000)
ax5.set_ylabel('Number of \n excited cells', fontsize=22, multialignment='center')
ax7.plot(results[3],below3,'g')
ax7.plot(results[3],above3,'r')
ax7.set_ylim(0,1)
ax7.fill_between(results[3],result_made_up_above,result_made_up_below,AF,color ='r')
ax7.fill_between(results[3],result_made_up_above,result_made_up_below,SR,color ='g')
ax7.text(750, 0.4, "SR",
         horizontalalignment='center', fontsize=40)
ax7.text(2500, 0.4, "AF",
         horizontalalignment='center', fontsize=40)
ax7.text(4000, 0.4, "SR",
         horizontalalignment='center', fontsize=40)
ax7.text(6000, 0.4, "AF",
         horizontalalignment='center', fontsize=40)
ax7.text(8850, 0.4, "SR",
         horizontalalignment='center', fontsize=40)
ax7.set_yticks([])#ax3.set_xaxis('off')
ax5.set_xticklabels([])
ax6.set_xticklabels([])
ax7.tick_params(axis='x', labelsize=20)
ax6.tick_params(axis='y', labelsize=20)
ax5.tick_params(axis='y', labelsize=20)