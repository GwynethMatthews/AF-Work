#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 10:03:26 2019

@author: Jack
"""

import numpy as np
import matplotlib.pyplot as plt

parameters = np.load('parameters.npy')

#### 20 p values so parameters[i even: i + 1] gives one full spectrum of p for 1 nu value

#### 16 nu values so parameters[i even: i + 31] gives one full spectrum of nu for 1 tau value

#full_data = np.load('full_data.npy')

full_data = []

for i in range(len(parameters)):
    full_data.append(np.load('PhaseData/onesr_data_' + str(i) + '.npy'))
    
#print(full_data[344][6][0][7])    ### 7th value is if it hits fail safe

for i in range(10):
    print(full_data[344][6][0][i])


    
#print(full_data[0][0][0][11])     ### Avg time in AF

refrac_tissue_frac_data = []   # Treating each repeat separately
prob_sustained_af = []     # Averaging over the 100 repeats to say which sustain and which don't (no way to define this prob from one repeat)

tau_vals_refrac_tissue = []    ### This is ugly but needs more tau values as not averaging over repeats
tau_vals_prob_af = []

tau_list = np.arange(50, 102, 2)
avg_refrac_tissue_for_tau_list = []

tau_af_count = 0
avg_refrac_tissue_for_tau = 0

for i in range(len(full_data)):    # Job number   0 - 832             
    for j in range(len(full_data[i])):      # p val within job   0 - 10
        
        sustained_af = 0
        af_count = 0
        
        for k in range(len(full_data[i][j])):      # Run number 0 - 100
            
            val = full_data[i][j][k][11]    #### Change 11 to do the same for other vals
                # Reset count of af for calculation of fraction of sustained af
            
            if val != None:
                af_count += 1
                tau_af_count += 1
                
                refrac_tissue_frac_data.append(val)        # Only append if there was AF
                avg_refrac_tissue_for_tau += val
                
                tau_vals_refrac_tissue.append(full_data[i][j][k][1])

                sustained_af += full_data[i][j][k][6]     # Check if there was AF above so add True (1) or False (0) if it was sustained or not
                
        
        if af_count > 0:
            
            prob = float(sustained_af) / af_count
            
            prob_sustained_af.append(prob)
            tau_vals_prob_af.append(full_data[i][j][0][1])
            
    if (i + 1) % 32 == 0:
        avg_refrac_tissue_for_tau_list.append(float(avg_refrac_tissue_for_tau) / tau_af_count)
        
#        avg_sustained_af = 
        
        tau_af_count = 0
        avg_refrac_tissue_for_tau = 0
            
    

            
            
#print(len(refrac_tissue_frac_data))
#print(len(prob_sustained_af))
#
#plt.figure(0)
#plt.plot(tau_vals_refrac_tissue, refrac_tissue_frac_data, 'rx')
#plt.plot(tau_list, avg_refrac_tissue_for_tau_list, 'bx')
#plt.ylabel('Fraction of Refractory Tissue')
#
#plt.figure(1)
#plt.plot(tau_vals_prob_af, prob_sustained_af, 'bx')
#plt.ylabel('Probability of Sustained AF')
#
#plt.figure(2)
#plt.plot(tau_list, avg_refrac_tissue_for_tau_list, 'x')


#print(np.shape(x[0][0]))

#print(x)