"""Atrium is the normal model (both Sq and Hex)"""
import Atrium as AC
import numpy as np
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as mpat
from matplotlib import collections

###############################################################################

# Initiating the Atrium

#A = AC.Atrium(hexagonal = False, model = 3, L = 200, v_para = 0.6, v_tran_1 = 0.6,
#               v_tran_2 = 0.6, d = 0.05, threshold_cells = 2, threshold = 0.25,
#               e = 0.05, s1 = 152, s2 = 32522, s3 = 335298, s4 = 74)
A = AC.Atrium(hexagonal = False, model = 1, L = 200, v_para = 1,
                 v_tran_1 = 0.1, v_tran_2 = 0.6, d = 0.05, threshold_cells = 1,
                 threshold = 0.25, e = 0.05, rp = 50, tot_time = 10**6, 
                 pace_rate = 220, s1 = 10, s2 = 4, s3 = 3, s4 = 4)

##############################################################################

# Animation function

def update1(frame_number, mat,A):
    """Next frame update for animation without ECG"""
    if frame_number in A.pace:
        A.SinusRhythm()
    A.Relaxing_ani()
    A.Conduct()
    A.TimeInAF()
    A.t += 1
    # WITH CONVOLUTION
    #convolution = gaussian_filter(A.phases.reshape([A.L, A.L]), sigma=1,
    #                              mode = ('wrap', 'constant'))
    #mat.set_data(convolution)
    
    # WITHOUT CONVOLUTION
    
    #if frame_number in [394,395,396,397,398,399]:
    #    print(A.t)
    #    print(A.phases[23280])
    #    print(A.resting[23280])
        #if 23280 in A.states[0]:
         #   print('fuck')
    A.phases[np.array([5490])] = 3
    data = A.phases.reshape([A.L, A.L])
    mat.set_data(data)
    
    return mat,

def update2(frame_number,collection,A):
    """Next frame update for animation without ECG"""
    if frame_number in A.pace:
        A.SinusRhythm()
    A.Relaxing_ani()
    A.Conduct()
    
    # WITH CONVOLUTION
    convolution = gaussian_filter(A.phases.reshape([A.L, A.L]), sigma=0.5,
                                  mode = ('wrap', 'constant'),cval = A.rp + 1)
    data = np.ravel(convolution)
    collection.set_array(data)
    
    # WITHOUT CONVOLUTION
    #collection.set_array(np.array(A.phases))
    
    return ax,

###############################################################################

#Running the Animation

if A.hexagonal == False:
    np.random.seed(A.seed_prop)
    
    fig1 = plt.figure(figsize = [15,15])
    ax = fig1.subplots(1, 1)
    mat1 = ax.matshow(A.phases.reshape([A.L, A.L]), cmap=plt.cm.jet_r)
    mat1.set_clim(0, A.rp)
    ax.set_axis_off()
    ani1 = animation.FuncAnimation(fig1, update1, frames = A.tot_time
                                   ,fargs = (mat1,A), interval=100, repeat = None)
    plt.axis([0, A.L, 0, A.L])

if A.hexagonal == True:
    np.random.seed(A.seed_prop)

    fig1 = plt.figure(figsize = [15,15])
    ax = fig1.subplots(1,1)
    patches = []
    offsets = []
    a = np.tan(np.pi/6)*0.5
    b = 0.5/np.cos(np.pi/6)
    c = 1-a-b
    for i in range(A.L):
        for j in range(A.L):
            if i%2 == 0 and j%2 == 0:
                offsets.extend([(j+0.5, i-(i*c))]) 
            elif i%2 == 0 and j%2 != 0:
                offsets.extend([(j+0.5, i-(i*c))]) 
            elif i%2 !=0 and j%2 == 0:
                offsets.extend([(j, i-(i*c))]) 
            else:
                offsets.extend([(j, i-(i*c))]) 
    for k in offsets:
        patches.extend([mpat.RegularPolygon(k, 6, radius = 0.5/np.cos(np.pi/6))])
        
    collection = collections.PatchCollection(patches, cmap=plt.cm.jet_r)
    ax.add_collection(collection, autolim=True)
    ax.axis('equal')
    ax.set_axis_off()
    ani1 = animation.FuncAnimation(fig1, update2, frames = A.tot_time
                                   ,fargs = (collection,A), interval=100, repeat = None)
    plt.axis([0, A.L + 1, 0, A.L + 1])