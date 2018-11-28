"""Atrium is the nuormal model (both Sq and Hex)"""
#import Atrium as AC1
"""Atrium2 is the excites if connected to enough cells model (both Sq and Hex)"""
#import Atrium2 as AC2
"""Atrium3 is the excites if the sum of all the 1/N values is greater than the threshold model (both Sq and Hex)"""
import Atrium3 as AC3
import numpy as np
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as mpat
from matplotlib import collections

###############################################################################

# Initiating the Atrium

#A = AC1.Atrium(hexagonal = False, v_para=0.6,v_tran_1= 0.6,v_tran_2= 0.4, e = 0.05, seed1=1520, seed2=250, seed3=230,seed4 = 204)
#A = AC2.Atrium(hexagonal = True, v_para=0.42,v_tran_1= 0.42,v_tran_2= 0.42, threshold = 2, e = 0.05, seed1=120, seed2=25, seed3=20,seed4 = 20)
A = AC3.Atrium(hexagonal = True, v_para=0.37,v_tran_1= 0.37,v_tran_2= 0.37, threshold = 0.2, seed1=1120, seed2=2024, seed3=23343)

##############################################################################

# Animation function

def update1(frame_number, mat,A):
    """Next frame update for animation without ECG"""
    if frame_number in A.pace:
        A.SinusRhythm()
    A.Relaxing_ani()
    A.Conduct()
    
    # WITH CONVOLUTION
    convolution = gaussian_filter(A.phases.reshape([A.size,A.size]), sigma=1,mode = ('wrap','constant'))
    mat.set_data(convolution)
    
    # WITHOUT CONVOLUTION
    #data = A.phases.reshape([A.size,A.size])
    #mat.set_data(data)
    
    return mat,

def update2(frame_number,collection,A):
    """Next frame update for animation without ECG"""
    if frame_number in A.pace:
        A.SinusRhythm()
    A.Relaxing_ani()
    A.Conduct()
    
    # WITH CONVOLUTION
    #convolution = gaussian_filter(A.phases.reshape([A.size,A.size]), 
    #                              sigma=0.65,mode = ('wrap','constant'),cval = A.rp)
    #data = np.ravel(convolution)
    #collection.set_array(data)
    
    # WITHOUT CONVOLUTION
    collection.set_array(np.array(A.phases))
    
    return ax,

###############################################################################

#Running the Animation

if A.hexagonal == False:
    np.random.seed(A.seed_prop)
    
    fig1 = plt.figure(figsize = [15,15])
    ax = fig1.subplots(1,1)
    mat1 = ax.matshow(A.phases.reshape([A.size,A.size]),cmap=plt.cm.jet_r)
    mat1.set_clim(0,A.rp)
    ax.set_axis_off()
    ani1 = animation.FuncAnimation(fig1, update1, frames = A.tot_time
                                   ,fargs = (mat1,A), interval=100, repeat = None)
    plt.axis([0,A.size,0,A.size])

if A.hexagonal == True:
    np.random.seed(A.seed_prop)

    fig1 = plt.figure(figsize = [15,15])
    ax = fig1.subplots(1,1)
    patches = []
    offsets = []
    a = np.tan(np.pi/6)*0.5
    b = 0.5/np.cos(np.pi/6)
    c = 1-a-b
    for i in range(A.size):
        for j in range(A.size):
            if i%2 == 0 and j%2 == 0:
                offsets.extend([(j+0.5,i-(i*c))]) 
            elif i%2 == 0 and j%2 != 0:
                offsets.extend([(j+0.5,i-(i*c))]) 
            elif i%2 !=0 and j%2 == 0:
                offsets.extend([(j,i-(i*c))]) 
            else:
                offsets.extend([(j,i-(i*c))]) 
    for k in offsets:
        patches.extend([mpat.RegularPolygon(k,6,radius = 0.5/np.cos(np.pi/6))])
        
    collection = collections.PatchCollection(patches, cmap=plt.cm.jet_r)
    ax.add_collection(collection, autolim=True)
    ax.axis('equal')
    ax.set_axis_off()
    ani1 = animation.FuncAnimation(fig1, update2, frames = A.tot_time
                                   ,fargs = (collection,A), interval=100, repeat = None)
    plt.axis([0,A.size+1,0,A.size+1])