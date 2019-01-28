"""Atrium is the normal model (both Sq and Hex)"""
import Atrium_Final as AC
import numpy as np
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as mpat

from matplotlib import collections


plt.rcParams['animation.ffmpeg_path']

###############################################################################
# Initiating the Atrium
convolve = True

seed1 = 1460926859
seed2 = 623486203
nu = 0.7

A = AC.SourceSinkModel(hexagonal=True, threshold=1, p_nonfire=0.05, pace_rate= 220,
                       L=100, tot_time= 1000, nu_para=nu, nu_trans=nu, rp = 30,
                       seed_connections=seed1, seed_prop=seed2)


###############################################################################
# Animation function

def update_square(frame_number, mat, A, convolve):
    """Next frame update for animation without ECG"""

    #print(A.t)
    # print(A.phases[0])

    if A.t % 100 == 0:
        ani.event_source.stop()
        A.change_connections(1, 1)
        ani.event_source.start()
    if A.t == 6000:
            A.p_nonfire = 0

    A.cmp_animation()

    ###### WITH CONVOLUTION ######

    if convolve:
        convolution = gaussian_filter(A.phases.reshape([A.L, A.L]), sigma=1.4,
                                  mode=('wrap', 'nearest'))
        

        mat.set_data(convolution)
    
    ###### WITHOUT CONVOLUTION ######
    else:
        data = A.phases.reshape([A.L, A.L])
        mat.set_data(data)
    
    return mat,


def update_hex(frame_number, collection, A, convolve):    # Frame number passed as default so needed
    """Next frame update for animation without ECG"""

    A.cmp_animation()

    # WITH CONVOLUTION

    if convolve:
        convolution = gaussian_filter(A.phases.reshape([A.L, A.L]), sigma=1.2,
                                      mode=('wrap', 'nearest'))

    
        data = np.ravel(convolution)
        collection.set_array(data)
        
    # WITHOUT CONVOLUTION
    else:
        collection.set_array(np.array(A.phases))


    ax.set_title('refractory period = %i, threshold = %0.2f, p not fire = %0.3f, \nseed connection = %i, seed propagation = %i, sigma = %0.1f \nnu = %0.3f, t = %i' % (A.rp, A.threshold, A.p_nonfire, A.seed_connections, A.seed_prop, 1.4, A.nu_para, A.t), fontsize=20)
    ax.title.set_position([0.5, 0.85])

    return ax,

###############################################################################

# Running the Animation

if not A.hexagonal:
    np.random.seed(A.seed_prop)
    


    fig1 = plt.figure(figsize=[5, 5])

    ax = fig1.subplots(1, 1)
    ax.tight_layout()

    mat1 = ax.matshow(A.phases.reshape([A.L, A.L]), cmap=plt.cm.jet_r)
    mat1.set_clim(0, A.rp)
    ax.set_axis_off()
    ani = animation.FuncAnimation(fig1, update_square, frames=A.tot_time,
                                  fargs=(mat1, A, convolve), interval=100,
                                  repeat=None)

    plt.axis([0, A.L, 0, A.L])
    plt.show()

if A.hexagonal:
    np.random.seed(A.seed_prop)


    fig1 = plt.figure(figsize = [15,15])
    ax = fig1.subplots(1,1)

    fig1.tight_layout()


    patches = []
    offsets = []
    a = np.tan(np.pi/6)*0.5
    b = 0.5/np.cos(np.pi/6)
    c = 1-a-b
    
    for i in range(A.L):
        for j in range(A.L):
            
            if i % 2 == 0 and j % 2 == 0:
                offsets.extend([(j+0.5, i-(i*c))]) 
                
            elif i % 2 == 0 and j % 2 != 0:
                offsets.extend([(j + 0.5, i - (i * c))])
                
            elif i % 2 != 0 and j % 2 == 0:
                offsets.extend([(j, i - (i * c))])
                
            else:
                offsets.extend([(j, i-(i * c))])
                
    for k in offsets:
        patches.extend([mpat.RegularPolygon(k, 6, radius=0.5/np.cos(np.pi/6))])
        
    collection = collections.PatchCollection(patches, cmap= plt.cm.jet_r)

    ax.add_collection(collection, autolim=True)

    #collection.set_edgecolor('face')
    collection.set_clim(0, A.rp)

    ax.axis('equal')
    ax.set_axis_off()
    # ax.set_title('nu = %f' % A.nu_para)
    ani = animation.FuncAnimation(fig1, update_hex, frames=A.tot_time,
                                  fargs=(collection, A, convolve),
                                  interval=50, repeat=None)

    plt.axis([-1, A.L + 1, -1, A.L + 1])
    plt.show()


#ani.save('Conduction block nu_0.5 p_0.25(p non fire is pinward_current).mp4', fps=25)