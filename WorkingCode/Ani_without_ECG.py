"""Atrium is the normal model (both Sq and Hex)"""
import Atrium_new_Jack as AC
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
# AC.Atrium(hexagonal=False, L=200, rp=50, tot_time=10**6, nu_para=0.6, nu_trans=0.6,
#                 pace_rate=220, p_nonfire=0.05, seed_connections=1, seed_prop=4)
A = AC.SourceSinkModel(hexagonal=True, L=100, tot_time=3000, nu_para=0.41, nu_trans=0.41, seed_connections=1, seed_prop=2)
# A = AC.DysfuncModel(hexagonal = True, L = 100,tot_time = 10**6)
# A = AC.Atrium(hexagonal = True, model = 2, L = 100, v_para = 0.6,
#                     v_tran_1 = 0.6, v_tran_2 = 0.5,
#                     threshold = 0.5, p = 0.25, rp = 50, tot_time = 10**6,
#                     pace_rate = 220, s2 = 10, s3 = 40, s4 = 30)
###############################################################################
# Animation function


def update_square(frame_number, mat, A, convolve):
    """Next frame update for animation without ECG"""
    # if A.model == 1:
    #    A.CMP2D_timestep_ani1()
    # else:
    # A.sinus_rhythm()
    print(A.t)
    # print(A.phases[0])
    if A.t == 1000:
        ani.event_source.stop()
        A.change_connections(1, 1)
        ani.event_source.start()

    A.cmp_animation()

    ###### WITH CONVOLUTION ######

    if convolve:
        convolution = gaussian_filter(A.phases.reshape([A.L, A.L]), sigma=1,
                                  mode=('wrap', 'constant'))
        
        mat.set_data(convolution)
    
    ###### WITHOUT CONVOLUTION ######
    else:
        data = A.phases.reshape([A.L, A.L])
        mat.set_data(data)
    
    return mat,


def update_hex(frame_number, collection, A, convolve):    # Frame number passed as default so needed
    """Next frame update for animation without ECG"""
    if A.t == 1000:
        ani.event_source.stop()
        A.change_connections(A.nu_para + 0.08, A.nu_trans + 0.08)
        ani.event_source.start()
        print(A.t)
        print(A.nu_para)
        
    if A.t == 2000:
        ani.event_source.stop()
        A.change_connections(A.nu_para + 0.08, A.nu_trans - 0.08)
        ani.event_source.start()
        print(A.t)
        print(A.nu_para)

    A.cmp_animation()
    
    # WITH CONVOLUTION
    if convolve:
        convolution = gaussian_filter(A.phases.reshape([A.L, A.L]), sigma=1.2,
                                      mode=('wrap', 'constant'), cval=A.rp/2)
    
        data = np.ravel(convolution)
        collection.set_array(data)
        
    # WITHOUT CONVOLUTION
    else:
        collection.set_array(np.array(A.phases))

    ax.set_title('nu = %f' % A.nu_para, fontsize=28)
    
    return ax,

###############################################################################

# Running the Animation


if not A.hexagonal:
    np.random.seed(A.seed_prop)
    
    fig1 = plt.figure(figsize=[15, 15])
    ax = fig1.subplots(1, 1)
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

    fig1 = plt.figure(figsize=[15, 15])
    ax = fig1.subplots(1, 1)
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
        
    collection = collections.PatchCollection(patches, cmap=plt.cm.jet_r)
    ax.add_collection(collection, autolim=True)
    # collection.set_edgecolor('face')
    collection.set_clim(0, A.rp)
    
    ax.axis('equal')
    ax.set_axis_off()
    # ax.set_title('nu = %f' % A.nu_para)
    ani = animation.FuncAnimation(fig1, update_hex, frames=A.tot_time,
                                  fargs=(collection, A, convolve),
                                  interval=5, repeat=None)

    plt.axis([-1, A.L + 1, -1, A.L + 1])
    plt.show()


ani.save('nu_increase_L100_rp30_t0.5_p0.25_sc1_sp2.mp4', fps=30, dpi=250, bitrate=5000)