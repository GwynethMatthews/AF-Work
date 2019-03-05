import Atrium_Final as AC
import numpy as np
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as mpat
from matplotlib import collections

#from pydrive.auth import GoogleAuth
#from pydrive.drive import GoogleDrive
#import os

#gauth = GoogleAuth()
## Try to load saved client credentials
#gauth.LoadCredentialsFile("mycreds.txt")
#
#if gauth.credentials is None:       # Authenticate if they're not there
#    gauth.LocalWebserverAuth()
#    
#elif gauth.access_token_expired:    # Refresh them if expired
#    gauth.Refresh()
#    
#else:         # Initialize the saved creds
#    gauth.Authorize()
#    
#gauth.SaveCredentialsFile("mycreds.txt")    # Save the current credentials to a file
#drive = GoogleDrive(gauth)

###############################################################################
# Initiating the Atrium
#data = []
convolve = True
grey_background = True
resting_cells = False

number_of_beats = 31

nu = 1
tau = 70
p = 0.27
seed1 = 23762876#53517399
seed2 = 2389723#543500048

A = AC.SourceSinkModel(hexagonal=True, threshold=1, p_nonfire=p, pace_rate=tau + 2,
                       Lx=70, Ly=100, tot_time= 11000, nu_para=nu, nu_trans=nu, rp=tau,
                       seed_connections=seed1, seed_prop=seed2, boundary = True, 
                       t_under = 1, t_under_on = False)

np.random.seed(A.seed_prop)

###############################################################################
# Animation function


#### Hex Lattice ####

def update_hex(frame_number, collection, A, convolve):    # Frame number passed as default so needed
    """Next frame update for animation without ECG"""

    if A.t < number_of_beats * A.pace_rate:
        A.pacing_with_change_of_rp_ani(time_between_pace_and_change_of_rp = 0,
                                 increment = -1)

    ### Change rp at end of pacing ###
#    elif A.t == (number_of_beats * A.pace_rate) + 2:
#        A.change_rp(-10)
#        A.cmp_animation()
#        print(A.rp)
#   
    else:
#        if A.t + 1 % 1000 == 0:
#            A.p_nonfire -= 0.002
        
        A.cmp_animation()
#    A.find_propagation_time()

    #print(len(A.states))
    ### CHANGING P_NONFIRE (smaller p_nonfire makes it more likely to fire) ###
#    if A.t in np.arange(number_of_beats * A.pace_rate, number_of_beats * A.pace_rate + 200*7, 200 ):
#        A.p_nonfire -= 0.01
    if len(A.states[0]) == 0:
        print(A.t)
#    ### Note if p is decreasing can have errors when it reaches 0, fixed if set time for p_nonfire = 0 ###    
#        ### CHANGING REFRACTORY PERIOD ###
#    if A.t == 10*A.pace_rate:
#        A.change_rp(130)
#        print(A.rp)
#        
#    ### CHANGING PACE_RATE ###
#    if A.t == 10*A.pace_rate:
#        A.pace_rate = 10**10
#    
#    ### CHANGING P_NONFIRE (smaller p_nonfire makes it more likely to fire) ###
#    if A.t in np.arange(10*110 + 200, 10*110 + 200*49, 200 ):
#        A.p_nonfire -= 0.001
#    ### Note if p is decreasing can have errors when it reaches 0, fixed if set time for p_nonfire = 0 ###
#    if A.t == 10900:
#        A.p_nonfire = 0
    
#    ### Ectopic Beat ###
#    if A.t % A.pace_rate == 0:
#        A.ectopic_beat([4950,4951,5049,5050,5051,5150,5151])
#    A.cmp_animation()
#    global data
#    data.extend([A.excitation_rate])
    #print('her')
    # WITH CONVOLUTION
    
    if convolve:
        if A.boundary == True:
            mode = ('wrap', 'nearest')
            
        elif A.boundary == False:
            mode = ('nearest', 'nearest')
            
        if grey_background:
        
            mx = np.array(A.phases.reshape([A.Ly, A.Lx]) == A.rp)
            mx1 = gaussian_filter(np.ma.masked_array(A.phases.reshape([A.Ly, A.Lx]), mx), sigma=1.2,
                                      mode=mode)
            
            a = max(mx1.flatten())

            mx2 = np.array(mx1 >= 0.9*a)
            mx3 = mx*mx2
            mx1 = np.ma.masked_array(mx1,mx3)

            data = np.ravel(mx1)
        else:
            convolution = gaussian_filter(A.phases.reshape([A.Ly, A.Lx]), sigma=1.2,
                                      mode=mode)
            
            data = np.ravel(convolution)
            
        collection.set_array(data)
        collection.set_clim(0, A.rp)
        
    # WITHOUT CONVOLUTION
    else:
        if grey_background:
            mx = np.array(A.phases.reshape([A.Ly, A.Lx]) == A.rp)
            data = np.ma.masked_array(A.phases,mx)
            collection.set_array(data)
        
        else: 
            collection.set_array(np.array(A.phases))


    ax1.set_title('refractory period = %i, threshold = %0.2f, \nseed connection = %i, seed propagation = %i, pace_rate = %i \nnu = %0.3f, p not fire = %0.3f,\nt = %i' % (A.rp, A.threshold, A.seed_connections, A.seed_prop, A.pace_rate, A.nu_para, A.p_nonfire, A.t), fontsize=20)
    ax1.title.set_position([0.5, 0.85])
    
    if resting_cells == True:
        A.resting_cells = np.roll(A.resting_cells, -1)
        A.resting_cells[-1] = len(A.resting[A.resting == True])
        ax2.clear()

        ax2.plot(A.time_for_graphs, A.resting_cells/(A.Lx * A.Ly))

        ax2.set_xlim(-500,0) 
        #ax2.set_ylim(-400,400)
        ax2.set_xlabel('Time')
        ax2.set_ylabel('Fraction of resting cells')
    
    return ax1,


def update_square(frame_number, mat, A, convolve):
    A.cmp_animation()

    ###### WITH CONVOLUTION ######

    if convolve:
        if A.boundary == False:
            convolution = gaussian_filter(A.phases.reshape([A.Lx, A.Ly]), sigma=1.4,
                                  mode=('wrap', 'nearest'))
        if A.boundary == True:
            convolution = gaussian_filter(A.phases.reshape([A.Lx, A.Ly]), sigma=1.4,
                                  mode=('wrap', 'nearest'))

        mat.set_data(convolution)
    
    ###### WITHOUT CONVOLUTION ######
    else:
        data = A.phases.reshape([A.Lx, A.Ly])
        mat.set_data(data)
    
    return mat,

###############################################################################

# Running the Animation

if not A.hexagonal:
    np.random.seed(A.seed_prop)
    
    fig1 = plt.figure(figsize=[8, 5])

    ax = fig1.subplots(1, 1)
    ax.tight_layout()

    mat1 = ax.matshow(A.phases.reshape([A.Lx, A.Ly]), cmap=plt.cm.jet_r)
    mat1.set_clim(0, A.rp)
    ax.set_axis_off()
    ani = animation.FuncAnimation(fig1, update_square, frames=A.tot_time,
                                  fargs=(mat1, A, convolve), interval=100,
                                  repeat=None)

    plt.axis([0, A.L, 0, A.L])
    plt.show()

if A.hexagonal:
    np.random.seed(A.seed_prop)

    fig1 = plt.figure(figsize = [7,7])
    
    if resting_cells == False:
        ax1 = fig1.subplots(1,1)
        
    if resting_cells == True:
        ax1 = plt.subplot2grid((6,6), (0,0), colspan=5, rowspan=5)
        ax2 = plt.subplot2grid((6,6), (5,0), colspan=5)
        
    patches = []
    offsets = []
    a = np.tan(np.pi/6)*0.5
    b = 0.5/np.cos(np.pi/6)
    c = 1-a-b
    
    for i in range(A.Ly):      # Works as Ly then Lx
        for j in range(A.Lx):
            
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
    
#    pcm = ax.pcolormesh(x, y, Z, vmin=0, vmax=max(A.phases), cmap=plt.cm.jet_r)
    collection = collections.PatchCollection(patches, cmap= plt.cm.jet_r)

    ax1.add_collection(collection, autolim=True)
    
    collection.set_clim(0, A.rp)
    collection.cmap.set_bad((169/600,169/600,169/600))
#    collection.set_edgecolor('face')

    ax1.axis('equal')
    ax1.set_axis_off()
    ani = animation.FuncAnimation(fig1, update_hex, frames=A.tot_time,
                                  fargs=(collection, A, convolve),
                                  interval=50, repeat=None)

    ax1.axis([-1, A.Lx + 1, -1, A.Ly + 1])
    plt.show()

###SAVING VIDEO###


#ani.save('No p change.mp4', fps=60)


#file_path = "NewVid.mp4"
#folder_id = "1zpBUFJO6XAkmoWuWnWGRsj6O6oJCwYI0"      ### Folder ID for AF_Stuff folder
#
#file_save = False
##ani.save(file_path, fps=30, dpi=250, bitrate=5000)
#
#if file_save == True:
#    ani.save(file_path, fps=30, dpi=250, bitrate=5000)
#    
#    file1 = drive.CreateFile({"parents": [{"kind": "drive#fileLink", "id": folder_id}]})
#    file1.SetContentFile(file_path)
#    file1.Upload()
#    
#    os.remove(file_path)
