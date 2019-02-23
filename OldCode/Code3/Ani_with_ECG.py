import Atrium_class as AC
import numpy as np
import matplotlib.gridspec as mgs
import matplotlib.pyplot as plt
import matplotlib.animation as animation
# Animation

def update2(frame_number, mat,A,denominator):
    """Next frame update for animation with ECG"""
    if frame_number in A.pace:
        A.SinusRhythm()
    A.Relaxing_ani()
    A.Conduct()
    data = A.phases.reshape([A.size,A.size])
    mat.set_data(data)
    ECG_value = A.ECG(denominator)
    A.potentials = np.roll(A.potentials,-1)
    A.potentials[-1] = ECG_value
    ax2.clear()
    ax2.plot(A.time_for_ECG,A.potentials)
    ax2.set_xlim(-500,0) 
    #ax2.set_ylim(-400,400)
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Potential')
    return mat,


A = AC.Atrium(v = 0.2)
LoP = [100.5,100.5]
np.random.seed(A.seed_prop)
fig2 = plt.figure()
mgs.GridSpec(4,5)
ax1 = plt.subplot2grid((4,5), (0,1), colspan=3, rowspan=3)
ax2 = plt.subplot2grid((4,5), (3,0), colspan=5)
#ax1.set_axis_off()
mat2 = ax1.matshow(A.phases.reshape([A.size,A.size]),cmap=plt.cm.gray_r)
mat2.set_clim(0,A.rp)
ani2 = animation.FuncAnimation(fig2, update2, frames = A.tot_time
                    ,fargs = (mat2,A,LoP), interval=100, repeat = None)

plt.show()