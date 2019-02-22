import Atrium2 as AC
import numpy as np
#import Atrium_Hex1 as AC
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def update1(frame_number, mat,A):
    """Next frame update for animation without ECG"""
    if frame_number in A.pace:
        A.SinusRhythm()
    A.Relaxing_ani()
    A.Conduct()
    data = A.phases.reshape([A.size,A.size])
    mat.set_data(data)
    return mat,


A = AC.Atrium()
np.random.rand(A.seed_prop)

fig1 = plt.figure(figsize = [15,15])
ax = fig1.subplots(1,1)
mat1 = ax.matshow(A.phases.reshape([A.size,A.size]),cmap=plt.cm.gray_r)

mat1.set_clim(0,A.rp)
ax.set_axis_off()
ani1 = animation.FuncAnimation(fig1, update1, frames = A.tot_time
                               ,fargs = (mat1,A), interval=100, repeat = None)
plt.axis([0,A.size,0,A.size])
