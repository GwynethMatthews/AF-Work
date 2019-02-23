import numpy as np
import Atrium_Hex1 as AC
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib as matplotlib
from matplotlib import collections, transforms

def update1(frame_number,collection,A):
    """Next frame update for animation without ECG"""
    if frame_number in A.pace:
        A.SinusRhythm()
    A.Relaxing_ani()
    A.Conduct()
    collection.set_array(np.array(A.phases))
    #ax.add_collection(collection, autolim=True)
    #ax.set_xlim(0,A.size) 
    #ax.set_ylim(0,A.size) 
    
    return ax,


A = AC.Atrium()
np.random.rand(A.seed_prop)

fig1 = plt.figure(figsize = [15,15])
ax = fig1.subplots(1,1)
patches = []
offsets = []
for i in A.x[0]:
    for j in A.y[:,0]:
        if i%2 ==0:
            offsets.extend([(j+0.5,i)]) 
        else:
            offsets.extend([(j+1,i)]) 
for k in offsets:
    l = matplotlib.patches.RegularPolygon(k,6,radius = 0.55)
    patches.extend([l])
collection = collections.PatchCollection(patches, cmap=plt.cm.gray_r)
ax.add_collection(collection, autolim=True)

ax.set_axis_off()
ani1 = animation.FuncAnimation(fig1, update1, frames = A.tot_time
                               ,fargs = (collection,A), interval=100, repeat = None)
plt.axis([0,A.size+1,0,A.size+1])
