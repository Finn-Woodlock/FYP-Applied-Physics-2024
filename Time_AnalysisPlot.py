import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
import sys
from scipy.stats import linregress
Dir = 'output/AlloyOFF_8_2_Timed/VariedOut/'
Dir = 'output/LRIMOFFTimed/VariedOut/'
Dir = 'output/AlloyOFFTimed/VariedOut/'

bigfont = False
pres = False

saveFormat = Dir+"Figure_{}_{}_"+".png"

plt.rcParams.update({'font.size': 45})
plt.rcParams.update({'font.family': 'Times New Roman'})

from matplotlib import colormaps as cm

c1 = 'k'
c2 = 'grey'
cmap = cm.get_cmap('gist_grey')
if pres:
    plt.rcParams.update({'font.size': 30})
    plt.rcParams.update({'font.family': 'Calibri'})
    plt.rcParams.update({'axes.facecolor':(255/255,248/255,197/255,1)})
    plt.rcParams.update({'legend.facecolor':(195/255,179/255,105/255,255/255)})
    plt.rcParams.update({'legend.edgecolor':'none'})
    plt.rcParams.update({'axes.edgecolor':(0/255,34/255,78/255,255/255)})
    plt.rcParams.update({'axes.labelcolor':(0/255,34/255,78/255,255/255)})
    plt.rcParams.update({'axes.titlecolor':(0/255,34/255,78/255,255/255)})
    plt.rcParams.update({'text.color':(0/255,34/255,78/255,255/255)})
    plt.rcParams.update({'ytick.color':(0/255,34/255,78/255,255/255)})
    plt.rcParams.update({'xtick.color':(0/255,34/255,78/255,255/255)})
    plt.rcParams.update({'ytick.labelcolor':(0/255,34/255,78/255,255/255)})
    plt.rcParams.update({'xtick.labelcolor':(0/255,34/255,78/255,255/255)})
    plt.rcParams.update({'axes.linewidth':1.6})
    c1 = (0/255,34/255,78/255,255/255)
    c2 = (82/255,89/255,109/255,255/255)
    cmap = cm.get_cmap('cividis')
if bigfont:
    plt.rcParams.update({'font.family': 'Calibri'})
    plt.rcParams.update({'font.size': 60})
    saveFormat  =  Dir+"BIG_" +"Figure_{}_{}_"+".png"

nskipT = 5
nskipR = 0      

Dirs = ['BruteForce','BruteForceOrdered','Modified','ModifiedOrdered']
files = [[]]*4



for i,dirx in enumerate(Dirs):
    arr = []
    paths = os.listdir(Dir+dirx+'/')
    paths = [path for path in paths if '.csv' in path]
    paths.sort(key = lambda x : int(x[:x.index('_')]))
    for path in paths:
        # check if current path is a file
        if os.path.isfile(os.path.join(Dir+dirx+'/', path)) and '.csv' in path:
            arr.append(str(os.path.join(Dir+dirx+'/', path)))
    files[i] = arr

s = [int(x[x[:x.rfind('cpp')-1].rfind('/')+1:x.rfind('cpp')-1])**2 for x in files[0]]
DBs = [[pd.read_csv(x,header=[0]).rename(columns=lambda x: x.strip()) for x in arr] for arr in files]

arrs = np.asarray([[x[["Times" , "Delta_Times"]].to_numpy() for x in arr] for arr in DBs])

print(linregress(np.asarray(s[:]),arrs[0,:,0,0]))
print(linregress(np.asarray(s[:]),arrs[1,:,0,0]))
print(linregress(np.asarray(s[:]),arrs[2,:,0,0]))
print(linregress(np.asarray(s[:]),arrs[3,:,0,0]))
print(linregress(np.asarray(s[6:]),arrs[0,6:,0,0]))
print(linregress(np.asarray(s[6:]),arrs[1,6:,0,0]))
print(linregress(np.asarray(s[6:]),arrs[2,6:,0,0]))
print(linregress(np.asarray(s[6:]),arrs[3,6:,0,0]))
print(linregress(np.log2(np.asarray(s[:])),np.log10(arrs[0,:,0,0])))
print(linregress(np.log2(np.asarray(s[:])),np.log10(arrs[1,:,0,0])))
print(linregress(np.log2(np.asarray(s[:])),np.log10(arrs[2,:,0,0])))
print(linregress(np.log2(np.asarray(s[:])),np.log10(arrs[3,:,0,0])))
print(linregress(np.log2(np.asarray(s[6:])),np.log10(arrs[0,6:,0,0])))
print(linregress(np.log2(np.asarray(s[6:])),np.log10(arrs[1,6:,0,0])))
print(linregress(np.log2(np.asarray(s[6:])),np.log10(arrs[2,6:,0,0])))
print(linregress(np.log2(np.asarray(s[6:])),np.log10(arrs[3,6:,0,0])))

print([f'{100*x:.1f}%' for x in (1-np.mean(arrs[3,:,:,0],1)/np.mean(arrs[1,:,:,0],1))])
print([f'{100*x:.1f}%' for x in (1-np.mean(arrs[2,:,:,0],1)/np.mean(arrs[0,:,:,0],1))])

fig, ax = plt.subplots(figsize=(16, 9),dpi=int(1000/9), frameon=False)
DeltaY = (max(np.max(arrs[0,:,:,0]+arrs[0,:,:,1]),np.max(arrs[2,:,:,0]+arrs[2,:,:,1])) - min(np.min(arrs[0,:,:,0]-arrs[0,:,:,1]),np.min(arrs[2,:,:,0]-arrs[2,:,:,1])))/10
ax.set(xlim=(min(s),max(s)), ylim=(np.clip(min(np.min(arrs[0,:,:,0]-arrs[0,:,:,1]),np.min(arrs[2,:,:,0]-arrs[2,:,:,1])),0,10),max(np.max(arrs[0,:,:,0]+arrs[0,:,:,1]),np.max(arrs[2,:,:,0]+arrs[2,:,:,1]))))
ax.set_xlabel("N : System Size")
ax.set_ylabel("t : Time taken (s)")
#ax.set_xscale('log', base=2)
#ax.set_yscale('log', base=10)
ax.xaxis.set_major_formatter('{x:.0e}')
ax.plot(s,np.mean(arrs[0,:,:,0],1),linestyle='--',marker='x',c=c2,label = 'Default')
ax.fill_between(s,np.max(arrs[0,:,:,0]+arrs[0,:,:,1],1),np.min(arrs[0,:,:,0]-arrs[0,:,:,1],1),color=c2,linestyle='--',alpha=0.3)
ax.plot(s,np.mean(arrs[2,:,:,0],1),marker='o',c=c1,label = 'Novel Method')
ax.fill_between(s,np.max(arrs[2,:,:,0]+arrs[2,:,:,1],1),np.min(arrs[2,:,:,0]-arrs[2,:,:,1],1),color=c1,alpha=0.3)
fig.tight_layout()
fig.subplots_adjust(wspace=None, hspace=None)
ax.spines[['right', 'top']].set_visible(False)
if not bigfont: ax.legend(loc=2)
fig.savefig(Dir+"TimeComparison.png")


fig, ax = plt.subplots(figsize=(16, 9),dpi=int(1000/9), frameon=False)
DeltaY = (max(np.max(arrs[1,:,:,0]+arrs[1,:,:,1]),np.max(arrs[3,:,:,0]+arrs[3,:,:,1])) - min(np.min(arrs[1,:,:,0]-arrs[1,:,:,1]),np.min(arrs[3,:,:,0]-arrs[3,:,:,1])))/10
ax.set(xlim=(min(s),max(s)), ylim=(np.clip(min(np.min(arrs[1,:,:,0]-arrs[1,:,:,1]),np.min(arrs[3,:,:,0]-arrs[3,:,:,1])),0,10),max(np.max(arrs[1,:,:,0]+arrs[1,:,:,1]),np.max(arrs[3,:,:,0]+arrs[3,:,:,1]))))
ax.set_xlabel("N : Number of points")
ax.set_ylabel("t : Time taken (s)")
#ax.set_xscale('log', base=2)
#x.set_yscale('log', base=10)
ax.xaxis.set_major_formatter('{x:.0e}')
#ax.xaxis.set_ticks(np.power(10,np.arange(int(np.floor(np.log10(max(s))))-0.5, 0.5+int(np.floor(np.log10(max(s)))),0.5)))
ax.plot(s,np.mean(arrs[1,:,:,0],1),linestyle='--',marker='x',c=c2,label = 'Default')
ax.fill_between(s,np.max(arrs[1,:,:,0]+arrs[1,:,:,1],1),np.min(arrs[1,:,:,0]-arrs[1,:,:,1],1),color=c2,linestyle='--',alpha=0.3)
ax.plot(s,np.mean(arrs[3,:,:,0],1),marker='o',c=c1,label = 'Novel Method')
ax.fill_between(s,np.max(arrs[3,:,:,0]+arrs[3,:,:,1],1),np.min(arrs[3,:,:,0]-arrs[3,:,:,1],1),color=c1,alpha=0.3)
fig.tight_layout()
fig.subplots_adjust(wspace=None, hspace=None)
ax.spines[['right', 'top']].set_visible(False)
if not bigfont: ax.legend(loc=2)
fig.savefig(Dir+"Ordered_TimeComparison.png")


fig, ax = plt.subplots(figsize=(16, 9),dpi=int(1000/9), frameon=False)
DeltaY = (max(np.max(arrs[0,:,:,0]+arrs[0,:,:,1]),np.max(arrs[2,:,:,0]+arrs[2,:,:,1])) - min(np.min(arrs[0,:,:,0]-arrs[0,:,:,1]),np.min(arrs[2,:,:,0]-arrs[2,:,:,1])))/10
ax.set(xlim=(min(s),max(s)), ylim=(np.clip(min(np.min(arrs[0,:,:,0]-arrs[0,:,:,1]),np.min(arrs[2,:,:,0]-arrs[2,:,:,1])),1e-5,10),max(np.max(arrs[0,:,:,0]+arrs[0,:,:,1]),np.max(arrs[2,:,:,0]+arrs[2,:,:,1]))))
ax.set_xlabel("N : System Size")
ax.set_ylabel("t : Time taken (s)")
ax.set_xscale('log', base=2)
ax.set_yscale('log', base=10)
ax.plot(s,np.mean(arrs[0,:,:,0],1),linestyle='--',marker='x',c=c2,label = 'Default')
ax.fill_between(s,np.max(arrs[0,:,:,0]+arrs[0,:,:,1],1),np.min(arrs[0,:,:,0]-arrs[0,:,:,1],1),color=c2,linestyle='--',alpha=0.3)
ax.plot(s,np.mean(arrs[2,:,:,0],1),marker='o',c=c1,label = 'Novel Method')
ax.fill_between(s,np.max(arrs[2,:,:,0]+arrs[2,:,:,1],1),np.min(arrs[2,:,:,0]-arrs[2,:,:,1],1),color=c1,alpha=0.3)
fig.tight_layout()
fig.subplots_adjust(wspace=None, hspace=None)
ax.spines[['right', 'top']].set_visible(False)
if not bigfont: ax.legend(loc=2)
fig.savefig(Dir+"TimeComparison_ll.png")

fig, ax = plt.subplots(figsize=(16, 9),dpi=int(1000/9), frameon=False)
DeltaY = (max(np.max(arrs[1,:,:,0]+arrs[1,:,:,1]),np.max(arrs[3,:,:,0]+arrs[3,:,:,1])) - min(np.min(arrs[1,:,:,0]-arrs[1,:,:,1]),np.min(arrs[3,:,:,0]-arrs[3,:,:,1])))/10
ax.set(xlim=(min(s),max(s)), ylim=(np.clip(min(np.min(arrs[1,:,:,0]-arrs[1,:,:,1]),np.min(arrs[3,:,:,0]-arrs[3,:,:,1])),1e-5,10),max(np.max(arrs[1,:,:,0]+arrs[1,:,:,1]),np.max(arrs[3,:,:,0]+arrs[3,:,:,1]))))
ax.set_xlabel("N : Number of points")
ax.set_ylabel("t : Time taken (s)")
ax.set_xscale('log', base=2)
ax.set_yscale('log', base=10)
ax.yaxis.set_major_formatter('{x:.0e}')
ax.xaxis.set_major_formatter('{x:.0e}')
ax.plot(s,np.mean(arrs[1,:,:,0],1),linestyle='--',marker='x',c=c2,label = 'Default')
ax.fill_between(s,np.max(arrs[1,:,:,0]+arrs[1,:,:,1],1),np.min(arrs[1,:,:,0]-arrs[1,:,:,1],1),color=c2,linestyle='--',alpha=0.3)
ax.plot(s,np.mean(arrs[3,:,:,0],1),marker='o',c=c1,label = 'Novel Method')
ax.fill_between(s,np.max(arrs[3,:,:,0]+arrs[3,:,:,1],1),np.min(arrs[3,:,:,0]-arrs[3,:,:,1],1),color=c1,alpha=0.3)
fig.tight_layout()
fig.subplots_adjust(wspace=None, hspace=None)
ax.spines[['right', 'top']].set_visible(False)
if not bigfont: ax.legend(loc=2)
fig.savefig(Dir+"Ordered_TimeComparison_ll.png")
plt.show()