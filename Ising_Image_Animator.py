import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col
import matplotlib.animation as animation
import pandas as pd
import os
import sys

OverWrite = True
fps = 32
n_frames = 32
skipN = 0
mirrors = 2
CombineTemps = True

#home_path = "PyOut\\LRIM4_Unordered_Segmented\\"
home_path = "output\\AlloyOFF32_8_2\\"
home_path = "output\\LRIMOFF32\\"
home_path = "output\\AlloyOFF32\\"

plt.rcParams.update({'font.size': 40})
plt.rcParams.update({'font.family': 'Calibri'})

from matplotlib import colormaps as cm

pres = True

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

#AllowedTemps = [1.,2.,4.,6.]
AllowedTemps = []

dir_path = "States\\"
out_path = "Images\\"
format_out = "Images\\Image_Animation_{:.2f}.gif"#_Mirror.gif"
writergif = animation.PillowWriter(fps=fps) 

res = {}
existing = []

dir_path = os.path.join(home_path, dir_path)
out_path = os.path.join(home_path, out_path)
format_out = os.path.join(home_path, format_out)
for path in os.listdir(dir_path):
    # check if current path is a file
    if os.path.isfile(os.path.join(dir_path, path)):
        if path.split('_')[1] in res.keys():
            res[path.split('_')[1]].append(os.path.join(dir_path, path))
        else:
            res[path.split('_')[1]] = [os.path.join(dir_path, path)]
for path in os.listdir(out_path):
    if os.path.isfile(os.path.join(out_path, path)):
        existing.append(os.path.join(out_path, path))


for k in res.keys():
    res[k].sort(key=lambda x: int(x[len(dir_path):].split('_')[-1].split('.')[0]))


def getImageMap(resShape,data,info,tmat,lims):
    Z = np.zeros((resShape[1],resShape[0]),dtype=int)
    x, y, c = data.T
    tinv = np.linalg.inv(tmat)
    if not tinv.all() == tmat.all():
        for I,i in enumerate(np.linspace(lims[0][0],lims[0][1],resShape[0])):
            for J,j in enumerate(np.linspace(lims[1][0],lims[1][1],resShape[1])):
                v = np.matmul(tinv,np.asarray([[i],[j]]).astype(float))
                past = np.asarray([-v[0]//info[0],-v[1]//info[1]]).reshape(2)
                v[0] = v[0] % info[0]
                v[1] = v[1] % info[1]
                Ps = np.asarray([[0,0,0,1,1,1,1,1,2,2,2,2,2],[0,1,2,2,1,0,-1,-2,2,1,0,-1,-2]])
                v2idx = np.asarray([list(Ps[:,i]+(past)) for i in range(len(Ps[0,:]))]+[list(-Ps[:,i+1]+(past)) for i in range(len(Ps[0,:-1]))])
                v2 = np.r_[[v + (v2idx[i]*info).reshape(v.shape) for i in range(len(v2idx))]].reshape(len(v2idx),2).T
                v2 = np.matmul(tmat,v2)
                xidx = [np.argmin(np.square((x-v2[0,k]))+np.square((y-v2[1,k]))) for k in range(len(v2idx))]
                xsqs = [np.square((x[xidx[k]]-v2[0,k]))+np.square((y[xidx[k]]-v2[1,k])) for k in range(len(v2idx))]
                past = np.asarray(v2idx[np.argmin(xsqs)])
                xidx = xidx[np.argmin(xsqs)]
                Z[-J-1,I] = xidx
    else:
        for I,i in enumerate(np.linspace(lims[0][0],lims[0][1],resShape[0])):
            for J,j in enumerate(np.linspace(lims[1][0],lims[1][1],resShape[1])):
                v = np.asarray([[i],[j]])
                past = np.asarray([-v[0]//info[0],-v[1]//info[1]]).reshape(2)
                Ps = np.asarray([[0,0,0,1,1,1,1,1,2,2,2,2,2],[0,1,2,2,1,0,-1,-2,2,1,0,-1,-2]])
                v2idx = np.asarray([list(Ps[:,i]+(past)) for i in range(len(Ps[0,:]))]+[list(-Ps[:,i+1]+(past)) for i in range(len(Ps[0,:-1]))])
                v2 = np.r_[[v + (v2idx[i]*info).reshape(v.shape) for i in range(len(v2idx))]].reshape(len(v2idx),2).T
                xidx = [np.argmin(np.square((x-v2[0,k]))+np.square((y-v2[1,k]))) for k in range(len(v2idx))]
                xsqs = [np.square((x[xidx[k]]-v2[0,k]))+np.square((y[xidx[k]]-v2[1,k])) for k in range(len(v2idx))]
                xidx = xidx[np.argmin(xsqs)]
                Z[-J-1,I] = xidx
    return Z

class AnimatedScatter(object):

    def __init__(self,files,nRes=5,mirrors = 1,text=[]):
        self.files = files
        self.stream = self.data_stream()
        self.fstream = self.file_stream()
        self.mirrors = mirrors
        self.text = text
        self.idx = -1
        self.setup_plot(nRes)

        self.ani = animation.FuncAnimation(self.fig, self.update, interval=1000/fps, init_func=self.init_plot, blit=True,frames=len(files))

    def setup_plot(self,nRes):
        self.init_data = next(self.stream)
        data, info,tmat = self.init_data
        x, y, c = data.T

        tinv = np.linalg.inv(tmat)
        rawPoints = np.matmul(tinv,data[:,:2].T)
        rx, ry = rawPoints
        rinfo = np.matmul(tinv,info[-2:].T)
        minrx,maxrx = 0, rinfo[0]
        minry,maxry = 0, rinfo[1]

        drat = (maxrx-minrx)/(maxry-ry.min())

        minx,maxx,diffx,meanx,adjx = 0, info[-2], info[-2], info[-2]/2, (maxrx-minrx)/np.sqrt(4*drat*len(x))
        miny,maxy,diffy,meany,adjy = 0, info[-1], info[-1], info[-1]/2, (maxry-minry)/np.sqrt((4/drat)*len(y))
        

        adjrx, adjry= (maxrx-minrx)/np.sqrt(4*drat*len(rx)), (maxry-minry)/np.sqrt((4/drat)*len(ry))
        adjx,adjy = np.matmul(tmat,[[adjrx],[0]]).flatten()[0], np.matmul(tmat,[[0],[adjry]]).flatten()[1]
        
        minrx,maxrx = minrx-adjrx, maxrx+adjrx
        minry,maxry = minry-adjry, maxry+adjry

        self.coord = [list(np.matmul(tmat,[[minrx],[minry]]).flatten()), list(np.matmul(tmat,[[maxrx+adjrx],[minry]]).flatten()),list(np.matmul(tmat,[[maxrx],[maxry]]).flatten()), list(np.matmul(tmat,[[minrx],[maxry]]).flatten())]
        self.coord.append(self.coord[0])

        self.lims = ((minx-2*adjx-((self.mirrors-1)*diffx)/2, maxx+2*adjx+((self.mirrors-1)*diffx)/2),(miny-2*adjy-((self.mirrors-1)*diffy)/2, maxy+2*adjy+((self.mirrors-1)*diffy)/2))

        self.res = (nRes*np.rint((self.mirrors*diffx+4*adjx)/adjx).astype(int),nRes*np.rint((self.mirrors*diffx+4*adjx)/adjy).astype(int))

        self.imMap = getImageMap(self.res,data.astype(float),info[3:],tmat,self.lims)

        self.fig, self.ax = plt.subplots(figsize=(9*(np.rint(diffx+4*adjx)/np.rint(diffy+4*adjy)), 9), dpi=self.res[1]/9, frameon=False)

        self.im = self.ax.imshow(data[self.imMap,2],cmap=cmap,aspect='auto', norm = col.Normalize(vmin=self.spinList.min(), vmax=self.spinList.max()),extent=[self.lims[0][0], self.lims[0][1], self.lims[1][0], self.lims[1][1]])
        


        self.ax.set(xlim=self.lims[0], ylim=self.lims[1])
        self.ax.set_aspect('equal', 'box')
        self.ax.set_yticklabels([])
        self.ax.set_xticklabels([])
        self.ax.set_xticks([])
        self.ax.set_yticks([]) 
        if len(self.text) > 0: 
            self.I,T = False,""
            for i,t in self.text:
                if self.idx >= i:
                    T = t
                    self.I = True
            if self.I:
                self.time_text = self.ax.text(self.lims[0][1]-(self.lims[0][1]-self.lims[0][0])*0.02, self.lims[1][1]-(self.lims[1][1]-self.lims[1][0])*0.02,"T = "+T,horizontalalignment='right', verticalalignment='top', color=c1,bbox=dict(facecolor=(255/255,248/255,197/255,1), alpha=0.5))
        self.fig.tight_layout()
        self.fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)


    def init_plot(self):
        data, info,tmat = self.init_data

        self.im.set_data(data[self.imMap,2])

        if len(self.text) > 0: 
            if not self.I:
                T = ""
                for i,t in self.text:
                    if self.idx >= i:
                        T = t
                        self.I = True
                if self.I:
                    self.time_text = self.ax.text(self.lims[0][1]-(self.lims[0][1]-self.lims[0][0])*0.02, self.lims[1][1]-(self.lims[1][1]-self.lims[1][0])*0.02,"T = "+T,horizontalalignment='right', verticalalignment='top', color=c1,bbox=dict(facecolor=(255/255,248/255,197/255,1), alpha=0.5))
            else:
                T = ""
                for i,t in self.text:
                    if self.idx >= i:
                        T = t
                self.time_text.set_text("T = "+T)
                return self.im, self.time_text,

        return self.im,

    def data_stream(self):
        while True:
            f=next(self.fstream)
            ndim = pd.read_csv(f, skiprows = lambda x: x not in [0],header=None).to_numpy()[0,1]
            df2 = pd.read_csv(f, skiprows = lambda x: x not in ([1,2,3]+[4+x for x in range(ndim)]),header=None)
            df = pd.read_csv(f, skiprows=[0,1,2,3]+[4+x for x in range(2*ndim+1)],header=0).rename(columns=lambda x: x.strip())
            tmat = pd.read_csv(f, skiprows=lambda x: x not in [4+ndim+x for x in range(ndim)],header=None).to_numpy().astype(float)
            self.spinList = pd.read_csv(f, skiprows=lambda x: x not in [4+2*ndim],header=None).to_numpy().astype(float)
            c = df[["Spin"]].to_numpy()
            xy = df[["X1","X2"]].to_numpy()
            yield (np.c_[xy[:,0], xy[:,1], c] , df2.to_numpy()[:,1].T,tmat)

    def file_stream(self):
        idx = -1
        while True:
            self.idx += 1
            idx += 1
            idx = idx % len(self.files)
            yield self.files[idx]

    def update(self, i):
        data,info,tmat = next(self.stream)

        self.im.set_data(data[self.imMap,2])

        if len(self.text) > 0: 
            if not self.I:
                T = ""
                for i,t in self.text:
                    if self.idx >= i:
                        T = t
                        self.I = True
                if self.I:
                    self.time_text = self.ax.text(self.lims[0][1]-(self.lims[0][1]-self.lims[0][0])*0.02, self.lims[1][1]-(self.lims[1][1]-self.lims[1][0])*0.02,"T = "+T,horizontalalignment='right', verticalalignment='top', color=c1,bbox=dict(facecolor=(255/255,248/255,197/255,1), alpha=0.5))
            else:
                T = ""
                for i,t in self.text:
                    if self.idx >= i:
                        T = t
                self.time_text.set_text("T = "+T)
                return self.im, self.time_text,
    
        return self.im,



if __name__ == '__main__':
    Ts = list(res.keys())
    Ts.sort(key=lambda x : float(x))
    F = []
    TT = []
    if CombineTemps:
        for i,t in enumerate(Ts[::-1]):
            if skipN > 0 and i % skipN != 0:
                continue
            if len(res[t]) < n_frames and n_frames > 0:
                continue
            if not (len(AllowedTemps) == 0 or float(t) in AllowedTemps):
                continue
            TT.append((len(F),f"{float(t):.2f}"))
            if n_frames > 0: F += (res[t][:n_frames])
            else: F += (res[t])
        if not format_out.format(float(-1)) in existing or OverWrite:
            print([(f"{TT[i-1][0]}-{TT[i][0]}",TT[i-1][1]) for i in range(1,len(TT))]+[(f"{TT[-1][0]}-{len(F)}",TT[-1][1])])
            a = AnimatedScatter(F,5,mirrors,TT)
            print("Saving : " + format_out.format(float(-1)))
            a.ani.save(format_out.format(float(-1)), writer=writergif)
            print("Saved  : " + format_out.format(float(-1)))
    else:
        for i,t in enumerate(Ts):
            if skipN > 0 and i % skipN != 0:
                continue
            if len(res[t]) < n_frames and n_frames > 0:
                continue
            if not (len(AllowedTemps) == 0 or float(t) in AllowedTemps):
                continue
            if not format_out.format(float(t)) in existing or OverWrite:
                a = AnimatedScatter(res[t][:n_frames if n_frames > 0 else len(res[t])],5,mirrors)
                print("Saving : " + format_out.format(float(t)))
                a.ani.save(format_out.format(float(t)), writer=writergif)
                print("Saved  : " + format_out.format(float(t)))
            else:
                continue
