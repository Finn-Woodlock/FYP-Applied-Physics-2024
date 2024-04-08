import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import labellines

Dir = 'FinalData/Gen6/AlloyOFF16_8_2/'
Dir = 'FinalData/Gen6/AlloyOFF16/'
Dir = 'FinalData/Gen6/LRIMOFF16/'

bigfont = False
autocorrect = False
inverse = False
pres = True

saveFormat = Dir+"Figure_{}_{}_"+f"{inverse}_{autocorrect}"+".png"

plt.rcParams.update({'font.size': 45})
plt.rcParams.update({'font.family': 'Times New Roman'})
from matplotlib import colormaps as cm
c1 = 'k'
c2 = 'grey'
cmap = cm.get_cmap('seismic')
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
    saveFormat  =  Dir+"BIG_" +"Figure_{}_{}_"+f"{inverse}_{autocorrect}"+".png"
nskipT = 5
nskipR = 0      

try:
    number = Dir[Dir.index('OFF')+3:].split('_')[0].split('/')[0]       
    float(number)
except:
    try:
        number = Dir[Dir.index('ON')+3:].split('_')[0].split('/')[0]    
        float(number)
    except:
        number = ""              
            
try:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
    Data = number+'_cppOut.csv'

    Data = pd.read_csv(Dir+Data,header=0).rename(columns=lambda x: x.strip())
    if inverse: Data['T'] = 1/Data['T']
    if nskipT > 0: Data=Data.drop(index=[len(Data)-1 -i for i in range(nskipT)])
    if autocorrect:
        try:     
            Data['EErr'] = Data['EErr'] * np.sqrt((1+Data['EAC'])/(1-Data['EAC']))
            Data['C'] = Data['C'] * ((1+Data['EAC'])/(1-Data['EAC']))
            Data['CErr'] = Data['CErr'] * ((1+Data['CAC'])/(1-Data['CAC']))* ((1+Data['EAC'])/(1-Data['EAC']))
            try:     
                Data['MErr'] = Data['MErr'] * np.sqrt((1+Data['MAC'])/(1-Data['MAC']))
                Data['|M|Err'] = Data['|M|Err'] * np.sqrt((1+Data['|M|AC'])/(1-Data['|M|AC']))
                Data['X'] = Data['X'] * ((1+Data['MAC'])/(1-Data['MAC']))
                Data['X(|M|)'] = Data['X(|M|)'] * ((1+Data['|M|AC'])/(1-Data['|M|AC']))
                Data['XErr'] = Data['XErr'] *  np.sqrt((1+Data['XAC'])/(1-Data['XAC']))* ((1+Data['MAC'])/(1-Data['MAC']))
                Data['X(|M|)Err'] = Data['X(|M|)Err'] *  np.sqrt((1+Data['X(|M|)AC'])/(1-Data['X(|M|)AC'])) * ((1+Data['|M|AC'])/(1-Data['|M|AC']))
            except: pass
        except: pass
    try:
        fig, ax = plt.subplots(figsize=(16, 9),dpi=int(1000/9), frameon=False)
        X, Y, err = Data['T'], Data['E'], Data['EErr']
        DeltaY = (max(Y+err)-min(Y+err))/10
        ax.set(xlim=(min(X),max(X)), ylim=(min(Y-err)-DeltaY,max(Y+err)+DeltaY))
        ax.set_xlabel("T : Temperature")
        ax.set_ylabel("E : Internal Energy")
        ax.plot(X,Y,marker='o',c=c1)
        ax.fill_between(X,Y+err,Y-err,color=c1,alpha=0.3)
        fig.tight_layout()
        fig.subplots_adjust(wspace=None, hspace=None)
        ax.spines[['right', 'top']].set_visible(False)
        fig.savefig(saveFormat.format("Energy","Plot"))

        fig, ax = plt.subplots(figsize=(16, 9),dpi=int(1000/9), frameon=False)
        X, Y, err = Data['T'], Data['C'], Data['CErr']
        DeltaY = (max(Y+err)-min(Y+err))/10
        ax.set(xlim=(min(X),max(X)), ylim=(min(Y-err)-DeltaY,max(Y+err)+DeltaY))
        ax.set_xlabel("T : Temperature")
        ax.set_ylabel("C : Heat Capacity")
        ax.plot(X,Y,marker='o',c=c1)
        ax.fill_between(X,Y+err,Y-err,color=c1,alpha=0.3)
        fig.tight_layout()
        fig.subplots_adjust(wspace=None, hspace=None)
        ax.spines[['right', 'top']].set_visible(False)
        fig.savefig(saveFormat.format("Heat_Capacity","Plot"))
        try:
            fig, ax = plt.subplots(figsize=(16, 9),dpi=int(1000/9), frameon=False)
            X, Y, err = Data['T'], Data['M'], Data['MErr']
            ax.set(xlim=(min(X),max(X)), ylim=(-1,1))
            ax.set_xlabel("T : Temperature")
            ax.set_ylabel("M : Signed Magnetisation")
            ax.plot(X,Y,marker='o',c=c1)
            ax.fill_between(X,Y+err,Y-err,color=c1,alpha=0.3)
            fig.tight_layout()
            fig.subplots_adjust(wspace=None, hspace=None)
            ax.spines[['right', 'top']].set_visible(False)
            fig.savefig(saveFormat.format("S_Magnetism","Plot"))
            
            fig, ax = plt.subplots(figsize=(16, 9),dpi=int(1000/9), frameon=False)
            X, Y, err = Data['T'], Data['X'], Data['XErr']
            DeltaY = (max(Y+err)-min(Y+err))/10
            ax.set(xlim=(min(X),max(X)), ylim=(min(Y-err)-DeltaY,max(Y+err)+DeltaY))
            ax.set_xlabel("T : Temperature")
            ax.set_ylabel("χ : Magnetic Susceptibility")
            ax.plot(X,Y,marker='o',c=c1)
            ax.fill_between(X,Y+err,Y-err,color=c1,alpha=0.3)
            fig.tight_layout()
            fig.subplots_adjust(wspace=None, hspace=None)
            ax.spines[['right', 'top']].set_visible(False)
            fig.savefig(saveFormat.format("S_Susceptibility","Plot"))

            fig, ax = plt.subplots(figsize=(16, 9),dpi=int(1000/9), frameon=False)
            X, Y, err = Data['T'], Data['|M|'], Data['|M|Err']
            ax.set(xlim=(min(X),max(X)), ylim=(0,1))
            ax.set_xlabel("T : Temperature")
            ax.set_ylabel("|M| : Unsigned Magnetisation")
            ax.plot(X,Y,marker='o',c=c1)
            ax.fill_between(X,Y+err,Y-err,color=c1,alpha=0.3)
            fig.tight_layout()
            fig.subplots_adjust(wspace=None, hspace=None)
            ax.spines[['right', 'top']].set_visible(False)
            fig.savefig(saveFormat.format("US_Magnetism","Plot"))
            
            fig, ax = plt.subplots(figsize=(16, 9),dpi=int(1000/9), frameon=False)
            X, Y, err = Data['T'], Data['X(|M|)'], Data['X(|M|)Err']
            DeltaY = (max(Y+err)-min(Y+err))/10
            ax.set(xlim=(min(X),max(X)), ylim=(min(Y-err)-DeltaY,max(Y+err)+DeltaY))
            ax.set_xlabel("T : Temperature")
            ax.set_ylabel("χ : Magnetic Susceptibility")
            ax.plot(X,Y,marker='o',c=c1)
            ax.fill_between(X,Y+err,Y-err,color=c1,alpha=0.3)
            fig.tight_layout()
            fig.subplots_adjust(wspace=None, hspace=None)
            ax.spines[['right', 'top']].set_visible(False)
            fig.savefig(saveFormat.format("US_Susceptibility","Plot"))
    
        except:
            pass
    except:
        pass

    try:
        fig, ax = plt.subplots(figsize=(16, 9),dpi=int(1000/9), frameon=False)
        ax.set_xlabel("T : Temperature")
        ax.set_ylabel("R : Serial Correlation")
        X, Y = Data['T'], Data['EAC']
        DeltaY = (max(Y)-min(Y))/10
        Xlims, Ylims = (min(X),max(X)) , (min(Y)-DeltaY,max(Y)+DeltaY)
        ax.plot(X,Y,marker='o',label = 'E : Internal Energy', c = c1, linestyle = '-')
        X, Y = Data['T'], Data['CAC']
        DeltaY = (max(Y)-min(Y))/10
        Xlims, Ylims = (min(min(X),Xlims[0]),max(max(X),Xlims[1])) , (min(min(Y)-DeltaY,Ylims[0]),max(max(Y)+DeltaY,Ylims[1]))
        ax.plot(X,Y,marker='x',label = 'C : Heat Capacity', c = c2, linestyle = '-')
        fig.tight_layout()
        fig.subplots_adjust(wspace=None, hspace=None)
        ax.spines[['right', 'top']].set_visible(False)
        ax.set(xlim=Xlims, ylim=Ylims)
        fig.legend(loc=1)
        
        fig.savefig(saveFormat.format("EnergyCapatity","AC"))
        try:
            fig, ax = plt.subplots(figsize=(16, 9),dpi=int(1000/9), frameon=False)
            ax.set_xlabel("T : Temperature")
            ax.set_ylabel("R : Serial Correlation")
            X, Y = Data['T'], Data['MAC']
            DeltaY = (max(Y)-min(Y))/10
            Xlims, Ylims = (min(X),max(X)) , (min(Y)-DeltaY,max(Y)+DeltaY)
            ax.plot(X,Y,marker='o',label = 'M : Signed Magnetisation', c = c1, linestyle = '-')
            X, Y = Data['T'], Data['XAC']
            DeltaY = (max(Y)-min(Y))/10
            Xlims, Ylims = (min(min(X),Xlims[0]),max(max(X),Xlims[1])) , (min(min(Y)-DeltaY,Ylims[0]),max(max(Y)+DeltaY,Ylims[1]))
            ax.plot(X,Y,marker='x',label = 'χ : Magnetic Susceptibility', c = c2, linestyle = '-')
            fig.tight_layout()
            fig.subplots_adjust(wspace=None, hspace=None)
            ax.spines[['right', 'top']].set_visible(False)
            ax.set(xlim=Xlims, ylim=Ylims)
            
            fig.savefig(saveFormat.format("S_MagSucc","AC"))

            fig, ax = plt.subplots(figsize=(16, 9),dpi=int(1000/9), frameon=False)
            ax.set_xlabel("T : Temperature")
            ax.set_ylabel("R : Serial Correlation")
            X, Y = Data['T'], Data['|M|AC']
            DeltaY = (max(Y)-min(Y))/10
            Xlims, Ylims = (min(X),max(X)) , (min(Y)-DeltaY,max(Y)+DeltaY)
            ax.plot(X,Y,marker='o',label = '|M| : Unsigned Magnetisation', c = c1, linestyle = '-')
            X, Y = Data['T'], Data['X(|M|)AC']
            DeltaY = (max(Y)-min(Y))/10
            Xlims, Ylims = (min(min(X),Xlims[0]),max(max(X),Xlims[1])) , (min(min(Y)-DeltaY,Ylims[0]),max(max(Y)+DeltaY,Ylims[1]))
            ax.plot(X,Y,marker='x',label = 'χ : Magnetic Susceptibility', c = c2, linestyle = '-')
            fig.tight_layout()
            fig.subplots_adjust(wspace=None, hspace=None)
            ax.spines[['right', 'top']].set_visible(False)
            ax.set(xlim=Xlims, ylim=Ylims)
            fig.legend(loc=1)
            
            fig.savefig(saveFormat.format("US_MagSucc","AC"))
            pass
        except:
            pass
    except:
        pass
except: pass

try:
    vals = np.ones((256, 4))
    vals[:, 0] = np.linspace(0/256, 73/256, 256)
    vals[:, 1] = np.linspace(39/256, 211/256, 256)
    vals[:, 2] = np.linspace(93/256, 102/256, 256)

    newcmap = ListedColormap(vals)
    Data = number+'_cppOutG.csv'

    Data = pd.read_csv(Dir+Data,header=0).rename(columns=lambda x: x.strip())
    Data = Data.astype(float)
    Cols = Data.columns[1:]
    R, Data = Data.to_numpy()[nskipR:,0], Data.to_numpy()[nskipR:,1:]
    if (Cols[0]+"Err" in Cols):
        if (Cols[0]+"AC" in Cols):
            Data,DataErr,DataAC = Data[:,:Data.shape[1]//3], Data[:,Data.shape[1]//3:2*Data.shape[1]//3], Data[:,2*Data.shape[1]//3:]
        else : Data,DataErr = Data[:,:Data.shape[1]//2], Data[:,Data.shape[1]//2:]
    if nskipT > 0:
        Data = Data[:,:-nskipT]
        if (Cols[0]+"Err" in Cols):
            DataErr = DataErr[:,:-nskipT]
            if (Cols[0]+"AC" in Cols):
                DataAC = DataAC[:,:-nskipT]
        Cols = Cols[:-nskipT]
    filter = np.logical_not(np.isnan(Data[:,0]))
    R = R[filter].astype(float)
    Data = Data[filter].astype(float)
    if (Cols[0]+"Err" in Cols):
        DataErr = np.abs(DataErr[filter].astype(float))
        if (Cols[0]+"AC" in Cols):
            DataAC = DataAC[filter].astype(float)
            
    if autocorrect:
        try:     
            DataErr = DataErr * np.sqrt((1+DataAC)/(1-DataAC))
        except: pass

    Cols = np.asarray([s[:s.rfind('0')] for s in Cols]).astype(float)[:Data.shape[1]]
    if inverse: Cols = 1/Cols
    Rx,Ry = np.meshgrid(Cols,R)

    fig, ax = plt.subplots(figsize=(16, 9),dpi=int(1000/9), frameon=False)
    X, Y = Rx, Ry
    ax.set(xlim=(np.min(X),np.max(X)), ylim=(np.min(Y),np.max(Y)))
    ax.set_xlabel("T : Temperature")
    ax.set_ylabel("r : Distance")
    CS = ax.contourf(Rx,Ry,Data,cmap=cmap,levels=100,vmin=-1.,vmax=1.)
    #CS = ax.contourf(Rx,Ry,Data,cmap=newcmap,levels=100,vmin=np.min(Data),vmax=np.max(Data))
    ax.contour(Rx,Ry,Data,levels=np.arange(-0.8,0.8,0.2),colors=[c1],vmin=-1.,vmax=1.,linestyles='dashed',linewidths=0.5)
    ax.contour(Rx,Ry,Data,levels=[0],colors=[c1],vmin=-1.,vmax=1.,linestyles='solid',linewidths=0.7)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(CS, cax=cax)
    cbar.set_label('g(r) : Spin-Spin Correlation', rotation=90)
    fig.tight_layout()
    fig.subplots_adjust(wspace=None, hspace=None)
    fig.savefig(saveFormat.format("Correlation","Plot"))

    fig, ax = plt.subplots(figsize=(16, 9),dpi=int(1000/9), frameon=False)
    X, Y = Rx, Ry
    ax.set(xlim=(np.min(X),np.max(X)), ylim=(np.min(Y),np.max(Y)))
    ax.set_xlabel("T : Temperature")
    ax.set_ylabel("r : Distance")
    CS = ax.contourf(Rx,Ry,DataErr,cmap=cmap,levels=100,vmin=-np.max(DataErr),vmax=np.max(DataErr))
    #CS = ax.contourf(Rx,Ry,DataErr,cmap=newcmap,levels=100,vmin=np.min(DataErr),vmax=np.max(DataErr))
    ax.contour(Rx,Ry,DataErr,levels=np.arange(-0.8,0.8,0.2),colors=[c1],vmin=-np.max(DataErr),vmax=np.max(DataErr),linestyles='dashed',linewidths=0.5)
    ax.contour(Rx,Ry,DataErr,levels=[0],colors=[c1],vmin=-np.max(DataErr),vmax=np.max(DataErr),linestyles='solid',linewidths=0.7)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(CS, cax=cax)
    cbar.set_label('σ : Standard Error', rotation=90)
    fig.tight_layout()
    fig.subplots_adjust(wspace=None, hspace=None)
    fig.savefig(saveFormat.format("Correlation","Err"))

    fig, ax = plt.subplots(figsize=(16, 9),dpi=int(1000/9), frameon=False)
    X, Y = Rx, Ry
    ax.set(xlim=(np.min(X),np.max(X)), ylim=(np.min(Y),np.max(Y)))
    ax.set_xlabel("T : Temperature")
    ax.set_ylabel("r : Distance")
    CS = ax.contourf(Rx,Ry,DataAC,cmap=cmap,levels=100,vmin=-1.,vmax=1.)
    #CS = ax.contourf(Rx,Ry,DataAC,cmap=newcmap,levels=100,vmin=np.min(DataAC),vmax=np.max(DataAC))
    ax.contour(Rx,Ry,DataAC,levels=np.arange(-0.8,0.8,0.2),colors=[c1],vmin=-1.,vmax=1.,linestyles='dashed',linewidths=0.5)
    ax.contour(Rx,Ry,DataAC,levels=[0],colors=[c1],vmin=-1.,vmax=1.,linestyles='solid',linewidths=0.7)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(CS, cax=cax)
    cbar.set_label('R : Serial Correlation', rotation=90)
    fig.tight_layout()
    fig.subplots_adjust(wspace=None, hspace=None)
    fig.savefig(saveFormat.format("Correlation","AC"))
except: pass
