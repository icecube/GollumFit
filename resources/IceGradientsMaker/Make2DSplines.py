import numpy as np

import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import matplotlib.colors as colors

from matplotlib import rcParams
from matplotlib.ticker import MultipleLocator,LogLocator

import photospline

import sys

rcParams['font.serif']            = 'Times new Roman'
rcParams['font.weight']           = 200
rcParams['font.size']             = 30
rcParams['text.usetex']           = True

rcParams['grid.color']            = 'black'
rcParams['grid.alpha']            = 0.10
rcParams['grid.linestyle']        = '-'

rcParams['axes.grid']             = False
rcParams['axes.linewidth']        = 1.5
rcParams['axes.labelpad']         = 14.0
rcParams['axes.labelsize']        = 30.0

rcParams['xtick.labelsize']       = 30
rcParams['ytick.labelsize']       = 30
rcParams['xtick.direction']       = 'in'
rcParams['ytick.direction']       = 'in'
rcParams['xtick.major.pad']       = 10
rcParams['ytick.major.pad']       = 10
rcParams['xtick.major.width']     = 1.5
rcParams['ytick.major.width']     = 1.5
rcParams['xtick.top']             = True
rcParams['ytick.right']           = True
rcParams['xtick.minor.visible']   = True
rcParams['ytick.minor.visible']   = True

rcParams['legend.frameon']        = False
rcParams['legend.fontsize']       = 30

def Spline(GradType,StartMode,EndMode,SmoothVal) :

    TwoGrads=np.loadtxt(ingrad+"/"+topology+"_"+GradType+"Grads.csv")

    fig, axs = plt.subplots(4, 5, sharex='col', sharey='row',figsize=(25,20))
    flatax=np.ndarray.flatten(axs)

    nKnots  = 50
    order   = [2,2]
    smooth  = [SmoothVal,SmoothVal]
    penalty = [2,2]

    for mode in range(StartMode,EndMode):

        print(GradType,"->",mode)

        mask=TwoGrads[:,0]==mode

        ZenithLwEdge = np.unique(TwoGrads[mask,1])
        ZenithUpEdge = np.unique(TwoGrads[mask,2])

        EnergyLwEdge = np.log10(np.unique(TwoGrads[mask,3]))
        EnergyUpEdge = np.log10(np.unique(TwoGrads[mask,4]))

        ZenithCenter = (ZenithLwEdge + ZenithUpEdge)/2
        EnergyCenter = (EnergyLwEdge + EnergyUpEdge)/2

        hval = np.zeros([len(ZenithCenter),len(EnergyCenter)])
        herr = np.zeros([len(ZenithCenter),len(EnergyCenter)])

        counter = 0
        for i in range(len(ZenithCenter)) :
            for j in range(len(EnergyCenter)) :
                hval[i,j] = TwoGrads[mask,5][counter]
                herr[i,j] = 1./np.abs(TwoGrads[mask,6][counter])
                counter += 1

        ZenithWidth = np.abs(ZenithLwEdge[0]-ZenithUpEdge[0])
        EnergyWidth = np.abs(EnergyLwEdge[0]-EnergyUpEdge[0])

        ZenithKnots = np.linspace(ZenithLwEdge[0]-ZenithWidth,ZenithUpEdge[-1]+ZenithWidth,nKnots)
        EnergyKnots = np.linspace(EnergyLwEdge[0]-EnergyWidth,EnergyUpEdge[-1]+EnergyWidth,nKnots)

        data, err = photospline.ndsparse.from_data(hval, herr)

        spline = photospline.glam_fit(data, err, [ZenithCenter,EnergyCenter], [ZenithKnots,EnergyKnots], order, smooth, penalty, verbose=False)


        ZenithKnotsCenter = (ZenithKnots[1:]+ZenithKnots[:-1])/2
        EnergyKnotsCenter = (EnergyKnots[1:]+EnergyKnots[:-1])/2

        hfit = np.zeros([len(ZenithKnotsCenter),len(EnergyKnotsCenter)])
        for i in range(len(ZenithKnotsCenter)) :
            for j in range(len(EnergyKnotsCenter)) :
                hfit[i,j] = spline.evaluate_simple([ZenithKnotsCenter[i],EnergyKnotsCenter[j]])


        img1 = flatax[mode].pcolormesh(np.append(ZenithLwEdge,ZenithUpEdge[-1]), np.power(10,np.append(EnergyLwEdge,EnergyUpEdge[-1])), hval.transpose(), norm=colors.Normalize(vmin=-0.05, vmax=0.05), cmap=plt.get_cmap('RdBu_r'))
        img2 = flatax[mode+10].pcolormesh(ZenithKnots, np.power(10,EnergyKnots), hfit.transpose(), norm=colors.Normalize(vmin=-0.05, vmax=0.05), cmap=plt.get_cmap('RdBu_r'))
        flatax[mode].set_xlim(ZenithLwEdge[0],ZenithUpEdge[-1])
        flatax[mode+10].set_xlim(ZenithLwEdge[0],ZenithUpEdge[-1])
        flatax[mode].set_ylim(np.power(10,EnergyLwEdge[0]),np.power(10,EnergyUpEdge[-1]))
        flatax[mode+10].set_ylim(np.power(10,EnergyLwEdge[0]),np.power(10,EnergyUpEdge[-1]))
        flatax[mode].set_yscale('log')
        flatax[mode+10].set_yscale('log')
        flatax[mode].set_xlabel(r'$\cos(\theta_{z}^{reco})$')
        flatax[mode+10].set_xlabel(r'$\cos(\theta_{z}^{reco})$')
        if mode==0 :
            flatax[mode].set_ylabel(r'$\textrm{E}_{dnn}/\textrm{GeV}$')
            flatax[mode+10].set_ylabel(r'$\textrm{E}_{dnn}/\textrm{GeV}$')
        
        if mode==EndMode-1 :
            cbar = fig.colorbar(img1, ax=flatax[mode], ticks=np.linspace(-0.05, 0.05, 5))
            cbar.ax.set_ylabel(r"$\textrm{Fractional Shift, 1 }\sigma$")
            cbar = fig.colorbar(img2, ax=flatax[mode+10], ticks=np.linspace(-0.05, 0.05, 5))
            cbar.ax.set_ylabel(r"$\textrm{Fractional Shift, 1 }\sigma$")

        xi, yi = np.meshgrid(ZenithCenter, EnergyCenter)
        flatax[mode].scatter(xi, np.power(10,yi), color="red", marker='.', s=1)

        spline.write(outspl+"/"+topology+"_"+GradType+"_"+str(mode)+".fits");


    plt.tight_layout()
    plt.savefig(outplots+"/"+topology+"_"+GradType+"Gradients"+str(sm)+".pdf")
    plt.close()


ingrad   = sys.argv[1]
topology = sys.argv[2] # Start, Through
sm = float(sys.argv[3]) # smooth 
outspl   = sys.argv[4]
outplots = sys.argv[5]


Spline("Amp",0,10,sm)
Spline("Phs",1,10,sm)

