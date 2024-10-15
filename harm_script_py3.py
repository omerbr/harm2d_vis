import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
from matplotlib import rc
from matplotlib.patches import Ellipse
from scipy.interpolate import interp1d
from matplotlib.gridspec import GridSpec
from matplotlib import cm,ticker
from numpy import sin, cos, tan, pi
import numpy as np
import matplotlib.animation as animation

#matplotlib.use('Agg') #so that it does ok with graphics in batch mode

if (np.size(sys.argv) == 1):
    dumps_path = "dumps/" # by default
else:
    dumps_path = sys.argv[1] + "/"
    
print("loading from " + dumps_path + "\n")

#choose Computer Modern Roman fonts by default
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'cmr10'
mpl.rcParams['font.sans-serif'] = 'cmr10'
mpl.rcParams['font.sans-serif'] = 'cmr10'
mpl.rcParams['axes.unicode_minus']=False

#font = { 'size'   : 20}
#rc('font', **font)
rc('xtick', labelsize=20) 
rc('ytick', labelsize=20) 
#rc('xlabel', **font) 
#rc('ylabel', **font) 

# legend = {'fontsize': 20}
# rc('legend',**legend)
axes = {'labelsize': 20}
rc('axes', **axes)
rc('mathtext',fontset='cm')
#use this, but at the expense of slowdown of rendering
#rc('text', usetex=True)
# #add amsmath to the preamble


#matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amssymb,amsmath}"]
import pdb

import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from scipy.interpolate import griddata
from scipy.interpolate import interp1d
from matplotlib import cm,ticker
from numpy import ma
import matplotlib.colors as colors
#use_math_text = True


NEUTRON_STAR = True

def mov_ns_accr(startn=0, endn=20000, dn = 1):
    rg("gdump")
    for i in range(startn,endn,dn):
        file = "ns_accr_%dx%dx%d_%04d.png" % (nx,ny,nz,i)
        if i != startn and os.path.isfile( file ):
            print("File %s exists, skipping..." % file)
            continue
        rd("dump%03d" % i)
        plco(lrho,isfilled=1,levels=np.linspace(-5,1,100),xy=1,xmax=20,ymax=10,dobh=0,cb=1,extend="both",pretty=1,cbxla=r"$\ \ \ \ \log_{10}\rho$",xla=r"$x\ [r_g]$",yla=r"$z\ [r_g]$");
        aphi = psicalc();
        if i == startn:
            maxaphi = np.abs(aphi).max()
        plc(aphi,isfilled=0,levels=np.linspace(0,maxaphi,20)[1:],xy=1,xmax=20,ymax=10,dobh=0,colors="w");
        plc(v1p,isfilled=0,levels=(0,),xy=1,xmax=20,ymax=10,dobh=0,colors="m",linewidths=2)
        mathify_axes_ticks(plt.gca())
        plt.title(r"$t=%g$"%np.round(t),fontsize=20)
        plt.draw()
        plt.savefig(file,bbox_inches='tight',pad_inches=0.02,  dpi=300)

def plot_nersc_charmer_scaling(fntsize=20,dosavefig=1):
    systems = [r"Hopper ($16\times16\times16$ per tile)",
               r"Hopper ($8\times16\times16$ per tile)",
               r"Edison ($32\times32\times32$ per tile)",
               r"Edison ($16\times16\times16$ per tile)"]
    plottypes = ["r--v","g-s","b-.^","m:o"]
    for whichsystem,plottype in zip(systems,plottypes):
        if whichsystem == r"Edison ($32\times32\times32$ per tile)":
            n_list = [24,        192,     1536,    12288]
            s_list = [133966, 130852,   124527,   119265 ]
        elif whichsystem == r"Edison ($16\times16\times16$ per tile)":
            n_list = [24,        192,     1536,    12288]
            s_list = [118656, 110196,     99272.9,    81561.1 ]
        elif whichsystem == r"Hopper ($16\times16\times16$ per tile)":
            n_list = [24,      192,          1536,    12288]
            s_list = [61455.4, 61390.9,   60288.5,  59552.2]
        elif whichsystem == r"Hopper ($8\times16\times16$ per tile)":
            n_list = [24,      192,          1536,    12288]
            s_list = [55879.9, 55068.6,   55187.3,  47352.6]
        n_list = np.array(n_list,dtype=np.float64)
        s_list = np.array(s_list,dtype=np.float64)
        n_array = 10**np.linspace(1,5,100)
        plt.plot(n_list, 100.*s_list/s_list[0],plottype,lw=2,label=whichsystem,markersize=10)
    plt.xscale("log")
    plt.yscale("linear")
    plt.xlim(10,1e5)
    plt.ylim(0,110)
    mathify_axes_ticks(plt.gca())
    plt.xlabel(r"${\rm Number\ of\ cores}$",fontsize=fntsize)
    plt.ylabel(r"${\rm Parallel\ efficiency\ [\%]}$",fontsize=fntsize)
    plt.title(r"${\rm CHARMER\ weak\ scaling\ at\ NERSC}$",fontsize=fntsize)
    plt.grid()
    leg = plt.legend(loc="lower right")
    ax = plt.gca()
    for label in ax.get_xticklabels() + ax.get_yticklabels() + leg.get_texts():
        label.set_fontsize(fntsize)
    if dosavefig:
        plt.savefig("nersc_scaling.eps",
                bbox_inches='tight',pad_inches=0.06,dpi=300)
    

def plot_mad_xy_slice():
    rg("gdump")
    rd("dump1332")
    #density
    plco(np.log10(rho),levels=np.linspace(-3,1.,100),isfilled=1,xy=2,symmx=1,xmax=25,ymax=25,cb=1,pretty=1,dobh=1,domathify=1,extend="both",cbxla=r"$\log_{10}\rho$")
    plc(np.log10(rho),levels=(-0.5,),isfilled=0,xy=2,symmx=1,xmax=25,ymax=25,pretty=1,dobh=1,domathify=1,cbxla=r"$\log_{10}\rho$",colors="k",linewidths=2)    
    plt.title(r"$a=0.5$ MAD, $t=%g$"%t,fontsize=20)
    plt.ylabel(r"$y\ [r_g]$")
    plt.xlabel(r"$x\ [r_g]$")
    plt.savefig("lrho_a05_MAD_t%g.png" % np.round(t),dpi=300,bbox_inches='tight',pad_inches=0.02)
    #bsq/rho
    plco(np.log10(bsq/rho),levels=np.linspace(-3,2,100),isfilled=1,xy=2,symmx=1,xmax=25,ymax=25,cb=1,dobh=1,cbxla=r"$\log_{10}b^2\!/\rho$",pretty=1,extend="both")
    plc(np.log10(rho),levels=(-0.5,),isfilled=0,xy=2,symmx=1,xmax=25,ymax=25,pretty=1,dobh=1,domathify=1,cbxla=r"$\log_{10}\rho$",colors="k",linewidths=2)
    plt.title(r"$t=%g$"%t)
    plt.title(r"$t=%g$"%t,fontsize=20)
    plt.title(r"$a=0.5$ MAD, $t=%g$"%t,fontsize=20)
    plt.ylabel(r"$y\ [r_g]$")
    plt.xlabel(r"$x\ [r_g]$")
    plt.savefig("bsqorho_a05_MAD_t%g.png" % np.round(t),dpi=300,bbox_inches='tight',pad_inches=0.02)
    #bsq/ug
    plco(np.log10(bsq/ug),levels=np.linspace(-2,2,100),isfilled=1,xy=2,symmx=1,xmax=25,ymax=25,cb=1,dobh=1,cbxla=r"$\log_{10}u_g$",pretty=1,extend="both")
    plc(np.log10(rho),levels=(-0.5,),isfilled=0,xy=2,symmx=1,xmax=25,ymax=25,pretty=1,dobh=1,domathify=1,cbxla=r"$\log_{10}u_g$",colors="k",linewidths=2)
    plt.title(r"$t=%g$"%t)
    plt.title(r"$t=%g$"%t,fontsize=20)
    plt.title(r"$a=0.5$ MAD, $t=%g$"%t,fontsize=20)
    plt.ylabel(r"$y\ [r_g]$")
    plt.xlabel(r"$x\ [r_g]$")
    plt.savefig("bsqoug_a05_MAD_t%g.png" % np.round(t),dpi=300,bbox_inches='tight',pad_inches=0.02)
    #Theta_e,c
    plco(np.log10(Tel4*1826.15),levels=np.linspace(0,4,100),isfilled=1,xy=2,symmx=1,xmax=25,ymax=25,cb=1,dobh=1,cbxla=r"$\log_{10}\Theta_{\rm e,c}$",pretty=1,extend="both")
    plc(np.log10(rho),levels=(-0.5,),isfilled=0,xy=2,symmx=1,xmax=25,ymax=25,pretty=1,dobh=1,domathify=1,cbxla=r"$\log_{10}u_g$",colors="k",linewidths=2)
    plt.title(r"$t=%g$"%t)
    plt.title(r"$t=%g$"%t,fontsize=20)
    plt.title(r"$a=0.5$ MAD, $t=%g$"%t,fontsize=20)
    plt.ylabel(r"$y\ [r_g]$")
    plt.xlabel(r"$x\ [r_g]$")
    plt.savefig("thetaec_a05_MAD_t%g.png" % np.round(t),dpi=300,bbox_inches='tight',pad_inches=0.02)
    #Theta_e
    plco(np.log10(Teldis*1826.15),levels=np.linspace(0,4,100),isfilled=1,xy=2,symmx=1,xmax=25,ymax=25,cb=1,dobh=1,cbxla=r"$\log_{10}\Theta_{\rm e,nc}$",pretty=1,extend="both")
    plc(np.log10(rho),levels=(-0.5,),isfilled=0,xy=2,symmx=1,xmax=25,ymax=25,pretty=1,dobh=1,domathify=1,cbxla=r"$\log_{10}u_g$",colors="k",linewidths=2)
    plt.title(r"$t=%g$"%t)
    plt.title(r"$t=%g$"%t,fontsize=20)
    plt.title(r"$a=0.5$ MAD, $t=%g$"%t,fontsize=20)
    plt.ylabel(r"$y\ [r_g]$")
    plt.xlabel(r"$x\ [r_g]$")
    plt.savefig("thetaenc_a05_MAD_t%g.png" % np.round(t),dpi=300,bbox_inches='tight',pad_inches=0.02)
    #u^r
    plco(np.log10(uu[1]*dxdxp[1,1]),levels=np.linspace(-1,1,100),isfilled=1,xy=2,symmx=1,xmax=25,ymax=25,cb=1,dobh=1,cbxla=r"$\log_{10}u^r$",pretty=1)
    plc(np.log10(rho),levels=(-0.5,),isfilled=0,xy=2,symmx=1,xmax=25,ymax=25,pretty=1,dobh=1,domathify=1,cbxla=r"$\log_{10}u_g$",colors="k",linewidths=2)
    plt.title(r"$t=%g$"%t)
    plt.title(r"$t=%g$"%t,fontsize=20)
    plt.title(r"$a=0.5$ MAD, $t=%g$"%t,fontsize=20)
    plt.ylabel(r"$y\ [r_g]$")
    plt.xlabel(r"$x\ [r_g]$")
    plt.savefig("ur_a05_MAD_t%g.png" % np.round(t),dpi=300,bbox_inches='tight',pad_inches=0.02)
    #heat flux
    heatflux = phi*(bu*ud[0]+uu*bd[0])/bsq**0.5
    plco(np.log10(heatflux[1]*dxdxp[1,1]),levels=np.linspace(-2,1,100),isfilled=1,xy=2,symmx=1,xmax=25,ymax=25,cb=1,dobh=1,cbyla=r"$\log_{10}(q^r\!u_t+u^r\!q_t)$",pretty=1,extend="both")
    plc(np.log10(rho),levels=(-0.5,),isfilled=0,xy=2,symmx=1,xmax=25,ymax=25,pretty=1,dobh=1,domathify=1,cbxla=r"$\log_{10}u_g$",colors="k",linewidths=2)
    plt.title(r"$t=%g$"%t)
    plt.title(r"$t=%g$"%t,fontsize=20)
    plt.title(r"$a=0.5$ MAD, $t=%g$"%t,fontsize=20)
    plt.ylabel(r"$y\ [r_g]$")
    plt.xlabel(r"$x\ [r_g]$")
    plt.savefig("radial_heat_flux_a05_MAD_t%g.png" % np.round(t),dpi=300,bbox_inches='tight',pad_inches=0.02)
    #passive floor scalar
    if DOFLR:
        plco((flr),levels=np.linspace(0,1,100),isfilled=1,xy=2,symmx=1,xmax=25,ymax=25,cb=1,dobh=1,cbxla=r"$\log_{10}$flr",pretty=1,extend="both")
        plc(np.log10(rho),levels=(-0.5,),isfilled=0,xy=2,symmx=1,xmax=25,ymax=25,pretty=1,dobh=1,domathify=1,colors="k",linewidths=2)
        plt.title(r"$t=%g$"%t)
        plt.title(r"$t=%g$"%t,fontsize=20)
        plt.title(r"$a=0.5$ MAD, $t=%g$"%t,fontsize=20)
        plt.ylabel(r"$y\ [r_g]$")
        plt.xlabel(r"$x\ [r_g]$")
        plt.savefig("floor_scalar_a05_MAD_t%g.png" % np.round(t),dpi=300,bbox_inches='tight',pad_inches=0.02)

def plot_cond_slice():
    os.chdir(os.path.join(os.environ["HOME"],"run/a05mad"))
    print(os.getcwd())
    rg("gdump")
    rd("dump1360")
    plt.plot(r[:,0,0],(Tel4*1836*1e-3)[:,ny/2,0],"r",label=r"$\times10^{-3}$, $\alpha_{\rm e}=10$, w/cond, no suppress cond")
    #alpha_e = 0.1
    os.chdir(os.path.join(os.environ["HOME"],"run/a05mad_alpha01"))
    print(os.getcwd())
    rd("dump1360")
    plt.plot(r[:,0,0],(Tel4*1836*1e-3)[:,ny/2,0],"orange",label=r"$\times10^{-3}$, $\alpha_{\rm e} = 0.1$, w/cond, no suppress cond")
    #bsq/rho > 20
    os.chdir(os.path.join(os.environ["HOME"],"run/a05madsavio_alpha10_suppressphi5"))
    print(os.getcwd())
    rd("dump1360")
    plt.plot(r[:,0,0],(Tel4*1836)[:,ny/2,0],"g",label=r"w/conduction, suppress cond $b^2\!/\rho>20$ or $b^2\!/u_g>200$")
    #bsq/rho > 10
    os.chdir(os.path.join(os.environ["HOME"],"run/a05madsavio_suppressphi3"))
    print(os.getcwd())
    rd("dump1360")
    plt.plot(r[:,0,0],(Tel4*1836)[:,ny/2,0],"b",label=r"w/conduction, suppress cond $b^2\!/\rho>10$ or $b^2\!/u_g>100$")
    #bsq/rho > 5
    os.chdir(os.path.join(os.environ["HOME"],"run/a05madsavio_alpha10_suppressphi4"))
    print(os.getcwd())
    rd("dump1360")
    plt.plot(r[:,0,0],(Tel4*1836)[:,ny/2,0],"m",label=r"w/conduction, suppress cond $b^2\!/\rho>5$ or $b^2\!/u_g>50$")
    #no conduction in previously floored regions
    os.chdir(os.path.join(os.environ["HOME"],"run/a05madsavio_flrscalar"))
    print(os.getcwd())
    rd("dump1360")
    plt.plot(r[:,0,0],(Tel4*1836)[:,ny/2,0],"c",label=r"no conduction in prev. floored regions via passive scalar")
    plt.plot(r[:,0,0],(Teldis*1836)[:,ny/2,0],"k",label=r"no conduction")
    plt.xlim(rhor,30)
    plt.ylim(0.1,100)
    plt.grid(b=1)
    plt.xscale("log")
    plt.yscale("log")
    plt.legend(loc="best",frameon=0)
    plt.ylabel(r"$\Theta_{\rm e}$",fontsize=20)
    plt.xlabel(r"$r\ [r_g]$",fontsize=20)
    plt.title(r"$a = %g$ MAD, $t = %g$, slice along $x-$axis" % (a, t),fontsize=20)
    plt.savefig("../conduction_suppression_comparison.png",dpi=200,bbox_inches='tight',pad_inches=0.02)

def eliot_theta_movie():
    firsttime = 1
    rg("gdump")
    for i in np.arange(0,1591):
        rd("dump%03d" % i); print (i); sys.stdout.flush()
        eliot_thetae_plot()
        if firsttime: 
            plt.tight_layout()
            firsttime = 0
        plt.gca().set_aspect("1.0")
        plt.draw()
        plt.savefig("frame%04d.png"%i,dpi=180,bbox_inches='tight',pad_inches=0.02)

def eliot_thetae_plot(sf = 0):
    mrat = 1836.15267245;
    Thetae = Tel4*mrat
    #reset Thetae to nans outside of the range
    Thetae[bsq/rho>10] = np.nan+Thetae[bsq/rho>10]
    res,cb=plco(np.log10(Thetae),xy=1,xmax=10,ymax=5,symmx=1,nc=100,isfilled=1,cb=1,pretty=1,levels=np.linspace(-1.,1.7,500),extend="both",cbticks=[-1,-0.5,0,0.5,1,1.5],cbxla=r"$\ \ \ \ \ \ \ \log_{10}\Theta_{\rm e}$",xla=r"$x\ [r_g]$",yla=r"$z\ [r_g]$")
    aphi = psicalc()
    plc(aphi,symmx=1,xy=-1,levels=np.linspace(0,0.3,20)[1:],colors="black",linewidths=1.)
    plt.xlim(-11,11)
    plt.ylim(-10,10)
    mathify_axes_ticks(plt.gca())
    if sf: plt.savefig("Thetae_a05_1loop_t15000.png",dpi=300,bbox_inches='tight',pad_inches=0.02)

def mathify_axes_ticks(ax,fontsize=20,xticks=None,yticks=None):
    if xticks is None:
        xticks = ax.get_xticks()
    if yticks is None:
        yticks = ax.get_yticks()
    if ax.get_xscale() != 'log': ax.set_xticklabels([(r'$%g$' % lab) for lab in xticks])
    if ax.get_yscale() != 'log': ax.set_yticklabels([(r'$%g$' % lab) for lab in yticks])
    if fontsize is not None:
        if ax.get_xscale() != 'log':
            for label in ax.get_xticklabels():
                label.set_fontsize(fontsize)
        if ax.get_yscale() != 'log':
            for label in ax.get_yticklabels():
                label.set_fontsize(fontsize)


def plot_thetae_beta_avg(fntsize=20,rx = 5,sf=0,fno=1,lw=3,docompareavgs=0,dosf=1,docond=1):
    plt.figure(fno,figsize=(6,8))
    plt.clf()
    runname = os.path.basename(os.getcwd())
    rg("gdump")
    #defalt with conduction, but if conduction is broken, can switch to no conduction
    if docond:
        whichTel = "Tel4"
        whichfel = "fel4"
    else:
        whichTel = "Teldis"
        whichfel = "feldis"
    if runname == "a0.5_mad_aphipow4_beta100_resetphi0ifnan":
        #NEEDS RECOMPUTE; STARTED MKAVG
        avgs = [
                np.load("avg2d_00190_00199_100.npz"),
                np.load("avg2d_00130_00140_100.npz"),
                np.load("avg2d_00110_00120_100.npz"),
                np.load("avg2d_00060_00070_100.npz")]
        title = r"Big loop, small torus, $a = %g$" % a
    elif runname == "a0.5_1loop":
        #NEEDS RECOMPUTE; STARTED MKAVG
        #this one has broken phi and ugel4
        avgs = [np.load("avg2d_00196_00205_100.npz"),
                np.load("avg2d_00160_00170_100.npz"),
                np.load("avg2d_00140_00150_100.npz"),
                np.load("avg2d_00110_00120_100.npz")]
        title = r"1-loop, small torus, $a = %g$" % a
    elif runname == "a05mad":
        avgs = [np.load("avg2d_00140_00200_100.npz")]
        title = r"MAD, large torus, $a = %g$" % a
    elif runname == "1loop_fevar_bsqou2500":
        #NEEDS RECOMPUTE; STARTED MKAVG
        avgs = [np.load("avg2d_00090_00100_100.npz"),
                np.load("avg2d_00110_00120_100.npz"),
                np.load("avg2d_00148_00158_100.npz")]
        title = r"1-loop, small torus, $a = %g$" % a
    elif runname == "2loop_fe05_avx_firstfloortry":
        avgs = [np.load("avg2d_00035_00057_100.npz")]
        title = r"$f_e$ = 0.5, 2-loop, small torus, a = %g" % a
    elif runname == "2loop_fevar_avx_firstfloortry_bsqou2500":
        avgs = [np.load("avg2d_00065_00075_100.npz",
                        "avg2d_00098_00108_100.npz",
                        "avg2d_00140_00170_100.npz")]
        title = r"2-loop, small torus, a = %g" % a
    title+=", %sconduction" % ("w/" if docond else "no ")
    mrat = 1836.152672
    cols = ["red", "orange", "green", "blue", "cyan", "magenta", "black"]
    lts  = [None, [10,5], [10,2,3,2], [2,2]]
    rhor = 1+(1-a**2)**0.5
    fracphi = dxdxp[3,3,0,0,0]*_dx3*nz/(2*np.pi)
    gs = GridSpec(3, 1)
    gs.update(left=0.15, right=0.99, top=0.95, bottom=0.08, wspace=0.06, hspace=0.06)
    for avg,col,lt in zip(avgs,cols,lts):
        #NOTE: need to use absB[0] to get B^r because B[0] = B^t is skipped
        PhiBH = 0.5*(4*pi)**0.5*(gdet*avg["absB"][0]*_dx2*_dx3*nz).mean(-1).sum(-1)/fracphi
        Mdot = (-gdet*avg["rhouu"][1]*_dx2*_dx3*nz).mean(-1).sum(-1)/fracphi
        ix = iofr(rx)
        ih = iofr(rhor)
        phibh = PhiBH/Mdot**0.5
        Mdot = Mdot[ix]
        phibh = phibh[ih]
        x = r[:,ny/2,0]
        y = (gdet*avg["rho%s"%whichTel]*mrat*_dx2*_dx3*nz).mean(-1).sum(-1)/(gdet*avg["rho"]*_dx2*_dx3*nz).mean(-1).sum(-1)
        vr = avg["uu"][1,:,ny/2,0]/avg["uu"][0,:,ny/2,0]*gcov[1,1,:,ny/2,0]**0.5
        tacc = -r[:,ny/2,0]/vr
        req = r[(avg["t"][-1]<2*tacc),ny/2,0][0]
        #
        #top panel (0,0): Theta_e
        #
        label = r"$[%g,%g), \phi_{\rm H} = %g$" % (
                        np.round(avg["t"][0]),
                        np.round(avg["t"][-1]*2-avg["t"][-2]),
                        np.round(phibh,decimals=1))
        if docompareavgs:
            label = r"$\langle \rho\Theta_{\rm e}\rangle_{\theta\varphi}/\langle\rho\rangle_{\theta\varphi}$"
        ax1 = plt.subplot(gs[0,0])
        l1, = plt.plot(x[x<=req],y[x<=req], lw = lw, color=col, label=label)
        l2, = plt.plot(x[x>=req],y[x>=req], lw = lw, color=col, alpha = 0.5)
        if lt is not None: [l.set_dashes(lt) for l in [l1,l2]]
        if docompareavgs:
            y1 = (4./3.-1)*avg["ug"+whichTel[1:]][:,ny/2,0]/avg["rho"][:,ny/2,0]*mrat
            y2 = avg[whichTel][:,ny/2]*mrat
            plt.plot(x,y1, lw = 1.5, color=col,
                     label=r"$\langle \rho \Theta_{\rm e}\rangle_{z=0}/\langle \rho \rangle_{z=0}$")
            plt.plot(x,y2, lw = 0.5, color=col,label=r"$\langle \Theta_{\rm e}\rangle_{z=0}$")
        #
        #middle panel (1,0): beta
        #
        y = (gdet*avg["rhobeta"]*_dx2*_dx3*nz).mean(-1).sum(-1)/(gdet*avg["rho"]*_dx2*_dx3*nz).mean(-1).sum(-1)
        ax2 = plt.subplot(gs[1,0])
        l1, = plt.plot(x[x<=req],y[x<=req],
                 lw = lw,
                 color=col,
                 label=r"$\langle\rho\beta\rangle_{\theta\varphi}/\langle\rho\rangle_{\theta\varphi}$")
        l2, = plt.plot(x[x>=req],y[x>=req],
                 lw = lw,
                 color=col,
                 alpha = 0.5)
        if lt is not None: [l.set_dashes(lt) for l in [l1,l2]]
        if docompareavgs:
            y2 = avg["rhobeta"][:,ny/2,0]/avg["rho"][:,ny/2,0] 
            y3 = avg["beta"][:,ny/2,0]
            y4 = (2*(gam-1)*avg["ug"]/avg["bsq"])[:,ny/2,0]
            plt.plot(x,y2, lw = 1.5, color=col,
                     label=r"$\langle \rho\beta\rangle_{z=0}/\langle \rho\rangle_{z=0}$")
            plt.plot(x,y3, lw = 0.5, color=col,
                     label=r"$\langle\beta\rangle_{z=0}$")
            l,=plt.plot(x,y4, lw = 0.5, color=col,
                     label=r"$\langle p_{\rm g}\rangle_{z=0}/\langle p_{\rm m}\rangle_{z=0}$")
            l.set_dashes([10,5])
            plt.legend(loc="best",frameon = 0,ncol=2)
        #
        #bottom panel (2,0): f_e
        #
        ax3 = plt.subplot(gs[2,0])
        y = (gdet*avg["rho%s" % whichfel]*_dx2*_dx3*nz).mean(-1).sum(-1)/(gdet*avg["rho"]*_dx2*_dx3*nz).mean(-1).sum(-1)
        l1, = plt.plot(x[x<=req],y[x<=req],
                 lw = lw,
                 color=col,
                 label=r"$\langle\rho f_{\rm e}\rangle_{\theta\varphi}/\langle\rho\rangle_{\theta\varphi}$")
        l2, = plt.plot(x[x>=req],y[x>=req],
                 lw = lw,
                 color=col,
                 alpha = 0.5)
        if lt is not None: [l.set_dashes(lt) for l in [l1,l2]]
        if docompareavgs:
            y2 = avg["rho%s" % whichfel][:,ny/2,0]/avg["rho"][:,ny/2,0]
            y3 = avg[whichfel][:,ny/2,0]
            plt.plot(x,y2,
                lw = 1.5,
                color=col,
                label=r"$\langle\rho f_{\rm e}\rangle_{z=0}/\langle\rho\rangle_{z=0}$")
            plt.legend(loc="best",frameon=0)
            plt.plot(x,y3,
                lw = 0.5,
                color=col,
                label=r"$\langle f_e\rangle_{z=0}$")
            plt.legend(loc="best",frameon=0)
        if docompareavgs: break
    #top panel: Theta_e
    plt.setp( ax1.get_xticklabels(), visible=False )
    ax1.set_ylim(0.5,20)
    ax1.set_ylabel(r"$\Theta_{\rm e}$",fontsize=fntsize,labelpad=0)
    ax1.legend(loc="best",frameon=0)
    #middle panel: beta
    plt.setp( ax2.get_xticklabels(), visible=False )
    ax2.set_ylabel(r"$\beta$",fontsize=fntsize,labelpad=0)
    ax2.set_ylim(1,1000)
    #bottom panel: f_e
    ax3.set_ylabel(r"$f_{\rm e}$",fontsize=fntsize,labelpad=0)
    ax3.set_xlabel(r"$r\ [r_g]$",fontsize=fntsize)
    ax3.set_ylim(1e-2,0.5)
    #all panels
    if runname == "a0.5_mad_aphipow4_beta100_resetphi0ifnan":
        ax1.set_ylim(0.8,30)
        ax2.set_ylim(1,1e3)
        ax3.set_ylim(0.05,1)
    if runname == "a05mad":
        if docond:
            ax1.set_ylim(50,2e3)
        else:
            ax1.set_ylim(2,3e2)
        ax2.set_ylim(1,100)
        ax3.set_ylim(0.05,2)
    if runname == "a0.5_1loop":
        ax1.set_ylim(0.1,10)
        ax2.set_ylim(1,3e3)
        ax3.set_ylim(0.005,1)
    if runname == "2loop_fe05_avx_firstfloortry":
        ax1.set_ylim(1,50)
    for ax in [ax1,ax2,ax3]:
        ax.set_xlim(rhor,30)
        ax.set_yscale("log")
        ax.set_xscale("log")
        ax.grid(b=1)
    for avg in avgs:
        avg.close()
    ax1.set_title(title)
    if dosf: plt.savefig(os.path.join("..","%s%s%s.png" % (runname, "compavgs" if docompareavgs else "", "nocond" if not docond else "")))
    

def convert_to_single_file(startn=0,endn=-1,ln=10,whichi=0,whichn=1,**kwargs):
    which = kwargs.pop("which","convert_file")
    rg("gdump")
    flist1 = np.sort(glob.glob( os.path.join(dumps_path, "dump[0-9][0-9][0-9]_0000") ) )
    flist2 = np.sort(glob.glob( os.path.join(dumps_path, "dump[0-9][0-9][0-9][0-9]_0000") ) )
    flist1.sort()
    flist2.sort()
    flist = np.concatenate((flist1,flist2))
    firsttime = 1
    for fldname in flist:
        #find the index of the file
        fldindex = np.int(fldname.split("_")[0].split("p")[-1])
        if fldindex < startn:
            continue
        if endn>=0 and fldindex >= endn:
            break
        if fldindex % whichn != whichi:
            #do every whichn'th snapshot starting with whichi'th snapshot
            continue
        #print( "Reading " + fldname + " ..." )
        fname = "dump%03d" % fldindex
        if os.path.isfile( fname ):
            print("File %s exists, skipping..." % fname)
            continue
        if not os.path.isfile( fname ):
            rd(fname)
        
#nuc macros begin

def plot_torii(sf=0):
    #
    #density plot
    #
    plt.figure(1)
    plt.clf()
    res = np.loadtxt("dens_ascii.txt")
    r_nr,th_nr,uknown,rho_nr = res.T
    plt.ylabel(r"$\rho$",fontsize=20)
    w = (th_nr<np.pi/2+0.02)&(th_nr>np.pi/2-0.02)
    #rodrigo density (normalize it to unity at max)
    plt.plot(r_nr[w],rho_nr[w],":",lw=2,label = "NR torus")
    #sasha density
    rd("dump000")
    plt.plot(r[:,ny/2,0],rho[:,ny/2,0],lw=2,label="GR torus")
    plt.xlim(1,100)
    plt.ylim(1e-15,1e-5)
    plt.xscale("log")
    plt.yscale("log")
    plt.legend(loc="best")
    plt.grid(b=1,which="major")
    plt.grid(b=1,axis="x",which="minor")
    plt.ylabel(r"$\rho\ [M_{\rm BH}/r_g^3]$",fontsize=20)
    plt.xlabel(r"$r\ [r_g]$",fontsize=20)
    if sf: plt.savefig("dens_vs_r.pdf")
    #
    #transverse density plot
    #
    plt.figure(2)
    plt.clf()
    res = np.loadtxt("dens_ascii.txt")
    r_nr,th,uknown,rho_nr = res.T
    plt.ylabel(r"$\rho$",fontsize=20)
    rx = 11.2
    w = (r_nr<rx+0.2)&(r_nr>rx-0.2)
    #rodrigo density (normalize it to unity at max)
    plt.plot(th_nr[w]-np.pi/2,rho_nr[w],":",lw=2,label = "NR torus")
    #sasha density
    rd("dump000")
    ix = iofr(rx)
    plt.plot(h[ix,:,0]-np.pi/2.,rho[ix,:,0],lw=2,label="GR torus")
    plt.xlim(-1,+1)
    plt.xscale("linear")
    plt.yscale("log")
    plt.legend(loc="best")
    plt.grid(b=1)
    plt.ylabel(r"$\rho\ [M_{\rm BH}/r_g^3]$",fontsize=20)
    plt.xlabel(r"$\theta-\pi/2$",fontsize=20)
    if sf: plt.savefig("dens_vs_theta.pdf")
    #
    #angular momentum plot
    #
    plt.figure(3)
    plt.clf()
    res = np.loadtxt("angz_ascii.txt")
    r_nr,th,ell = res.T
    #rodrigo angular momentum
    w = (th_nr<np.pi/2+0.02)&(th_nr>np.pi/2-0.02)
    plt.plot(r_nr[w],ell[w],":",lw=2,label="NR torus")
    #sasha angular momentum
    rd("dump000")
    plt.plot(r[:,ny/2,0],(-ud[3]/ud[0])[:,ny/2,0],lw=2,label="GR torus")
    plt.plot(r[:,ny/2,0],ellk(a,r[:,ny/2,0]),"--",lw=2,label="Keplerian angular momentum")
    plt.ylabel(r"$\ell \approx rv_{\rm \varphi}\ [r_g c]$",fontsize=20)
    plt.xlabel(r"$r\ [r_g]$",fontsize=20)
    plt.xlim(1,100)
    plt.legend(loc="best")
    plt.grid(b=1,which="both")
    plt.xscale("log")
    plt.yscale("log")
    plt.ylim(1,10)
    if sf: plt.savefig("angz.pdf")
    #
    #Temperature plot
    #
    plt.figure(4)
    plt.clf()
    res = np.loadtxt("temp_ascii.txt")
    r_nr,th_nr,T_nr = res.T
    #rodrigo angular momentum
    w = (th_nr<np.pi/2+0.02)&(th_nr>np.pi/2-0.02)
    plt.plot(r_nr[w],T_nr[w],":",lw=2,label="NR torus")
    #sasha angular momentum
    rd("dump000")
    plt.plot(r[:,ny/2,0],Tnuc_cgs[:,ny/2,0],lw=2,label="GR torus")
    plt.ylabel(r"$T\ {\rm [K]}$",fontsize=20)
    plt.xlabel(r"$r\ [r_g]$",fontsize=20)
    plt.xlim(1,100)
    plt.legend(loc="best")
    plt.grid(b=1,which="both")
    plt.xscale("log")
    plt.yscale("log")
    plt.ylim(5e8,5e10)
    if sf: plt.savefig("temp.pdf")
    #
    #Degeneracy plot
    #
    plt.figure(5)
    plt.clf()
    res = np.loadtxt("etae_ascii.txt")
    r_nr,th_nr,etae_nr = res.T
    #rodrigo angular momentum
    w = (th_nr<np.pi/2+0.02)&(th_nr>np.pi/2-0.02)
    plt.plot(r_nr[w],etae_nr[w],":",lw=2,label="NR torus")
    #sasha angular momentum
    rd("dump000")
    plt.plot(r[:,ny/2,0],etae[:,ny/2,0],lw=2,label="GR torus")
    plt.ylabel(r"$\eta_{\rm e}$",fontsize=20)
    plt.xlabel(r"$r\ [r_g]$",fontsize=20)
    plt.xlim(1,100)
    plt.legend(loc="best")
    plt.grid(b=1,which="both")
    plt.xscale("log")
    plt.yscale("log")
    plt.ylim(1e-1,10)
    if sf: plt.savefig("etae.pdf")
    #
    #Q plot (nuclear heating)
    #
    plt.figure(6)
    plt.clf()
    res = np.loadtxt("ncfn_ascii.txt")
    r_nr,th_nr,etae_nr = res.T
    w = (th_nr<np.pi/2+0.02)&(th_nr>np.pi/2-0.02)
    #rodrigo
    plt.plot(r_nr[w],-etae_nr[w],":",lw=2,label="NR torus")
    print(etae_nr[w])
    #sasha
    rd("dump001")
    plt.plot(r[:,ny/2,0],-Q[:,ny/2,0],lw=2,label="GR torus")
    plt.ylabel(r"$-Q$",fontsize=20)
    plt.xlabel(r"$r\ [r_g]$",fontsize=20)
    plt.xlim(1,100)
    plt.legend(loc="best")
    plt.grid(b=1,which="both")
    plt.xscale("log")
    plt.yscale("log")
    plt.ylim(1e-10,1e-4)
    if sf: plt.savefig("Q.pdf")
    #
    #G plot (neutrino emission)
    #
    plt.figure(7)
    plt.clf()
    res = np.loadtxt("gamn_ascii.txt")
    r_nr,th_nr,etae_nr = res.T
    w = (th_nr<np.pi/2+0.02)&(th_nr>np.pi/2-0.02)
    #rodrigo
    plt.plot(r_nr[w],etae_nr[w],":",lw=2,label="NR torus")
    print(etae_nr[w])
    #sasha
    rd("dump001")
    plt.plot(r[:,ny/2,0],G[:,ny/2,0],lw=2,label="GR torus")
    plt.ylabel(r"$\Gamma$",fontsize=20)
    plt.xlabel(r"$r\ [r_g]$",fontsize=20)
    plt.xlim(1,100)
    plt.legend(loc="best")
    plt.grid(b=1,which="both")
    plt.xscale("log")
    plt.yscale("log")
    plt.ylim(1e-7,1e-5)
    if sf: plt.savefig("G.pdf")

#nuc macros end

def ellk(a,r):
    ekval = ek(a,r)
    lkval = lk(a,r)
    return(lkval/ekval)

def ek(a,r):
    #-u_t, I suspect
    ek = (r**2-2*r+a*r**0.5)/(r*(r**2-3*r+2*a*r**0.5)**0.5)
    return(ek)

def lk(a,r):
    udphi = r**0.5*(r**2-2*a*r**0.5+a**2)/(r*(r**2-3*r+2*a*r**0.5)**0.5)
    return( udphi )

def felcalc(ugel=None):
    if "game4" not in globals(): game4 = 4./3.
    SMALL = 1e-20
    mrat = 1836.152672
    #Calculations for fel
    c1 = .92 #heating constant
        
    Tpr = (gam-1.)*ug/rho;
    if ugel is None:
        ugel = ugel4
    Tel = (game4-1.)*ugel/rho;
    Tel[Tel<=SMALL] = Tel[Tel<=SMALL]*0 + SMALL;
    Trat = np.abs(Tpr/Tel) ;
    
    ppres = rho*Tpr #Proton pressure
    
    beta = ppres/(bsq+SMALL) * 2. ;
    beta[beta > 1.e20] = beta[beta > 1.e20]*0 + 1.e20;
    mbeta = 2.-.2*np.log10(Trat);
    
    #Trat<=1
    c2 = 1.6 / Trat ;
    c3 = 18. + 5.*np.log10(Trat);
    #Trat>1
    c2[Trat>1] = 1.2 / Trat[Trat>1]
    c3[Trat>1] = c3[Trat>1]*0 + 18
    
    c22 = c2**2;
    c32 = c3**2;
    
    qrat = c1 * (c22+(beta**mbeta))/(c32 + (beta**mbeta)) * np.exp(-1./beta)*(mrat*Trat)**.5 ;
    
    fel = 1./(1.+qrat);
    
    return fel;

#effective radius for Noble-like cooling
def Rz_relcorr(a,r):
    risco=Risco(a)
    #cap r from below so does not drop below risco
    rcap = amax(r,risco)
    lkcap = lk(a,rcap)
    ekcap = ek(a,rcap)
    Rz = (lkcap**2-a**2*(ekcap**2-1.))/r
    return( Rz )

#returns T = pgas/rho of a NT disk
def Tstar(a,r,hor):
    return( Rz_relcorr(a,r) * hor**2 / r )

def Tnt(a,r,hor):
    return Tstar(a,r,hor)

def Risco(ain):
    eps = np.finfo(np.float64).eps
    a = np.minimum(ain,1.)
    Z1 = 1 + (1. - a**2)**(1./3.) * ((1. + a)**(1./3.) + (1. - a)**(1./3.))
    Z2 = (3*a**2 + Z1**2)**(1./2.)
    risco = 3 + Z2 - np.sign(a)* ( (3 - Z1)*(3 + Z1 + 2*Z2) )**(1./2.)
    return(risco)

def Ebind(r,a):
    #1+u_t, I suspect
    Eb = 1 - (r**2-2*r+a*r**0.5)/(r*(r**2-3*r+2*a*r**0.5)**0.5)
    return( Eb )

def etaNT(a):
    return( Ebindisco(a) )

def Ebindisco(a):
    eps = np.finfo(np.float64).eps
    a0 = 0.99999 #1.-1e8*eps
    if a > a0: 
        a = a0
        Eb = Ebind( Risco(a), a )
        return((a-a0)/(1.-a0)*(1.-3.**(-0.5)) + (1.-a)/(1.-a0)*Eb)
    Eb = Ebind( Risco(a), a)
    #Eb = (1.-3.**(-0.5))*a**2
    return( Eb )

def mkmov_simple(starti=0,endi=400,di=1):
    for i in range(starti,endi+1,di):
        plt.clf()
        rd("dump%03d" % i);
        aphi=psicalc()
        gs = GridSpec(1, 1)
        gs.update(left=0.108, right=1.02, top=0.95, bottom=0.15, wspace=0.0, hspace=0.0)
        ax = plt.subplot(gs[0,0])
        if i == starti: aphimax = aphi.max()
        mkfrm_simple(aphimax=aphimax,ax=ax)
        plt.draw();
        plt.savefig("frame%04d.png"%i,dpi=225)

def mkfrm_simple_2panel(ax=None,ln=None,lnx=None,lny=None,aphimax=None,fntsize=20,fig=None,
                        xmax=20,ymax=20,asp=1.,vmin=-6,vmax=0,OmegaNS=1./20.,freeze_rotation=1):
        if ax is None:
            plt.clf()
            gs = GridSpec(1, 2)
            gs.update(left=0.08, right=0.92, top=0.94, bottom=0.14, wspace=0.15, hspace=0.08)
            #plt.title(r"$t=$%05g$\,r_g/c$"%np.round(t)); 
            ax1 = plt.subplot(gs[0,0])
        cs = plc(lrho,xy=1,xmax=xmax,ymax=ymax,cmap=cm.jet,isfilled=1,levels=np.linspace(vmin,vmax,100),
                 cb=0,symmx=1,extend="both",pretty=1,dobh=0,ax=ax1)
        ax1.set_xlabel(r"$x\ [r_g]$",fontsize=fntsize,labelpad=0)
        ax1.set_ylabel(r"$z\ [r_g]$",fontsize=fntsize,labelpad=-10)
        aphi = psicalc()
        if aphimax is None:
            aphimax = np.max(aphi)
        plc(aphi,levels=np.linspace(-aphimax*0.7,aphimax*0.7,21),colors="k",
            linewidths=1,xy=-1,symmx=1,alpha=0.5,dobh=0,ax=ax1)
        ax1.set_aspect(asp)
        #Ellipse
        el = Ellipse((0,0), 1.1*2*Rin, 1.1*2*Rin, facecolor='grey', alpha=1)
        art=ax1.add_artist(el)
        art.set_zorder(0)
        #
        #right panel
        #
        ax2 = plt.subplot(gs[0,1])
        cs = plc(lrho,xy=2,xmax=xmax,ymax=ymax,cmap=cm.jet,isfilled=1,levels=np.linspace(vmin,vmax,100),
                 cb=0,symmx=1,extend="both",pretty=1,dobh=0,ax=ax2,delta_phi=-OmegaNS*t)
        ax2.set_xlabel(r"$x\ [r_g]$",fontsize=fntsize,labelpad=0)
        ax2.set_ylabel(r"$y\ [r_g]$",fontsize=fntsize,labelpad=2)
        #cb.ax.set_xlabel(r"$\log\rho$",fontsize=fntsize,ha="left")
        #
        #plc(aphi,levels=np.linspace(-aphimax*0.7,aphimax*0.7,21),colors="k",linewidths=1,xy=-1,symmx=1,alpha=0.5,dobh=0,ax=ax1)
        ax2.set_aspect(asp)
        #
        #Ellipse
        #
        el = Ellipse((0,0), 1.1*2*Rin, 1.1*2*Rin, facecolor='grey', alpha=1)
        art=ax2.add_artist(el)
        art.set_zorder(0)
        #
        if freeze_rotation:
            ax2.plot([0, Rin*np.cos(OmegaNS*0)], [0, Rin*np.sin(OmegaNS*0)], color='k', linestyle='-', linewidth=2)
        else:
            ax2.plot([0, Rin*np.cos(OmegaNS*t)], [0, Rin*np.sin(OmegaNS*t)], color='k', linestyle='-', linewidth=2)
        #
        if fig is None:
            fig = plt.gcf()
        label = r"$\log\rho$"
        cx1,cb1 = mkvertcolorbar(ax2,fig,gap=0.02,width=0.02,vmin=vmin,vmax=vmax,loc="right",
                   label=label,fntsize=fntsize*0.8,cmap=cm.jet,extend="both")
        cb1.ax.set_xlabel(label,fontsize=fntsize,ha="center",x=1.5,labelpad=0)
        plt.setp( ax2.get_yticklabels(), visible=False )
        plt.annotate(r"$t=$%05g$\,r_g/c$"%np.round(t),(0.5,1),xytext=(0,0),xycoords="figure fraction",textcoords="offset points",fontsize=20,rotation=0,rotation_mode="anchor",ha="center",va="top",color="black")

        #Ellipse
        #el = Ellipse((0,0), 1.1*2*Rin, 1.1*2*Rin, facecolor='grey', alpha=1)
        #art=ax1.add_artist(el)
        #art.set_zorder(0)
        
        #print(i);

def mkfrm_simple(ax=None,ln=None,lnx=None,lny=None,aphimax=None,fntsize=20,fig=None,xmax=25,ymax=15,asp=1.):
        if ax is None:
            plt.clf()
            gs = GridSpec(1, 1)
            gs.update(left=0.108, right=1.02, top=0.95, bottom=0.15, wspace=0.0, hspace=0.0)
            ax = plt.subplot(gs[0,0])
        cs, cb = plco(lrho,xy=1,xmax=xmax,ymax=ymax,cmap=cm.jet,isfilled=1,levels=np.linspace(-6,0,100),cb=1,symmx=1,extend="both",pretty=1,dobh=0,ax=ax)
        ax.set_xlabel(r"$R\ [r_g]$",fontsize=fntsize,labelpad=0)
        ax.set_ylabel(r"$z\ [r_g]$",fontsize=fntsize,labelpad=-10)
        cb.ax.set_xlabel(r"$\log\rho$",fontsize=fntsize,ha="left")
        aphi = psicalc()
        plc(aphi,levels=np.linspace(-aphimax*0.7,aphimax*0.7,21),colors="k",linewidths=1,xy=-1,symmx=1,alpha=0.5,dobh=0,ax=ax)
        ax.set_aspect(asp)
        #Ellipse
        el = Ellipse((0,0), 1.1*2*Rin, 1.1*2*Rin, facecolor='grey', alpha=1)
        art=ax.add_artist(el)
        art.set_zorder(0)
        #print(i);
        plt.title(r"$t=$%05g$\,r_g/c$"%np.round(t)); 


def convert_wrapper(**kwargs):
    if len(sys.argv[2:])==2 and sys.argv[2].isdigit() and sys.argv[3].isdigit():
        whichi = int(sys.argv[2])
        whichn = int(sys.argv[3])
    else:
        print( "Usage: %s %s <whichi> <whichn>" % (sys.argv[0], sys.argv[1]) )
        return
    convert_to_single_file(whichi = whichi, whichn = whichn, **kwargs)

def mkmov_wrapper(**kwargs):
    if len(sys.argv[2:])>=2 and sys.argv[2].isdigit() and sys.argv[3].isdigit():
        whichi = int(sys.argv[2])
        whichn = int(sys.argv[3])
        if len(sys.argv) >= 5:
            startn = int(sys.argv[4])
        else:
            startn = 0
        if len(sys.argv) >= 6: 
            endn = int(sys.argv[4])
        else:
            endn = -1
    else:
        print( "Usage: %s %s <whichi> <whichn> <startn> <endn>" % (sys.argv[0], sys.argv[1]) )
        return
    mkmov(whichi = whichi, whichn = whichn, startn = startn, endn = endn, **kwargs)

def mkmov(startn=0,endn=-1,ln=10,whichi=0,whichn=1,**kwargs):
    which = kwargs.pop("which","mkfrm8panel")
    dosavefig = kwargs.pop("dosavefig",1)
    print("Doing %s movie" % which)
    rg("gdump")
    #compute the total magnetic flux at t = 0
    rd("dump000")
    aphi=psicalc()
    aphimax = aphi.max()
    #construct file list
    flist1 = np.sort(glob.glob( os.path.join(dumps_path, "dump[0-9][0-9][0-9]") ) )
    flist2 = np.sort(glob.glob( os.path.join(dumps_path, "dump[0-9][0-9][0-9][0-9]") ) )
    flist1.sort()
    flist2.sort()
    flist = np.concatenate((flist1,flist2))
    if len(flist) == 0:
        flist1 = np.sort(glob.glob( os.path.join(dumps_path, "dump[0-9][0-9][0-9]_0000") ) )
        flist2 = np.sort(glob.glob( os.path.join(dumps_path, "dump[0-9][0-9][0-9][0-9]_0000") ) )
        flist1.sort()
        flist2.sort()
        flist = np.concatenate((flist1,flist2))
    firsttime = 1
    dpi = 135
    for fldname in flist:
        #find the index of the file
        fldindex = np.int(fldname.split("_")[0].split("p")[-1])
        if fldindex < startn:
            continue
        if endn>=0 and fldindex >= endn:
            break
        if fldindex % whichn != whichi:
            #do every whichn'th snapshot starting with whichi'th snapshot
            continue
        if dosavefig:
            fname = "%s%04d.png" % (which,fldindex)
            if os.path.isfile( fname ):
                print("File %s exists, skipping..." % fname)
                continue
        #print( "Reading " + fldname + " ..." )
        rd("dump%03d" % fldindex);
        if which == "mkfrm8panel":
            if firsttime:
                firsttime = 0
                fig = plt.figure(figsize=(8.*16./9.,8))
                plt.clf()
            mkfrm8panel(fig=fig,aphimax = aphimax)
        elif which == "mkfrm8p3panel":
            if firsttime:
                firsttime = 0
                fig = plt.figure(figsize=(14.4,10.8))
                dpi = 100
                plt.clf()
            mkfrm8p3panel(fig=fig,aphimax = aphimax)
        elif which == "mkfrm8p3panelnocond":
            if firsttime:
                firsttime = 0
                fig = plt.figure(figsize=(14.4,10.8))
                dpi = 100
                plt.clf()
            mkfrm8p3panel(fig=fig,aphimax = aphimax,docond = 0)
        elif which == "mkfrm4panel":
            if firsttime:
                firsttime = 0
                fig = plt.figure(figsize=(12,8))
                plt.clf()
            mkfrm4panel(fig=fig,ax=None,aphimax = aphimax, lnx=17, lny=10)
        elif which == "mkfrm_simple":
            dpi = 225.
            if firsttime:
                firsttime = 0
                fig = plt.figure(figsize=(4.8*16./9.,4.8))
                plt.clf()
            mkfrm_simple(fig=fig,aphimax = aphimax)
        elif which == "mkfrm_simple_2panel":
            dpi = 225.
            if firsttime:
                firsttime = 0
                fig = plt.figure(1,figsize=(4.8*18./9.,4.8))
                plt.clf()
            mkfrm_simple_2panel(fig=fig,aphimax = aphimax)
        else:
            print("Unknown movie type: %s" % which)
            return
        print(fldindex)
        plt.draw()
        if dosavefig:
            plt.savefig(fname,dpi = dpi)

            
#############
def mkfrm4panel(ln=10,lnx=None,lny=None,aphimax=None,fntsize=20,fig=None,xmax=15,ymax=5,asp=1.):
    global cx, cb
    if lnx is None: lnx = ln
    if lny is None: lny = ln
    if fig is None: fig = plt.gcf();
    if "game4" not in globals(): game4 = 4./3.
    aphi=psicalc()
    #fig = plt.gcf()
    if aphimax is None: aphimax = aphi.max()
    gs = GridSpec(2, 2)
    gs.update(left=0.04, right=0.99, top=0.89, bottom=0.09, wspace=0.06, hspace=0.06)
    #
    #top left panel (0,0)
    #
    ax1 = plt.subplot(gs[0,0])
    cmap = cm.jet
    label = None
    #
    # Density
    #
    vmin = -6; vmax = 0; tcks = np.arange(vmin,vmax+1)[::2]
    cs1 = plc(np.log10(rho),ax=ax1,levels=np.linspace(vmin,vmax,100),
             mirrorx=False,k = nz/2,
             isfilled=1,xy=1,dobh=1,xmax=xmax,ymax=ymax,pretty=1,extend="both")
    plc(aphi,ax=ax1,levels=np.linspace(-aphimax,aphimax,10)[1:-1],colors="k",linewidths=1,xy=-1)
    ax1.set_aspect(asp)
    ax1.set_xlim(-lnx,lnx)
    ax1.set_ylim(-lny,lny)
    #ax1.set_xlabel(r"$x\ [r_g]$",fontsize=20)
    ax1.set_ylabel(r"$z\ [r_g]$",fontsize=20,labelpad=-20)
    cx1,cb1 = mkvertcolorbar(ax1,fig,gap=0.02,width=0.05,vmin=vmin,vmax=vmax,loc="top right",
                   label=label,ticks=tcks,fntsize=fntsize*0.8,cmap=cmap,extend="both")
    label = r"$\log\rho$"
    cb1.ax.set_xlabel(label,fontsize=fntsize,ha="center",x=0.5,labelpad=-fntsize)
    plt.setp( ax1.get_xticklabels(), visible=False )
    #
    # Internal energy
    #
    vmin = -7; vmax = -1; tcks = np.arange(vmin,vmax+1)[::2]
    cs2 = plc(np.log10(ug),ax=ax1,levels=np.linspace(vmin,vmax,100),
             mirrorx=True,k = nz/2,
             isfilled=1,xy=-1,dobh=1,pretty=1,extend="both")
    plc(aphi,ax=ax1,levels=np.linspace(-aphimax,aphimax,10)[1:-1],colors="k",linewidths=1,xy=-1,
        mirrorx=True)
    cx2,cb2 = mkvertcolorbar(ax1,fig,gap=0.02,width=0.05,vmin=vmin,vmax=vmax,loc="top left",
                   label=label,ticks=tcks,fntsize=fntsize*0.8,cmap=cmap,extend="both")
    label = r"$\log u$"
    cb2.ax.set_xlabel(label,fontsize=fntsize,ha="center",x=0.5,labelpad=-fntsize)
    #cb2.ax.yaxis.set_ticks_position('left')
    #cb2.ax.tick_params(axis='x',direction='in',labeltop='on')
    #cb2.ax.xaxis.set_ticks_position('top')
    # for label in cb2.ax.get_xticklabels() + cb2.ax.get_yticklabels():
    #     label.set_fontsize(fntsize)
    #
    #bottom left panel (1,0)
    #
    ax1 = plt.subplot(gs[1,0])
    cmap = cm.jet
    label = None
    #
    # Density
    #
    vmin = -6; vmax = 0; tcks = np.arange(vmin,vmax+1)
    cs1 = plc(np.log10(rho[:,ny/2,:]),ax=ax1,levels=np.linspace(vmin,vmax,100),
             mirrorx=False,
             xcoord=(r*cos(ph-pi/2))[:,ny/2,:],ycoord=(r*sin(ph-pi/2))[:,ny/2,:],
             isfilled=1,xy=-1,dobh=1,pretty=1,extend="both")
    ax1.set_aspect(asp)
    ax1.set_xlim(-lnx,lnx)
    ax1.set_ylim(-lny,lny)
    ax1.set_xlabel(r"$x\ [r_g]$",fontsize=20)
    ax1.set_ylabel(r"$y\ [r_g]$",fontsize=20,labelpad=-20)
    #disable color bar here since duplicates top one
    # cx1,cb1 = mkvertcolorbar(ax1,fig,gap=0,vmin=vmin,vmax=vmax,
    #                label=label,ticks=tcks,fntsize=fntsize,cmap=cmap,extend="both")
    # label = r"$\log\rho$"
    # cb1.ax.set_xlabel(label,fontsize=fntsize,ha="left",x=0)
    #
    # Internal energy
    #
    vmin = -7; vmax = -1; tcks = np.arange(vmin,vmax+1)[::2]
    cs2 = plc(np.log10(ug)[:,ny/2,:],ax=ax1,levels=np.linspace(vmin,vmax,100),
             xcoord=(r*cos(ph+pi/2))[:,ny/2,:],ycoord=(r*sin(ph+pi/2))[:,ny/2,:],
             mirrorx=False,
             isfilled=1,xy=-1,dobh=1,pretty=1,extend="both")
    ###########################################################################
    #top right panel (0,2)
    ###########################################################################
    v0 = 1
    ax1 = plt.subplot(gs[0,1])
    cmap = cm.jet
    label = None
    #
    # 4Tec/Tg-1
    #
    vmin = -v0; vmax = v0; tcks = np.arange(vmin,vmax+v0,v0)
    cs1 = plc(4*(game4-1)*ugel4/(gam-1)/ug-1,ax=ax1,levels=np.linspace(vmin,vmax,100),
             mirrorx=False,k = nz/2,
             isfilled=1,xy=1,dobh=1,extend="both",xmax=xmax,ymax=ymax,pretty=1)
    plc(aphi,ax=ax1,levels=np.linspace(-aphimax,aphimax,10)[1:-1],colors="k",linewidths=1,xy=-1)
    ax1.set_aspect(asp)
    ax1.set_xlim(-lnx,lnx)
    ax1.set_ylim(-lny,lny)
    #ax1.set_xlabel(r"$x\ [r_g]$",fontsize=20)
    #ax1.set_ylabel(r"$z\ [r_g]$",fontsize=20,labelpad=-20)
    cx1,cb1 = mkvertcolorbar(ax1,fig,gap=0.02,width=0.05,vmin=vmin,vmax=vmax,loc="top right",
                   label=label,ticks=tcks,fntsize=fntsize*0.8,cmap=cmap,extend="both")
    label = r"$4T_{\rm e,c}/T_{\rm p}{-}1$"
    cb1.ax.set_xlabel(label,fontsize=fntsize,ha="center",x=0.5,labelpad=-fntsize)
    plt.setp( ax1.get_xticklabels(), visible=False )
    plt.setp( ax1.get_yticklabels(), visible=False )
    #
    # 4Te/Tg-1
    #
    vmin = -v0; vmax = v0; tcks = np.arange(vmin,vmax+v0,v0)
    cs1 = plc(4*(game4-1)*ugeldis/(gam-1)/ug-1,ax=ax1,levels=np.linspace(vmin,vmax,100),
             mirrorx=True,k = nz/2,
             isfilled=1,xy=-1,dobh=1,extend="both",pretty=1)
    plc(aphi,ax=ax1,levels=np.linspace(-aphimax,aphimax,10)[1:-1],colors="k",linewidths=1,xy=-1,
        mirrorx=True)
    cx2,cb2 = mkvertcolorbar(ax1,fig,gap=0.02,width=0.05,vmin=vmin,vmax=vmax,loc="top left",
                   label=label,ticks=tcks,fntsize=fntsize*0.8,cmap=cmap,extend="both")
    label = r"$4T_{\rm e}/T_{\rm p}{-}1$"
    cb2.ax.set_xlabel(label,fontsize=fntsize,ha="center",x=0.5,labelpad=-fntsize)
    #
    #bottom right panel (1,2)
    #
    ax1 = plt.subplot(gs[1,1])
    cmap = cm.jet
    label = None
    #
    # 4Te/Tg - 1
    #
    vmin = -v0; vmax = v0; tcks = np.arange(vmin,vmax+v0,v0)
    cs1 = plc((4*(game4-1)*ugel4/(gam-1)/ug-1)[:,ny/2,:],
              ax=ax1,levels=np.linspace(vmin,vmax,100),
              mirrorx=False,
              xcoord=(r*cos(ph-pi/2))[:,ny/2,:],ycoord=(r*sin(ph-pi/2))[:,ny/2,:],
              isfilled=1,xy=-1,dobh=1,pretty=1,extend="both")
    ax1.set_aspect(asp)
    ax1.set_xlim(-lnx,lnx)
    ax1.set_ylim(-lny,lny)
    ax1.set_xlabel(r"$x\ [r_g]$",fontsize=20)
    plt.setp( ax1.get_yticklabels(), visible=False )
    #
    # 4Te,c/Tg-1
    #
    vmin = -v0; vmax = v0; tcks = np.arange(vmin,vmax+v0,v0)
    cs2 = plc((4*(game4-1)*ugeldis/(gam-1)/ug-1)[:,ny/2,:],
              ax=ax1,levels=np.linspace(vmin,vmax,100),
              xcoord=(r*cos(ph+pi/2))[:,ny/2,:],ycoord=(r*sin(ph+pi/2))[:,ny/2,:],
              mirrorx=False,
              isfilled=1,xy=-1,dobh=1,pretty=1,extend="both")
    plt.annotate(r"t=%g" % np.round(t),(0,0),xytext=(5,5),xycoords="figure fraction",textcoords="offset points",fontsize=20,rotation=0,rotation_mode="anchor",ha="left",va="bottom",color="black")

############
            
def mkfrm8panel(ln=10,aphimax=None,fntsize=20,fig=None):
    global cx, cb
    if fig is None: fig = plt.gcf();
    if "game4" not in globals(): game4 = 4./3.
    aphi=psicalc()
    #fig = plt.gcf()
    if aphimax is None: aphimax = aphi.max()
    gs = GridSpec(2, 4)
    gs.update(left=0.04, right=0.99, top=0.89, bottom=0.09, wspace=0.06, hspace=0.06)
    #
    #top left panel (0,0)
    #
    ax1 = plt.subplot(gs[0,0])
    cmap = cm.jet
    label = None
    #
    # Density
    #
    vmin = -6; vmax = 0; tcks = np.arange(vmin,vmax+1)[::2]
    cs1 = plc(np.log10(rho),ax=ax1,levels=np.linspace(vmin,vmax,100),
             mirrorx=False,k = nz/2,
             isfilled=1,xy=1,dobh=1,xmax=10,ymax=5,pretty=1,extend="both")
    plc(aphi,ax=ax1,levels=np.linspace(-aphimax,aphimax,10)[1:-1],colors="k",linewidths=1,xy=-1)
    ax1.set_aspect("equal")
    ax1.set_xlim(-ln,ln)
    ax1.set_ylim(-ln,ln)
    #ax1.set_xlabel(r"$x\ [r_g]$",fontsize=20)
    ax1.set_ylabel(r"$z\ [r_g]$",fontsize=20,labelpad=-20)
    cx1,cb1 = mkvertcolorbar(ax1,fig,gap=0.02,width=0.05,vmin=vmin,vmax=vmax,loc="top right",
                   label=label,ticks=tcks,fntsize=fntsize*0.8,cmap=cmap,extend="both")
    label = r"$\log\rho$"
    cb1.ax.set_xlabel(label,fontsize=fntsize,ha="center",x=0.5,labelpad=-fntsize)
    plt.setp( ax1.get_xticklabels(), visible=False )
    #
    # T_g
    #
    vmin = -4; vmax = -1; tcks = np.arange(vmin,vmax+1)
    cs2 = plc(np.log10(0.5*(gam-1)*ug/rho),ax=ax1,levels=np.linspace(vmin,vmax,100),
             mirrorx=True,k = nz/2,
             isfilled=1,xy=-1,dobh=1,pretty=1,extend="both")
    plc(aphi,ax=ax1,levels=np.linspace(-aphimax,aphimax,10)[1:-1],colors="k",linewidths=1,xy=-1,
        mirrorx=True)
    cx2,cb2 = mkvertcolorbar(ax1,fig,gap=0.02,width=0.05,vmin=vmin,vmax=vmax,loc="top left",
                   label=label,ticks=tcks,fntsize=fntsize*0.8,cmap=cmap,extend="both")
    label = r"$\log T_g/2$"
    cb2.ax.set_xlabel(label,fontsize=fntsize,ha="center",x=0.5,labelpad=-fntsize)
    #cb2.ax.yaxis.set_ticks_position('left')
    #cb2.ax.tick_params(axis='x',direction='in',labeltop='on')
    #cb2.ax.xaxis.set_ticks_position('top')
    # for label in cb2.ax.get_xticklabels() + cb2.ax.get_yticklabels():
    #     label.set_fontsize(fntsize)
    #
    #bottom left panel (1,0)
    #
    ax1 = plt.subplot(gs[1,0])
    cmap = cm.jet
    label = None
    #
    # Density
    #
    vmin = -6; vmax = 1; tcks = np.arange(vmin,vmax+1)
    cs1 = plc(np.log10(rho[:,ny/2,:]),ax=ax1,levels=np.linspace(vmin,vmax,100),
             mirrorx=False,
             xcoord=(r*cos(ph-pi/2))[:,ny/2,:],ycoord=(r*sin(ph-pi/2))[:,ny/2,:],
             isfilled=1,xy=-1,dobh=1,pretty=1,extend="both")
    ax1.set_aspect("equal")
    ax1.set_xlim(-ln,ln)
    ax1.set_ylim(-ln,ln)
    ax1.set_xlabel(r"$x\ [r_g]$",fontsize=20)
    ax1.set_ylabel(r"$y\ [r_g]$",fontsize=20,labelpad=-20)
    #disable color bar here since duplicates top one
    # cx1,cb1 = mkvertcolorbar(ax1,fig,gap=0,vmin=vmin,vmax=vmax,
    #                label=label,ticks=tcks,fntsize=fntsize,cmap=cmap,extend="both")
    # label = r"$\log\rho$"
    # cb1.ax.set_xlabel(label,fontsize=fntsize,ha="left",x=0)
    #
    # Internal energy
    #
    vmin = -4; vmax = -1; tcks = np.arange(vmin,vmax+1,1)
    cs2 = plc(np.log10(0.5*(gam-1)*ug/rho)[:,ny/2,:],ax=ax1,levels=np.linspace(vmin,vmax,100),
             xcoord=(r*cos(ph+pi/2))[:,ny/2,:],ycoord=(r*sin(ph+pi/2))[:,ny/2,:],
             mirrorx=False,
             isfilled=1,xy=-1,dobh=1,pretty=1,extend="both")
    # cx2,cb2 = mkvertcolorbar(ax1,fig,gap=0.04,vmin=vmin,vmax=vmax,loc="top left",
    #                label=label,ticks=tcks,fntsize=fntsize,cmap=cmap,extend="both")
    # label = r"$\log u_g$"
    # cb2.ax.set_xlabel(label,fontsize=fntsize,ha="right",x=1.3)
    # cb2.ax.yaxis.set_ticks_position('left')
    # for label in cb2.ax.get_xticklabels() + cb2.ax.get_yticklabels():
    #     label.set_fontsize(fntsize)
    ###########################################################################
    #top middle left panel (0,1)
    ###########################################################################
    ax1 = plt.subplot(gs[0,1])
    cmap = cm.jet
    label = None
    #
    # log Tecond
    #
    vmin = -4; vmax = -1; tcks = np.arange(vmin,vmax+1,1)
    cs1 = plc(np.log10((game4-1)*ugel4/rho),ax=ax1,levels=np.linspace(vmin,vmax,100),
             mirrorx=False,k = nz/2,
             isfilled=1,xy=1,dobh=1,extend="both",xmax=10,ymax=5,pretty=1)
    plc(aphi,ax=ax1,levels=np.linspace(-aphimax,aphimax,10)[1:-1],colors="k",linewidths=1,xy=-1)
    ax1.set_aspect("equal")
    ax1.set_xlim(-ln,ln)
    ax1.set_ylim(-ln,ln)
    #ax1.set_xlabel(r"$x\ [r_g]$",fontsize=20)
    #ax1.set_ylabel(r"$z\ [r_g]$",fontsize=20,labelpad=-20)
    cx1,cb1 = mkvertcolorbar(ax1,fig,gap=0.02,width=0.05,vmin=vmin,vmax=vmax,loc="top right",
                   label=label,ticks=tcks,fntsize=fntsize*0.8,cmap=cmap,extend="both")
    label = r"$\log T_{\rm e,c}$"
    cb1.ax.set_xlabel(label,fontsize=fntsize,ha="center",x=0.5,labelpad=-fntsize)
    plt.setp( ax1.get_xticklabels(), visible=False )
    plt.setp( ax1.get_yticklabels(), visible=False )
    #
    # Te/Tg
    #
    vmin = -4; vmax = -1; tcks = np.arange(vmin,vmax+1,1)
    cs2 = plc(np.log10((game4-1)*ugeldis/rho),ax=ax1,levels=np.linspace(vmin,vmax,100),
             mirrorx=True,k = nz/2,
             isfilled=1,xy=-1,dobh=1,extend="both",pretty=1)
    plc(aphi,ax=ax1,levels=np.linspace(-aphimax,aphimax,10)[1:-1],colors="k",linewidths=1,xy=-1,
        mirrorx=True)
    cx2,cb2 = mkvertcolorbar(ax1,fig,gap=0.02,width=0.05,vmin=vmin,vmax=vmax,loc="top left",
                   label=label,ticks=tcks,fntsize=fntsize*0.8,cmap=cmap,extend="both")
    label = r"$\log T_{\rm e}$"
    cb2.ax.set_xlabel(label,fontsize=fntsize,ha="center",x=0.5,labelpad=-fntsize)
    # cb2.ax.yaxis.set_ticks_position('left')
    # for label in cb2.ax.get_xticklabels() + cb2.ax.get_yticklabels():
    #     label.set_fontsize(fntsize)
    #
    #bottom middle right panel (1,1)
    #
    ax1 = plt.subplot(gs[1,1])
    cmap = cm.jet
    label = None
    #
    # Te/Tg
    #
    vmin = -4; vmax = -1; tcks = np.arange(vmin,vmax+1,1)
    cs1 = plc(np.log10((game4-1)*ugel4/rho)[:,ny/2,:],
              ax=ax1,levels=np.linspace(vmin,vmax,100),
              mirrorx=False,
              xcoord=(r*cos(ph-pi/2))[:,ny/2,:],ycoord=(r*sin(ph-pi/2))[:,ny/2,:],
              isfilled=1,xy=-1,dobh=1,pretty=1,extend="both")
    ax1.set_aspect("equal")
    ax1.set_xlim(-ln,ln)
    ax1.set_ylim(-ln,ln)
    ax1.set_xlabel(r"$x\ [r_g]$",fontsize=20)
    #ax1.set_ylabel(r"$y\ [r_g]$",fontsize=20,labelpad=-20)
    plt.setp( ax1.get_yticklabels(), visible=False )
    # cx1,cb1 = mkvertcolorbar(ax1,fig,gap=0,vmin=vmin,vmax=vmax,
    #                label=label,ticks=tcks,fntsize=fntsize,cmap=cmap,extend="both")
    # label = r"$\log\frac{T_{\rm e,c}}{T_g}$"
    # cb1.ax.set_xlabel(label,fontsize=fntsize,ha="left",x=0)
    #
    # Tecond/Tg
    #
    vmin = -4; vmax = -1; tcks = np.arange(vmin,vmax+1,1)
    cs2 = plc(np.log10((game4-1)*ugeldis/rho)[:,ny/2,:],
              ax=ax1,levels=np.linspace(vmin,vmax,100),
              xcoord=(r*cos(ph+pi/2))[:,ny/2,:],ycoord=(r*sin(ph+pi/2))[:,ny/2,:],
              mirrorx=False,
              isfilled=1,xy=-1,dobh=1,pretty=1,extend="both")
    # cx2,cb2 = mkvertcolorbar(ax1,fig,gap=0.04,vmin=vmin,vmax=vmax,isleft=1,
    #                label=label,ticks=tcks,fntsize=fntsize,cmap=cmap,extend="both")
    # label = r"$\log\frac{T_{\rm e}}{T_g}$"
    # cb2.ax.set_xlabel(label,fontsize=fntsize,ha="right",x=1.3)
    # cb2.ax.yaxis.set_ticks_position('left')
    # for label in cb2.ax.get_xticklabels() + cb2.ax.get_yticklabels():
    #     label.set_fontsize(fntsize)
    ###########################################################################
    #top middle right panel (0,2)
    ###########################################################################
    ax1 = plt.subplot(gs[0,2])
    cmap = cm.jet
    label = None
    #
    # Tecond/Tg
    #
    vmin = -1; vmax = 1; tcks = np.arange(vmin,vmax+1,1)
    cs1 = plc(np.log10((game4-1)*ugel4/(gam-1)/ug),ax=ax1,levels=np.linspace(vmin,vmax,100),
             mirrorx=False,k = nz/2,
             isfilled=1,xy=1,dobh=1,extend="both",xmax=10,ymax=5,pretty=1)
    plc(aphi,ax=ax1,levels=np.linspace(-aphimax,aphimax,10)[1:-1],colors="k",linewidths=1,xy=-1)
    ax1.set_aspect("equal")
    ax1.set_xlim(-ln,ln)
    ax1.set_ylim(-ln,ln)
    #ax1.set_xlabel(r"$x\ [r_g]$",fontsize=20)
    #ax1.set_ylabel(r"$z\ [r_g]$",fontsize=20,labelpad=-20)
    cx1,cb1 = mkvertcolorbar(ax1,fig,gap=0.02,width=0.05,vmin=vmin,vmax=vmax,loc="top right",
                   label=label,ticks=tcks,fntsize=fntsize*0.8,cmap=cmap,extend="both")
    label = r"$\log T_{\rm e,c}/T_{\rm g}$"
    cb1.ax.set_xlabel(label,fontsize=fntsize,ha="center",x=0.5,labelpad=-fntsize)
    plt.setp( ax1.get_xticklabels(), visible=False )
    plt.setp( ax1.get_yticklabels(), visible=False )
    #
    # Te/Tg
    #
    vmin = -1; vmax = 1; tcks = np.arange(vmin,vmax+1,1)
    cs2 = plc(np.log10((game4-1)*ugeldis/(gam-1)/ug),ax=ax1,levels=np.linspace(vmin,vmax,100),
             mirrorx=True,k = nz/2,
             isfilled=1,xy=-1,dobh=1,extend="both",pretty=1)
    plc(aphi,ax=ax1,levels=np.linspace(-aphimax,aphimax,10)[1:-1],colors="k",linewidths=1,xy=-1,
        mirrorx=True)
    cx2,cb2 = mkvertcolorbar(ax1,fig,gap=0.02,width=0.05,vmin=vmin,vmax=vmax,loc="top left",
                   label=label,ticks=tcks,fntsize=fntsize*0.8,cmap=cmap,extend="both")
    label = r"$\log T_{\rm e}/T_{\rm g}$"
    cb2.ax.set_xlabel(label,fontsize=fntsize,ha="center",x=0.5,labelpad=-fntsize)
    # cb2.ax.yaxis.set_ticks_position('left')
    # for label in cb2.ax.get_xticklabels() + cb2.ax.get_yticklabels():
    #     label.set_fontsize(fntsize)
    #
    #bottom middle right panel (1,2)
    #
    ax1 = plt.subplot(gs[1,2])
    cmap = cm.jet
    label = None
    #
    # Te/Tg
    #
    vmin = -1; vmax = 1; tcks = np.arange(vmin,vmax+0.5,0.5)
    cs1 = plc(np.log10((game4-1)*ugel4/(gam-1)/ug)[:,ny/2,:],
              ax=ax1,levels=np.linspace(vmin,vmax,100),
              mirrorx=False,
              xcoord=(r*cos(ph-pi/2))[:,ny/2,:],ycoord=(r*sin(ph-pi/2))[:,ny/2,:],
              isfilled=1,xy=-1,dobh=1,pretty=1,extend="both")
    ax1.set_aspect("equal")
    ax1.set_xlim(-ln,ln)
    ax1.set_ylim(-ln,ln)
    ax1.set_xlabel(r"$x\ [r_g]$",fontsize=20)
    #ax1.set_ylabel(r"$y\ [r_g]$",fontsize=20,labelpad=-20)
    plt.setp( ax1.get_yticklabels(), visible=False )
    # cx1,cb1 = mkvertcolorbar(ax1,fig,gap=0,vmin=vmin,vmax=vmax,
    #                label=label,ticks=tcks,fntsize=fntsize,cmap=cmap,extend="both")
    # label = r"$\log\frac{T_{\rm e,c}}{T_g}$"
    # cb1.ax.set_xlabel(label,fontsize=fntsize,ha="left",x=0)
    #
    # Tecond/Tg
    #
    vmin = -1; vmax = 1; tcks = np.arange(vmin,vmax+0.5,0.5)
    cs2 = plc(np.log10((game4-1)*ugeldis/(gam-1)/ug)[:,ny/2,:],
              ax=ax1,levels=np.linspace(vmin,vmax,100),
              xcoord=(r*cos(ph+pi/2))[:,ny/2,:],ycoord=(r*sin(ph+pi/2))[:,ny/2,:],
              mirrorx=False,
              isfilled=1,xy=-1,dobh=1,pretty=1,extend="both")
    # cx2,cb2 = mkvertcolorbar(ax1,fig,gap=0.04,vmin=vmin,vmax=vmax,isleft=1,
    #                label=label,ticks=tcks,fntsize=fntsize,cmap=cmap,extend="both")
    # label = r"$\log\frac{T_{\rm e}}{T_g}$"
    # cb2.ax.set_xlabel(label,fontsize=fntsize,ha="right",x=1.3)
    # cb2.ax.yaxis.set_ticks_position('left')
    # for label in cb2.ax.get_xticklabels() + cb2.ax.get_yticklabels():
    #     label.set_fontsize(fntsize)
    ###########################################################################
    #top right panel (0,2)
    ###########################################################################
    v0 = 1
    ax1 = plt.subplot(gs[0,3])
    cmap = cm.jet
    label = None
    #
    # Tec/Tg-0.25
    #
    vmin = -v0; vmax = v0; tcks = np.arange(vmin,vmax+v0,v0)
    cs1 = plc(4*(game4-1)*ugel4/(gam-1)/ug-1,ax=ax1,levels=np.linspace(vmin,vmax,100),
             mirrorx=False,k = nz/2,
             isfilled=1,xy=1,dobh=1,extend="both",xmax=10,ymax=5,pretty=1)
    plc(aphi,ax=ax1,levels=np.linspace(-aphimax,aphimax,10)[1:-1],colors="k",linewidths=1,xy=-1)
    ax1.set_aspect("equal")
    ax1.set_xlim(-ln,ln)
    ax1.set_ylim(-ln,ln)
    #ax1.set_xlabel(r"$x\ [r_g]$",fontsize=20)
    #ax1.set_ylabel(r"$z\ [r_g]$",fontsize=20,labelpad=-20)
    cx1,cb1 = mkvertcolorbar(ax1,fig,gap=0.02,width=0.05,vmin=vmin,vmax=vmax,loc="top right",
                   label=label,ticks=tcks,fntsize=fntsize*0.8,cmap=cmap,extend="both")
    label = r"$4T_{\rm e,c}/T_{\rm g}{-}1$"
    cb1.ax.set_xlabel(label,fontsize=fntsize,ha="center",x=0.5,labelpad=-fntsize)
    plt.setp( ax1.get_xticklabels(), visible=False )
    plt.setp( ax1.get_yticklabels(), visible=False )
    #
    # Te/Tg
    #
    vmin = -v0; vmax = v0; tcks = np.arange(vmin,vmax+v0,v0)
    cs2 = plc(4*(game4-1)*ugeldis/(gam-1)/ug-1,ax=ax1,levels=np.linspace(vmin,vmax,100),
             mirrorx=True,k = nz/2,
             isfilled=1,xy=-1,dobh=1,extend="both",pretty=1)
    plc(aphi,ax=ax1,levels=np.linspace(-aphimax,aphimax,10)[1:-1],colors="k",linewidths=1,xy=-1,
        mirrorx=True)
    cx2,cb2 = mkvertcolorbar(ax1,fig,gap=0.02,width=0.05,vmin=vmin,vmax=vmax,loc="top left",
                   label=label,ticks=tcks,fntsize=fntsize*0.8,cmap=cmap,extend="both")
    label = r"$4T_{\rm e}/T_{\rm g}{-}1$"
    cb2.ax.set_xlabel(label,fontsize=fntsize,ha="center",x=0.5,labelpad=-fntsize)
    # cb2.ax.yaxis.set_ticks_position('left')
    # for label in cb2.ax.get_xticklabels() + cb2.ax.get_yticklabels():
    #     label.set_fontsize(fntsize)
    #
    #bottom right panel (1,2)
    #
    ax1 = plt.subplot(gs[1,3])
    cmap = cm.jet
    label = None
    #
    # 4Tecond/Tg - 1
    #
    vmin = -v0; vmax = v0; tcks = np.arange(vmin,vmax+v0,v0)
    cs1 = plc((4*(game4-1)*ugel4/(gam-1)/ug-1)[:,ny/2,:],
              ax=ax1,levels=np.linspace(vmin,vmax,100),
              mirrorx=False,
              xcoord=(r*cos(ph-pi/2))[:,ny/2,:],ycoord=(r*sin(ph-pi/2))[:,ny/2,:],
              isfilled=1,xy=-1,dobh=1,pretty=1,extend="both")
    ax1.set_aspect("equal")
    ax1.set_xlim(-ln,ln)
    ax1.set_ylim(-ln,ln)
    ax1.set_xlabel(r"$x\ [r_g]$",fontsize=20)
    plt.setp( ax1.get_yticklabels(), visible=False )
    #ax1.set_ylabel(r"$y\ [r_g]$",fontsize=20,labelpad=-20)
    # cx1,cb1 = mkvertcolorbar(ax1,fig,gap=0,vmin=vmin,vmax=vmax,
    #                label=label,ticks=tcks,fntsize=fntsize,cmap=cmap,extend="both")
    # label = r"$\log\frac{T_{\rm e,c}}{T_g}$"
    # cb1.ax.set_xlabel(label,fontsize=fntsize,ha="left",x=0)
    #
    # Te/Tg
    #
    vmin = -v0; vmax = v0; tcks = np.arange(vmin,vmax+v0,v0)
    cs2 = plc((4*(game4-1)*ugeldis/(gam-1)/ug-1)[:,ny/2,:],
              ax=ax1,levels=np.linspace(vmin,vmax,100),
              xcoord=(r*cos(ph+pi/2))[:,ny/2,:],ycoord=(r*sin(ph+pi/2))[:,ny/2,:],
              mirrorx=False,
              isfilled=1,xy=-1,dobh=1,pretty=1,extend="both")
    # cx2,cb2 = mkvertcolorbar(ax1,fig,gap=0.04,vmin=vmin,vmax=vmax,isleft=1,
    #                label=label,ticks=tcks,fntsize=fntsize,cmap=cmap,extend="both")
    # label = r"$\log\frac{T_{\rm e}}{T_g}$"
    # cb2.ax.set_xlabel(label,fontsize=fntsize,ha="right",x=1.3)
    # cb2.ax.yaxis.set_ticks_position('left')
    # for label in cb2.ax.get_xticklabels() + cb2.ax.get_yticklabels():
    #     label.set_fontsize(fntsize)
    plt.annotate(r"t=%g" % np.round(t),(0,0),xytext=(5,5),xycoords="figure fraction",textcoords="offset points",fontsize=20,rotation=0,rotation_mode="anchor",ha="left",va="bottom",color="black")


############
            
def mkfrm8p3panel(aphimax=None,fntsize=20,fig=None,**kwargs):
    global cx, cb
    labelpad = -25
    xlabelpad = -fntsize
    alpha = 1
    bbox = dict(boxstyle="round,pad=0.1", fc="w", ec="w", alpha=1.)
    iti = kwargs.get("iti",None)
    itf = kwargs.get("itf",None)
    fti = kwargs.get("fti",None)
    ftf = kwargs.get("ftf",None)
    dic = get_titf(iti,itf,fti,ftf)
    #copy extra parameters over
    for key in dic: kwargs[key] = dic[key]
    aphimax = kwargs.get("maxaphi",None)
    ncont = kwargs.get("ncont",10)
    ln = kwargs.get("plotlen",10)
    docond = kwargs.get("docond",1)
    cbw = 0.03
    gap = 0.01
    SMALL = 1e-20
    if fig is None: fig = plt.gcf();
    if "game4" not in globals(): game4 = 4./3.
    aphi=psicalc()
    if aphimax is None: aphimax = aphi.max()
    #fig = plt.gcf()
    gs = GridSpec(2, 4)
    left=0.04; right=0.99; top=0.93; bottom=0.32; wspace=0.06; hspace=0.06
    gs.update(left=left, right=right, top=top, bottom=bottom, wspace=wspace, hspace=hspace)
    #
    #top left panel (0,0)
    #
    ax1 = plt.subplot(gs[0,0])
    cmap = cm.jet
    label = None
    #
    # Density
    #
    vmin = kwargs["vmin"]; vmax = kwargs["vmax"]; tcks = np.arange(vmin,vmax+1)[::2]
    cs1 = plc(np.log10(rho),symmx=1,ax=ax1,levels=np.linspace(vmin,vmax,100),
             mirrorx=False,
             isfilled=1,xy=1,dobh=1,xmax=10,ymax=5,pretty=1,extend="both")
    plc(aphi,symmx=1,ax=ax1,levels=np.linspace(-aphimax,aphimax,ncont)[1:-1],colors="k",linewidths=1,xy=-1)
    ax1.set_aspect("equal")
    ax1.set_xlim(-ln,ln)
    ax1.set_ylim(-ln,ln)
    #ax1.set_xlabel(r"$x\ [r_g]$",fontsize=fntsize,labelpad = xlabelpad,alpha=alpha,bbox=bbox,zorder=20)
    ax1.set_ylabel(r"$z\ [r_g]$",fontsize=fntsize,labelpad=labelpad)
    cx1,cb1 = mkvertcolorbar(ax1,fig,gap=gap,width=cbw,vmin=vmin,vmax=vmax,loc="top center",
                   label=label,ticks=tcks,fntsize=fntsize*0.8,cmap=cmap,extend="both")
    label = r"$\log\rho$"
    cb1.ax.set_xlabel(label,fontsize=fntsize,ha="center",x=0.5,labelpad=-1.05*fntsize)
    plt.setp( ax1.get_xticklabels(), visible=False )
    mathify_axes_ticks(ax1)
    #
    #bottom left panel (1,0)
    #
    ax1 = plt.subplot(gs[1,0])
    cmap = cm.jet
    label = None
    #
    # Density
    #
    cs1 = plc(np.log10(rho),xy=-2,symmx=1,ax=ax1,levels=np.linspace(vmin,vmax,100),
             isfilled=1,dobh=1,pretty=1,extend="both")
    ax1.set_aspect("equal")
    ax1.set_xlim(-ln,ln)
    ax1.set_ylim(-ln,ln)
    ax1.set_xlabel(r"$x\ [r_g]$",fontsize=fntsize,labelpad = xlabelpad,alpha=alpha,bbox=bbox,zorder=20)
    ax1.set_ylabel(r"$y\ [r_g]$",fontsize=fntsize,labelpad=labelpad)
    mathify_axes_ticks(ax1)
    ###########################################################################
    #top middle left panel (0,1)
    ###########################################################################
    ax1 = plt.subplot(gs[0,1])
    cmap = cm.jet
    label = None
    #
    # log beta
    #
    vmin = -1; vmax = 3; tcks = np.arange(vmin,vmax+1,1)
    beta = 2*(gam-1)*ug/(bsq+SMALL)
    cs1 = plc(np.log10(beta),ax=ax1,levels=np.linspace(vmin,vmax,100),
             symmx=1,
             isfilled=1,xy=1,dobh=1,extend="both",pretty=1)
    plc(aphi,symmx=1,ax=ax1,levels=np.linspace(-aphimax,aphimax,ncont)[1:-1],colors="k",linewidths=1,xy=-1)
    ax1.set_aspect("equal")
    ax1.set_xlim(-ln,ln)
    ax1.set_ylim(-ln,ln)
    cx1,cb1 = mkvertcolorbar(ax1,fig,gap=gap,width=cbw,vmin=vmin,vmax=vmax,loc="top center",
                   label=label,ticks=tcks,fntsize=fntsize*0.8,cmap=cmap,extend="both")
    label = r"$\log \beta$"
    cb1.ax.set_xlabel(label,fontsize=fntsize,ha="center",x=0.5,labelpad=-1.05*fntsize)
    plt.setp( ax1.get_xticklabels(), visible=False )
    plt.setp( ax1.get_yticklabels(), visible=False )
    mathify_axes_ticks(ax1)
    #
    #bottom middle right panel (1,1)
    #
    ax1 = plt.subplot(gs[1,1])
    cmap = cm.jet
    label = None
    cs1 = plc(np.log10(beta),
              ax=ax1,levels=np.linspace(vmin,vmax,100),
              symmx=1,
              isfilled=1,xy=-2,dobh=1,pretty=1,extend="both")
    ax1.set_aspect("equal")
    ax1.set_xlim(-ln,ln)
    ax1.set_ylim(-ln,ln)
    ax1.set_xlabel(r"$x\ [r_g]$",fontsize=fntsize,labelpad = xlabelpad,alpha=alpha,bbox=bbox,zorder=20)
    #ax1.set_ylabel(r"$y\ [r_g]$",fontsize=fntsize,labelpad=labelpad)
    plt.setp( ax1.get_yticklabels(), visible=False )
    mathify_axes_ticks(ax1)
    ###########################################################################
    #top middle right panel (0,2)
    ###########################################################################
    v0 = 1
    ax1 = plt.subplot(gs[0,2])
    cmap = cm.jet
    label = None
    #
    # f_e
    #
    if not np.isnan(phi.max()) and docond:
        fe = felcalc(ugel4)+SMALL
        label = r"$\log f_{\rm e}$"
    else:
        fe = felcalc(ugeldis)+SMALL
        label = r"$\log f_{\rm e}\ ({\rm no\ cond})$"
    vmin = -3; vmax = 0; tcks = np.arange(vmin,vmax+v0,v0)
    cs1 = plc(np.log10(fe),ax=ax1,levels=np.linspace(vmin,vmax,100),
             symmx = 1,
             isfilled=1,xy=1,dobh=1,extend="both",xmax=10,ymax=5,pretty=1)
    plc(aphi,symmx=1,ax=ax1,levels=np.linspace(-aphimax,aphimax,ncont)[1:-1],colors="k",linewidths=1,xy=-1)
    ax1.set_aspect("equal")
    ax1.set_xlim(-ln,ln)
    ax1.set_ylim(-ln,ln)
    cx1,cb1 = mkvertcolorbar(ax1,fig,gap=gap,width=cbw,vmin=vmin,vmax=vmax,loc="top center",
                   label=label,ticks=tcks,fntsize=fntsize*0.8,cmap=cmap,extend="both")
    cb1.ax.set_xlabel(label,fontsize=fntsize,ha="center",x=0.5,labelpad=-1.05*fntsize)
    plt.setp( ax1.get_xticklabels(), visible=False )
    plt.setp( ax1.get_yticklabels(), visible=False )
    mathify_axes_ticks(ax1)
    #
    #bottom middle right panel (1,2)
    #
    ax1 = plt.subplot(gs[1,2])
    cmap = cm.jet
    label = None
    cs1 = plc(np.log10(fe),
              ax=ax1,levels=np.linspace(vmin,vmax,100),
              symmx=1,
              isfilled=1,xy=-2,dobh=1,pretty=1,extend="both")
    ax1.set_aspect("equal")
    ax1.set_xlim(-ln,ln)
    ax1.set_ylim(-ln,ln)
    ax1.set_xlabel(r"$x\ [r_g]$",fontsize=fntsize,labelpad = xlabelpad,alpha=alpha,bbox=bbox,zorder=20)
    plt.setp( ax1.get_yticklabels(), visible=False )
    mathify_axes_ticks(ax1)
    ###########################################################################
    #top right panel (0,3)
    ###########################################################################
    ax1 = plt.subplot(gs[0,3])
    cmap = cm.jet
    label = None
    #
    # Theta_e
    #
    mrat = 1836.152672
    if not np.isnan(phi.max()) and docond:
        Thetae = (game4-1)*ugel4/rho*mrat
        label = r"$\log \Theta_{\rm e}$"
    else:
        Thetae = (game4-1)*ugeldis/rho*mrat
        label = r"$\log \Theta_{\rm e}\ ({\rm no\ cond})$"
    thetaemin = kwargs.get("thetaemin",0)
    thetaemax = kwargs.get("thetaemax",2)
    vmin = thetaemin; vmax = thetaemax; tcks = np.arange(vmin,vmax+1,1)
    cs1 = plc(np.log10(Thetae),ax=ax1,levels=np.linspace(vmin,vmax,100),
             symmx = 1,
             isfilled=1,xy=1,dobh=1,extend="both",pretty=1)
    plc(aphi,symmx=1,ax=ax1,levels=np.linspace(-aphimax,aphimax,ncont)[1:-1],colors="k",linewidths=1,xy=-1)
    ax1.set_aspect("equal")
    ax1.set_xlim(-ln,ln)
    ax1.set_ylim(-ln,ln)
    cx1,cb1 = mkvertcolorbar(ax1,fig,gap=gap,width=cbw,vmin=vmin,vmax=vmax,loc="top center",
                   label=label,ticks=tcks,fntsize=fntsize*0.8,cmap=cmap,extend="both")
    cb1.ax.set_xlabel(label,fontsize=fntsize,ha="center",x=0.5,labelpad=-1.05*fntsize)
    plt.setp( ax1.get_xticklabels(), visible=False )
    plt.setp( ax1.get_yticklabels(), visible=False )
    mathify_axes_ticks(ax1)
    #
    #bottom right panel (1,3)
    #
    ax1 = plt.subplot(gs[1,3])
    cmap = cm.jet
    label = None
    cs1 = plc(np.log10(Thetae),
              ax=ax1,levels=np.linspace(vmin,vmax,100),
              symmx=1,
              isfilled=1,xy=-2,dobh=1,pretty=1,extend="both")
    ax1.set_aspect("equal")
    ax1.set_xlim(-ln,ln)
    ax1.set_ylim(-ln,ln)
    ax1.set_xlabel(r"$x\ [r_g]$",fontsize=fntsize,labelpad = xlabelpad,alpha=alpha,bbox=bbox,zorder=20)
    plt.setp( ax1.get_yticklabels(), visible=False )
    mathify_axes_ticks(ax1)
    #
    # bottom time-dependence panels
    #
    topts = 0.28
    v = np.load("qty.npz")
    #choose large labelpad to avoid the x-label
    plot_ts_panels(v,left=0.06,right=0.96,top=topts,bottom=0.03,xlabelpad=100,fntsize=fntsize,doannotate=0,**kwargs)
    plt.annotate(r"$t=%g\ r_g/c$" % np.round(t),(0.5,topts),xytext=(-40,-2),xycoords="figure fraction",textcoords="offset points",fontsize=fntsize,rotation=0,rotation_mode="anchor",ha="left",va="top",color="black",zorder=20)
    v.close()


def mkvertcolorbar(ax,fig,vmin=0,vmax=1,label=None,ylabel=None,ticks=None,fntsize=20,cmap=mpl.cm.jet,gap=0.03,width=0.02,extend="neither",loc="right"):
    box = ax.get_position()
    #pdb.set_trace()
    # cpos = [box.x0,box.y0+box.height+0.05,box.width,0.03]
    locs = loc.split()
    loc0 = locs[0]
    if len(locs)>1:
        loc1 = locs[1]
    else:
        loc1 = None
    if loc0 == "left":
        cpos = box.x0-gap-width,box.y0,width,box.height
    elif loc0 == "right":
        cpos = box.x0+box.width+gap,box.y0,width,box.height
    elif loc0 == "top":
        if loc1 == "right":
            cpos = box.x0+box.width*0.55,box.y0+box.height+gap,box.width*0.45,width
        elif loc1 == "left":
            cpos = box.x0+box.width*0.0,box.y0+box.height+gap,box.width*0.45,width
        else:
            cpos = box.x0,box.y0+box.height+gap,box.width,width
    ax1 = fig.add_axes(cpos)
    #cmap = mpl.cm.jet
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    if loc0 == "left" or loc0 == "right":
        ori = "vertical"
    else:
        ori = "horizontal"
    if ticks is not None:
        cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
                                        norm=norm,
                                        orientation=ori,
                                        ticks=ticks,
                                        extend=extend)
    else:
        cb1 = mpl.colorbar.ColorbarBase(ax1,  cmap=cmap,
                                        norm=norm,
                                        orientation=ori,
                                        extend=extend)
    if loc0 == "top":
        cb1.ax.xaxis.set_ticks_position('top')
        mathify_axes_ticks(cb1.ax,fontsize=fntsize,xticks=ticks)
    elif loc0 == "left":
        cb1.ax.yaxis.set_ticks_position('left')
        mathify_axes_ticks(cb1.ax,fontsize=fntsize,yticks=ticks)
    if label is not None:
        ax1.set_xlabel(label,fontsize=fntsize)
    if ylabel is not None:
        ax1.set_ylabel(ylabel,fontsize=fntsize)
    for label in ax1.get_xticklabels() + ax1.get_yticklabels():
        label.set_fontsize(fntsize)
    return ax1,cb1
            
def compute_udrift():
    #set velocity to drift velocity
    SMALL = 1.0e-20
    GAMMAMAX=250.0
    beta = -bu[0]/((bsq+SMALL)*uu[0])
    betasq = beta*beta*bsq
    betasqmax = 1.-1./(GAMMAMAX*GAMMAMAX)
    betasq[betasq>betasqmax] = betasqmax + 0*betasq[betasq>betasqmax]
    gamma = 1./sqrt(1.-betasq)
    uudr = gamma*(uu+beta*bu)
    return uudr

def Qmri(dir=2):
    """
    APPROXIMATELY Computes number of theta cells resolving one MRI wavelength
    """
    global bu,rho,uu,_dx2,_dx3
    #cvel()
    #corrected this expression to include both 2pi and dxdxp[3][3]
    #also corrected defition of va^2 to contain bsq+gam*ug term
    #need to figure out how to properly measure this in fluid frame
    vaudir = np.abs(bu[dir])/np.sqrt(rho+bsq+gam*ug)
    omega = dxdxp[3][3]*uu[3]/uu[0]+1e-15
    lambdamriudir = 2*np.pi * vaudir / omega
    if dir == 2:
        res=lambdamriudir/_dx2
    elif dir == 3:
        res=lambdamriudir/_dx3
    return(res)

def goodlabs(fntsize=20):
    ax = plt.gca()
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(fntsize)


def plotdth(rx=5,sf=0,hor=0.3):
    ix = iofr(20)
    x = (h[ix,:,0]-pi/2)/hor
    y = (dxdxp[2,2]*_dx2)[ix,:,0]/hor
    plt.clf()
    plt.plot(x,y,"k",lw=2)
    ax1 = plt.gca()
    goodlabs()
    plt.grid()
    plt.xlabel(r"$(\theta-\pi/2)/(h/r)$",fontsize=20)
    plt.ylabel(r"${\rm d}\theta/(h/r)$",fontsize=20)
    plt.twinx()
    ax2 = plt.gca()
    ix = iofr(rx)
    l,=plt.plot(x,0.1*rho.mean(-1)[ix],"r",lw=2)
    l.set_dashes([10,5])
    plt.xlim(0,x.max())
    ax2.set_ylabel(r"$\rho$",fontsize=20)
    #plt.draw()
    ax1ticks = ax1.get_yticks()
    nticks = len(ax1ticks)
    ax2ticks = ax2.get_yticks()
    ax2ymax = ax2ticks[-1]
    ax2ticks = np.linspace(0,ax2ymax,nticks)
    ax1.set_yticks(ax1ticks)
    ax2.set_yticks(ax2ticks)
    goodlabs()
    if sf: plt.savefig("dthvsth.pdf",bbox_inches='tight',pad_inches=0.02)


def iofr(rval):
    rval = np.array(rval)
    if np.max(rval) < r[0,0,0]:
        return 0
    res = interp1d(r[:,0,0], ti[:,0,0], kind='linear', bounds_error = False, fill_value = 0)(rval)
    if len(res.shape)>0 and len(res)>0:
        res[rval<r[0,0,0]]*=0
        res[rval>r[nx-1,0,0]]=res[rval>r[nx-1,0,0]]*0+nx-1
    else:
        res = np.float64(res)
    return(np.floor(res+0.5).astype(int))


def horcalc1d(which=1):
    """
    Compute root mean square deviation of disk body from equatorial plane
    """
    tiny=np.finfo(rho.dtype).tiny
    up=(gdet*rho*(h-np.pi/2)*which).sum(axis=1)
    dn=(gdet*rho*which).sum(axis=1)
    thetamid2d=up/(dn+tiny)+np.pi/2
    thetamid3d=np.empty((nx,ny,nz),dtype=h.dtype)
    hoverr3d=np.empty((nx,ny,nz),dtype=h.dtype)
    up=(gdet*rho*(h-thetamid2d[:,None,:])**2*which).sum(-1).sum(-1)
    dn=(gdet*rho*which).sum(-1).sum(-1)
    hoverr1d= (up/(dn+tiny))**0.5
    return(hoverr1d)

def plotTel():
    #internal energy ratio map
    res,cb = plco(ugel4.mean(-1)/ug.mean(-1),cb=1,cmap=cm.jet,levels=np.linspace(0.1,2.1,100),isfilled=1,nc=100,xy=1,xmax=24.95,ymax=12.5,cbxla=r"$u_{g,e}/u_g$",xla=r"$x\ [r_g]$",yla="$y\ [r_g]$",extend="both")
    cb.set_ticks(np.arange(0.1,2.3,0.2))
    plt.gca().set_aspect("equal")
    #grid (every 16th gridline)
    plc(ti,levels=np.arange(0,nx+16,16)[1:],xy=-1,colors="w")
    plc(tj,levels=np.arange(0,ny,16),xy=-1,colors="w")
    #compute disk thickness
    hor = horcalc1d()
    #show disk h/r with red lines
    plc(abs(np.pi/2-h[:,:,0])-hor[:,None],levels=(0,),xy=-1,linewidths=2,colors="red")
    plt.savefig("ugeoug_fe05_wgridhor.png",dpi=300,bbox_inches='tight',pad_inches=0.02)

#read in a dump file
def rd(dump):
    read_file(dump,type="dump")


#read in a grid file
def rg(dump,**kwargs):
    read_file(dump,type="gdump",**kwargs)

#high-level function that reads either MPI or serial gdump's
def read_file(dump,type=None,savedump=True,saverdump=False):
    if type is None:
        if dump.startswith("dump"):
            type = "dump"
            print("Reading a dump file %s ..." % dump)
        elif dump.startswith("gdump"):
            type = "gdump"
            print("Reading a gdump file %s ..." % dump)
        elif dump.startswith("rdump"):
            type = "rdump"
            print("Reading a rdump file %s ..." % dump)
        elif dump.startswith("fdump"):
            type = "fdump"
            print("Reading a fdump file %s ..." % dump)
        else:
            print("Couldn't guess dump type; assuming it is a data dump")
            type = "dump"
    #normal dump
    if os.path.isfile( dumps_path + dump ):
        headerline = read_header(dumps_path + dump, returnheaderline = True)
        gd = read_body(dumps_path + dump,nx=N1+2*N1G,ny=N2+2*N2G,nz=N3+2*N3G,noround=1)
        if saverdump:
            #if the full dump file does not exist, create it
            dumpfullname = dumps_path + dump + "flr"
            if type == "rdump" and not os.path.isfile(dumpfullname):
                sys.stdout.write("Saving rdump to %s..." % dumpfullname)
                sys.stdout.flush()
                fout = open( dumpfullname, "wb" )
                #join header items with " " (space) as a glue
                #see http://stackoverflow.com/questions/12377473/python-write-versus-writelines-and-concatenated-strings
                #write it out with a new line char at the end
                fout.write(headerline)
                fout.flush()
                os.fsync(fout.fileno())
                #reshape the dump content
                gd1 = gd.transpose(1,2,3,0)
                #at this point gd1 is exactly the unmodified rdump
                #now, modify it by adding flr at its end and inserting kel4x in the middle
                #initialize kel4a-kel4e with keldis
                kel4x = gd1[...,0:5]*0 + gd1[...,11:12]
                flr = gd1[...,0:1]*0
                gd1 = np.concatenate((gd1[...,:10],kel4x,gd1[...,10:],flr),axis=-1)
                gd1.tofile(fout)
                fout.close()
                print( " done!" )
        res = data_assign(myfloat(gd),type=type,nx=N1+2*N1G,ny=N2+2*N2G,nz=N3+2*N3G)
        return res
    #MPI-type dump that is spread over many files
    else:
        flist = np.sort(glob.glob( dumps_path + dump + "_[0-9][0-9][0-9][0-9]" ))
        if len(flist) == 0:
            print( "Could not find %s or its MPI counterpart" % dump )
            return
        sys.stdout.write( "Reading %s (%d files)" % (dump, len(flist)) )
        sys.stdout.flush()
        ndots = 10
        dndot = len(flist)/ndots
        if dndot == 0: dndot = 1
        for i,fname in enumerate(flist):
            #print( "Reading file %d out of %d..." % (i,len(flist)) )
            #header for each file might be different, so read each
            header = read_header(fname,issilent=1)
            #headerline = read_header(fname, issilent=1,returnheaderline = True)
            #header = headerline.split()
            if header is None:
                print( "Error reading header of %s, aborting..." % fname )
                return
            lgd = read_body(fname,nx=N1+2*N1G,ny=N2+2*N2G,nz=N3+2*N3G)
            #this gives an array of dimensions (-1,N1,N2,N3)+potentially ghost cells
            if 0 == i:
                #create full array: of dimensions, (-1,nx,ny,nz)
                fgd = np.zeros( (lgd.shape[0], nx+2*N1G, ny+2*N2G, nz+2*N3G), dtype=np.float32)
            if not type == "rdump" and not type == "fdump":
                #construct full indices: ti, tj, tk
                #fti,ftj,ftk = mgrid[0:nx,0:ny,0:nz]
                lti,ltj,ltk = lgd[0:3,:,:].view();
                lti = np.int64(lti)
                ltj = np.int64(ltj)
                ltk = np.int64(ltk)
                fgd[:,lti+N1G,ltj+N2G,ltk+N3G] = lgd[:,:,:,:]
            else:
                print(starti,startj,startk)
                fgd[:,starti:starti+N1+2*N1G,startj:startj+N2+2*N2G,startk:startk+N3+2*N3G] = lgd[:,:,:,:]
            del lgd
            if i%dndot == 0:
                sys.stdout.write(".")
                sys.stdout.flush()
        res = data_assign(fgd,type=type,nx=nx+2*N1G,ny=ny+2*N2G,nz=nz+2*N3G)
        if savedump:
            #if the full dump file does not exist, create it
            dumpfullname = dumps_path + dump
            if (type == "dump" or type == "gdump") and not os.path.isfile(dumpfullname):
                sys.stdout.write("Saving full dump to %s..." % dumpfullname)
                sys.stdout.flush()
                header[1] = header[4] #N1 = nx
                header[2] = header[5] #N2 = ny
                header[3] = header[6] #N3 = nz
                fout = open( dumpfullname, "wb" )
                #join header items with " " (space) as a glue
                #see http://stackoverflow.com/questions/12377473/python-write-versus-writelines-and-concatenated-strings
                #write it out with a new line char at the end
                hline = ""
                for ii in range(len(header)):
                    hline += header[ii].decode("utf-8") + " "
                hline += "\n"
                headerline = bytes(hline, 'utf-8')
                fout.write(headerline)
                #fout.write(" ".join(header) + "\n")
                fout.flush()
                os.fsync(fout.fileno())
                #reshape the dump content
                gd1 = fgd.transpose(1,2,3,0)
                gd1.tofile(fout)
                fout.close()
                print( " done!" )
                if res is not None:
                    return res
        return res

#read in a header
def read_header(dump,issilent=True,returnheaderline=False):
    global t,nx,ny,nz,N1,N2,N3,N1G,N2G,N3G,starti,startj,startk,_dx1,_dx2,_dx3,a,gam,Rin,Rout,hslope,R0,ti,tj,tk,x1,x2,x3,r,h,ph,gcov,gcon,gdet,drdx,gn3,gv3,guu,gdd,dxdxp, games, startx1, startx2, startx3, x10, x20, game, game4, game5, tf, NPR, DOKTOT, eHEAT, eCOND, DONUCLEAR, DOFLR, DOCYLINDRIFYCOORDS
    global fractheta
    global fracphi
    global rbr
    global npow2
    global cpow2
    global global_x10
    global global_x20
    global global_fracdisk
    global global_fracjet
    global global_r0disk
    global global_rdiskend
    global global_r0jet
    global global_rjetend
    global global_jetnu
    global global_rsjet
    global global_r0grid
    global BL
    #read image
    fin = open( dump, "rb" )
    headerline = fin.readline()
    header = headerline.split()
    nheadertot = len(header)
    fin.close()
    if not dump.startswith(dumps_path + "rdump"):
        if not issilent: print( "dump header: len(header) = %d" % len(header) )
        nheader = 37
        n = 0
        t = myfloat(np.float64(header[n])); n+=1
        #per tile resolution
        N1 = int(header[n]); n+=1
        N2 = int(header[n]); n+=1
        N3 = int(header[n]); n+=1
        #total resolution
        nx = int(header[n]); n+=1
        ny = int(header[n]); n+=1
        nz = int(header[n]); n+=1
        #numbers of ghost cells
        N1G = int(header[n]); n+=1
        N2G = int(header[n]); n+=1
        N3G = int(header[n]); n+=1
        startx1 = myfloat(float(header[n])); n+=1
        startx2 = myfloat(float(header[n])); n+=1
        startx3 = myfloat(float(header[n])); n+=1
        _dx1=myfloat(float(header[n])); n+=1
        _dx2=myfloat(float(header[n])); n+=1
        _dx3=myfloat(float(header[n])); n+=1
        tf=myfloat(float(header[n])); n+=1
        nstep=myfloat(float(header[n])); n+=1
        a=myfloat(float(header[n])); n+=1
        gam=myfloat(float(header[n])); n+=1
        cour=myfloat(float(header[n])); n+=1
        DTd=myfloat(float(header[n])); n+=1
        DTl=myfloat(float(header[n])); n+=1
        DTi=myfloat(float(header[n])); n+=1
        DTr=myfloat(float(header[n])); n+=1
        DTr01=myfloat(float(header[n])); n+=1
        dump_cnt=myfloat(float(header[n])); n+=1
        image_cnt=myfloat(float(header[n])); n+=1
        rdump_cnt=myfloat(float(header[n])); n+=1
        rdump01_cnt=myfloat(float(header[n])); n+=1
        dt=myfloat(float(header[n])); n+=1
        lim=myfloat(float(header[n])); n+=1
        failed=myfloat(float(header[n])); n+=1
        Rin=myfloat(float(header[n])); n+=1
        Rout=myfloat(float(header[n])); n+=1
        hslope=myfloat(float(header[n])); n+=1
        R0=myfloat(float(header[n])); n+=1
        if n<len(header):
            NPR=int(header[n]); n+=1
            DOKTOT=int(header[n]); n+=1
            eHEAT=int(header[n]); n+=1
            eCOND=int(header[n]); n+=1
            DONUCLEAR=int(header[n]); n+=1
            nheader = 42
        else:
            NPR=-1
            DOKTOT=-1
            eHEAT=-1
            eCOND=-1
            DONUCLEAR=0
            DOFLR = 0
        if n<len(header):
          DOFLR=int(header[n]); n+=1
          nheader = 43
        else:
          DOFLR = 0
        if n<len(header):
            nheader = 60
            DOCYLINDRIFYCOORDS = myfloat(header[n]); n+=1
            fractheta = myfloat(header[n]); n+=1
            fracphi   = myfloat(header[n]); n+=1
            rbr       = myfloat(header[n]); n+=1
            npow2     = myfloat(header[n]); n+=1
            cpow2     = myfloat(header[n]); n+=1
            global_x10 = myfloat(header[n]); n+=1
            global_x20 = myfloat(header[n]); n+=1
            global_fracdisk   = myfloat(header[n]); n+=1
            global_fracjet    = myfloat(header[n]); n+=1
            global_r0disk     = myfloat(header[n]); n+=1
            global_rdiskend   = myfloat(header[n]); n+=1
            global_r0jet      = myfloat(header[n]); n+=1
            global_rjetend    = myfloat(header[n]); n+=1
            global_jetnu      = myfloat(header[n]); n+=1
            global_rsjet      = myfloat(header[n]); n+=1
            global_r0grid     = myfloat(header[n]); n+=1
        if n<len(header):
            nheader = 61
            BL = myfloat(header[n]); n+=1
    else:
        print("rdump header")
        nheader = 51
        n = 0
        #per tile resolution
        N1 = int(header[n]); n+=1
        N2 = int(header[n]); n+=1
        N3 = int(header[n]); n+=1
        #total resolution
        nx = int(header[n]); n+=1
        ny = int(header[n]); n+=1
        nz = int(header[n]); n+=1
        #numbers of ghost cells
        N1G = int(header[n]); n+=1
        N2G = int(header[n]); n+=1
        N3G = int(header[n]); n+=1
        #starting indices
        starti = int(header[n]); n+=1
        startj = int(header[n]); n+=1
        startk = int(header[n]); n+=1
        t = myfloat(header[n]); n+=1
        tf = myfloat(header[n]); n+=1
        nstep = int(header[n]); n+=1
        a = myfloat(header[n]); n+=1
        gam = myfloat(header[n]); n+=1
        game = myfloat(header[n]); n+=1
        game4 = myfloat(header[n]); n+=1
        game5 = myfloat(header[n]); n+=1
        cour = myfloat(header[n]); n+=1
        DTd = myfloat(header[n]); n+=1
        DTl = myfloat(header[n]); n+=1
        DTi = myfloat(header[n]); n+=1
        DTr = myfloat(header[n]); n+=1
        DTr01 = myfloat(header[n]); n+=1
        dump_cnt = myfloat(header[n]); n+=1
        image_cnt = myfloat(header[n]); n+=1
        rdump_cnt = myfloat(header[n]); n+=1
        rdump01_cnt=myfloat(float(header[n])); n+=1
        dt = myfloat(header[n]); n+=1
        lim = myfloat(header[n]); n+=1
        failed = myfloat(header[n]); n+=1
        Rin = myfloat(header[n]); n+=1
        Rout = myfloat(header[n]); n+=1
        hslope = myfloat(header[n]); n+=1
        R0 = myfloat(header[n]); n+=1
        fractheta = myfloat(header[n]); n+=1
        fracphi = myfloat(header[n]); n+=1
        rbr = myfloat(header[n]); n+=1
        npow2 = myfloat(header[n]); n+=1
        cpow2 = myfloat(header[n]); n+=1
        x10 = myfloat(header[n]); n+=1
        x20 = myfloat(header[n]); n+=1
        mrat = myfloat(header[n]); n+=1
        fel0 = myfloat(header[n]); n+=1
        felfloor = myfloat(header[n]); n+=1
        tdump = myfloat(header[n]); n+=1
        trdump = myfloat(header[n]); n+=1
        timage = myfloat(header[n]); n+=1
        tlog  = myfloat(header[n]); n+=1
    if n < len(header):
        nheader = 60
        global_fracdisk   = myfloat(header[n]); n+=1
        global_fracjet    = myfloat(header[n]); n+=1
        global_r0disk     = myfloat(header[n]); n+=1
        global_rdiskend   = myfloat(header[n]); n+=1
        global_r0jet      = myfloat(header[n]); n+=1
        global_rjetend    = myfloat(header[n]); n+=1
        global_jetnu      = myfloat(header[n]); n+=1
        global_rsjet      = myfloat(header[n]); n+=1
        global_r0grid     = myfloat(header[n]); n+=1
    if n != nheader or n != nheadertot:
        print("Wrong number of elements in header: nread = %d, nexpected = %d, nototal = %d: incorrect format?"
              % (n, nheader, nheadertot) )
        return headerline
    if returnheaderline:
        return headerline
    else:
        return header
            
def read_body(dump,nx=None,ny=None,nz=None,noround=False):
        fin = open( dump, "rb" )
        header = fin.readline()
        if dump.startswith(dumps_path + "rdump"):
            dtype = np.float64
            body = np.fromfile(fin,dtype=dtype,count=-1)
            gd = body.view().reshape((nx,ny,nz,-1), order='C')
            if noround:
                gd=gd.transpose(3,0,1,2)
            else:
                gd=myfloat(gd.transpose(3,0,1,2))
        elif dump.startswith(dumps_path + "fdump"):
            dtype = np.int64
            body = np.fromfile(fin,dtype=dtype,count=-1)
            gd = body.view().reshape((-1,nz,ny,nx), order='F')
            gd=myfloat(gd.transpose(0,3,2,1))
        else:
            dtype = np.float32
            body = np.fromfile(fin,dtype=dtype,count=-1)
            gd = body.view().reshape((-1,nz,ny,nx), order='F')
            gd=myfloat(gd.transpose(0,3,2,1))
        return gd

def data_assign(gd,type=None,**kwargs):
    if type is None:
        print("Please specify data type")
        return
    if type == "gdump":
        gdump_assign(gd,**kwargs)
        return None
    elif type == "dump":
        dump_assign(gd,**kwargs)
        return None
    elif type == "rdump":
        gd = rdump_assign(gd,**kwargs)
        return gd
    elif type == "fdump":
        gd = fdump_assign(gd,**kwargs)
        return gd
    else:
        print("Unknown data type: %s" % type)
        return gd
    
def gdump_assign(gd,**kwargs):
    global t,nx,ny,nz,N1,N2,N3,_dx1,_dx2,_dx3,a,gam,Rin,Rout,hslope,R0,ti,tj,tk,x1,x2,x3,r,h,ph,gcov,gcon,gdet,drdx,gn3,gv3,guu,gdd,dxdxp, games
    nx = kwargs.pop("nx",nx)
    ny = kwargs.pop("ny",ny)
    nz = kwargs.pop("nz",nz)
    ti,tj,tk,x1,x2,x3,r,h,ph = gd[0:9,:,:].view();  n = 9
    gv3 = gd[n:n+16].view().reshape((4,4,nx,ny,nz),order='F').transpose(1,0,2,3,4); n+=16
    gn3 = gd[n:n+16].view().reshape((4,4,nx,ny,nz),order='F').transpose(1,0,2,3,4); n+=16
    gcov = gv3
    gcon = gn3
    guu = gn3
    gdd = gv3
    gdet = gd[n]; n+=1
    drdx = gd[n:n+16].view().reshape((4,4,nx,ny,nz),order='F').transpose(1,0,2,3,4); n+=16
    dxdxp = drdx
    if n != gd.shape[0]:
        print("rd: WARNING: nread = %d < ntot = %d: incorrect format?" % (n, gd.shape[0]) )
        return 1
    return 0

#read in a dump file
def dump_assign(gd,**kwargs):
    global t,nx,ny,nz,_dx1,_dx2,_dx3,gam,hslope,a,R0,Rin,Rout,ti,tj,tk,x1,x2,x3,r,h,ph,rho,ug,vu,B,pg,cs2,Sden,U,gdetB,divb,uu,ud,bu,bd,v1m,v1p,v2m,v2p,gdet,bsq,gdet,alpha,rhor, qdot, ktot, Ttot, game, qisosq, pflag, qisodotb, kel, uelvar, Tel4, Tel5,Teldis, Tels, kel4, kel5,ugel,ugeldis, ugcon, sel, ugscon, ugel4, ugel5,stot, uelvar, Telvar, Tsel, sel, ugels, games, phi, keldis, phihat,csphib,lrho, Tnuc, Tnuc_cgs, etae, flr, G, Q
    global kel4a, kel4b, kel4c, kel4d, kel4e, ugel4a, ugel4b, ugel4c, ugel4d, ugel4e
    global Tel4a, Tel4b, Tel4c, Tel4d, Tel4e
    global rhonp, rhoalpha, rhofloor, Ye
    
    global flrfrac

    global pg
    
    nx = kwargs.pop("nx",nx)
    ny = kwargs.pop("ny",ny)
    nz = kwargs.pop("nz",nz)
    ti,tj,tk,x1,x2,x3,r,h,ph,rho,ug = gd[0:11,:,:].view(); n = 11
    pg = (gam-1)*ug
    lrho=np.log10(rho)
    vu=np.zeros_like(gd[0:4])
    B=np.zeros_like(gd[0:4])
    vu[1:4] = gd[n:n+3]; n+=3
    B[1:4] = gd[n:n+3]; n+=3
    #if electrons are evolved
    if gd.shape[0] == 51 or eCOND == 1 or eHEAT == 1:
      doel = 1
      ktot = gd[n]; n+=1
      kel4 = gd[n]; n+=1
      ugel4 = 3.*kel4*rho**(4./3.)
      Tel4 = kel4*rho**(4./3.)/rho
      if NPR == 18 + DOFLR:
          kel4a = gd[n]; n+=1
          kel4b = gd[n]; n+=1
          kel4c = gd[n]; n+=1
          kel4d = gd[n]; n+=1
          kel4e = gd[n]; n+=1
          ugel4a = 3.*kel4a*rho**(4./3.)
          ugel4b = 3.*kel4b*rho**(4./3.)
          ugel4c = 3.*kel4c*rho**(4./3.)
          ugel4d = 3.*kel4d*rho**(4./3.)
          ugel4e = 3.*kel4e*rho**(4./3.)
          Tel4a = kel4a*rho**(4./3.)/rho
          Tel4b = kel4b*rho**(4./3.)/rho
          Tel4c = kel4c*rho**(4./3.)/rho
          Tel4d = kel4d*rho**(4./3.)/rho
          Tel4e = kel4e*rho**(4./3.)/rho
      kel5 = gd[n]; n+=1
      ugel5 = 3./2.*kel5*rho**(5./3.)
      Tel5 = kel5*rho**(5./3.)/rho
      keldis = gd[n];n+=1
      ugeldis =3.*keldis*rho**(4./3.)
      Teldis =keldis*rho**(4./3.)/rho
      phi = gd[n];n+=1
      if DOFLR: flr = gd[n]; n+=1
      qdot = gd[n]; n+=1
      Ttot = ug*(gam-1.)/rho
      qisosq = gd[n]; n+=1
      qisodotb = gd[n]; n+=1
      pflag = gd[n]; n+=1
      uelvar = gd[n]; n+=1
    else:
      doel = 0
      if DOKTOT == 1:
        ktot = gd[n]; n+=1
      if DOFLR == 1:
        flr = gd[n]; n+=1
      if NEUTRON_STAR == True:
        flrfrac = gd[n]; n+=1
      '''
        print("Using neutron star")
      else:
        print("Not using neutron star")
      '''
    if DONUCLEAR == 1:
      rhonp = gd[n]; n+=1
      rhoalpha = gd[n]; n+=1
      rhofloor = gd[n]; n+=1
      Ye = gd[n]; n+=1
      Tnuc = gd[n]; n+=1
      mneutron_cgs = 1.6749286e-24; c = 2.99792458e10; Mbh_cgs = 3 * 1.99e33; k_cgs = 1.380658e-16
      Tnuc_cgs = Tnuc * mneutron_cgs*c**2 / k_cgs
      etae = gd[n]; n+=1
      G = gd[n]; n+=1
      Q = gd[n]; n+=1
    divb = gd[n]; n+=1
    uu = gd[n:n+4]; n+=4
    ud = gd[n:n+4]; n+=4
    bu = gd[n:n+4]; n+=4
    bd = gd[n:n+4]; n+=4
    bsq = mdot(bu,bd)
    v1m,v1p,v2m,v2p,v3m,v3p=gd[n:n+6]; n+=6
    gdet=gd[n]; n+=1

    rhor = 1+(1-a**2)**0.5
    if 0 and doel:
      phihat = uu-uu;
      phihat[3] = uu[0]/uu[0];
      phihat = phihat/np.sqrt(mdot(mdot(gcov,phihat),phihat))
      csphib = spatialdot(phihat,bu)/np.sqrt(spatialdot(bu,bu)*spatialdot(phihat,phihat))
    if "guu" in globals():
        #lapse
        alpha = (-guu[0,0])**(-0.5)
    if n != gd.shape[0]:
        print("rd: WARNING: nread = %d < ntot = %d: incorrect format?" % (n, gd.shape[0]) )
        return 1
    return 0

def rdump_assign(gd,**kwargs):
    global t,nx,ny,nz,_dx1,_dx2,_dx3,gam,hslope,a,R0,Rin,Rout,ti,tj,tk,x1,x2,x3,r,h,ph,rho,ug,vu,B,pg,cs2,Sden,U,gdetB,divb,uu,ud,bu,bd,v1m,v1p,v2m,v2p,gdet,bsq,gdet,alpha,rhor, qdot, ktot, Ttot, game, qisosq, pflag, qisodotb, kel, uelvar, Tel4, Tel5,Teldis, Tels, kel4, kel5,ugel,ugeldis, ugcon, sel, ugscon, ugel4, ugel5,stot, uelvar, Telvar, Tsel, sel, ugels, games, phi, keldis, phihat,csphib,lrho
    nx = kwargs.pop("nx",nx)
    ny = kwargs.pop("ny",ny)
    nz = kwargs.pop("nz",nz)
    n = 0
    rho = gd[n]; n+=1
    ug = gd[n]; n+=1
    vu=np.zeros_like(gd[0:4])
    B=np.zeros_like(gd[0:4])
    vu[1:4] = gd[n:n+3]; n+=3
    B[1:4] = gd[n:n+3]; n+=3
    # if n != gd.shape[0]:
    #     print("rd: WARNING: nread = %d < ntot = %d: incorrect format?" % (n, gd.shape[0]) )
    #     return 1
    return gd

def fdump_assign(gd,**kwargs):
    global t,nx,ny,nz,_dx1,_dx2,_dx3,gam,hslope,a,R0,Rin,Rout,ti,tj,tk,x1,x2,x3,r,h,ph,rho,ug,vu,B,pg,cs2,Sden,U,gdetB,divb,uu,ud,bu,bd,v1m,v1p,v2m,v2p,gdet,bsq,gdet,alpha,rhor, qdot, ktot, Ttot, game, qisosq, pflag, qisodotb, kel, uelvar, Tel4, Tel5,Teldis, Tels, kel4, kel5,ugel,ugeldis, ugcon, sel, ugscon, ugel4, ugel5,stot, uelvar, Telvar, Tsel, sel, ugels, games, phi, keldis, phihat,csphib,lrho,fail
    nx = kwargs.pop("nx",nx)
    ny = kwargs.pop("ny",ny)
    nz = kwargs.pop("nz",nz)
    fail = gd
    return gd

def re():
    global tvec, rmed, Etot, pp,K_mid,ugmid,ldot,edot,m_dot


    ener = np.loadtxt('ener.out')
    tvec = ener[:,0]
    rmed = ener[:,1]
    pp = ener[:,2]
    Etot = ener[:,3]
    K_mid = ener[:,4]
    ug_mid = ener[:,5]
    m_dot = ener[:,6]
    edot = ener[:,7]
    ldot= ener[:,8]





def mdot(a,b):
    """
    Computes a contraction of two tensors/vectors.  Assumes
    the following structure: tensor[m,n,i,j,k] OR vector[m,i,j,k], 
    where i,j,k are spatial indices and m,n are variable indices. 
    """
    if (a.ndim == 3 and b.ndim == 3) or (a.ndim == 4 and b.ndim == 4):
          c = (a*b).sum(0)
    elif a.ndim == 5 and b.ndim == 4:
          c = np.empty(np.maximum(a[:,0,:,:,:].shape,b.shape),dtype=b.dtype)
          for i in range(a.shape[0]):
                c[i,:,:,:] = (a[i,:,:,:,:]*b).sum(0)
    elif a.ndim == 4 and b.ndim == 5:
          c = np.empty(np.maximum(b[0,:,:,:,:].shape,a.shape),dtype=a.dtype)
          for i in range(b.shape[1]):
                c[i,:,:,:] = (a*b[:,i,:,:,:]).sum(0)
    elif a.ndim == 5 and b.ndim == 5:
          c = np.empty((a.shape[0],b.shape[1],a.shape[2],a.shape[3],max(a.shape[4],b.shape[4])),dtype=a.dtype)
          for i in range(c.shape[0]):
                for j in range(c.shape[1]):
                      c[i,j,:,:,:] = (a[i,:,:,:,:]*b[:,j,:,:,:]).sum(0)
    elif a.ndim == 5 and b.ndim == 6:
          c = np.empty((a.shape[0],b.shape[1],b.shape[2],max(a.shape[2],b.shape[3]),max(a.shape[3],b.shape[4]),max(a.shape[4],b.shape[5])),dtype=a.dtype)
          for mu in range(c.shape[0]):
              for k in range(c.shape[1]):
                  for l in range(c.shape[2]):
                      c[mu,k,l,:,:,:] = (a[mu,:,:,:,:]*b[:,k,l,:,:,:]).sum(0)
    else:
           raise Exception('mdot', 'wrong dimensions')
    return c

def utoT(uel,rho):
    '''
    Converts internal energy to temperture in variable gamma case
    '''
    
    mrat = 1836.152672
    theta = np.zeros(uel.shape)
    for i in range(uel.shape[0]):
        for j in range(uel.shape[1]):
            aco = rho[i,j,0]
            bco = uel[i,j,0]*mrat
            if(bco>=0.):
                theta[i,j,0] = (-6.*aco+5.*bco+np.sqrt(36.*aco*aco+180.*aco*bco+25.*bco*bco))/(30.*aco)
            else:
                theta[i,j,0] =(6.*aco+5.*bco-np.sqrt(36.*aco*aco-180.*aco*bco+25.*bco*bco))/(30.*aco)
    return theta/mrat

        
def psicalc(B1=None):
    """
    Computes the field vector potential
    """
    global B
    if B1 is None: B1 = B[1]
    daphi = -(gdet*B1).mean(-1)*_dx2
    aphi=daphi[:,::-1].cumsum(axis=1)[:,::-1]
    aphi-=0.5*daphi #correction for half-cell shift between face and center in theta
    return(aphi)

def bgradTecalc():
    """
    Computes the Temperature gradient along bu
    """
    dTx = np.diff(Teldis,axis = 0)/_dx1
    z = np.zeros((1,ny,1))
    z[0,:,0] = dTx[nx-2,:,0]
    dTx = np.append(dTx,z,axis = 0)
    dTy = np.diff(Teldis,axis = 1)/_dx2
    z = np.zeros((nx,1,1))
    z[:,0,0] = dTy[:,ny-2,0]
    dTy = np.append(dTy,z, axis = 1)

    return (bu[1]*dTx+bu[2]*dTy)/np.sqrt(bsq)

def gradTemagcalc():
    """
    Computes the Temperature gradient magnitude
    """
    dTx = np.diff(Teldis,axis = 0)/_dx1
    z = np.zeros((1,ny,1))
    z[0,:,0] = dTx[nx-2,:,0]
    dTx = np.append(dTx,z,axis = 0)
    dTy = np.diff(Teldis,axis = 1)/_dx2
    z = np.zeros((nx,1,1))
    z[:,0,0] = dTy[:,ny-2,0]
    dTy = np.append(dTy,z, axis = 1)


    return np.sqrt(dTx**2.+dTy**2.)


def myfloat(f,acc=1):
    """ acc=1 means np.float32, acc=2 means np.float64 """
    if acc==1:
        return( np.float32(f) )
    else:
        return( np.float64(f) )

def get_fracphi():
    fracphi = dxdxp[3,3,0,0,0]*_dx3*nz/(2*np.pi)
    return( fracphi )

def plco(myvar,**kwargs):
    global r,h,ph
    ax = kwargs.setdefault("ax",None)
    if ax is not None:
        ax.cla()
    else:
        plt.clf()
    return plc(myvar,**kwargs)

def plc(myvar,**kwargs): #plc
    global r,h,ph
    #xcoord = kwargs.pop('x1', None)
    #ycoord = kwargs.pop('x2', None)
    if(np.min(myvar)==np.max(myvar)):
        print("The quantity you are trying to plot is a constant = %g." % np.min(myvar))
        return
    cb = kwargs.pop('cb', False)
    cbpad = kwargs.pop('cbpad', 0.03)
    nc = kwargs.pop('nc', 15)
    k = kwargs.pop('k',0)
    mirrorx = kwargs.pop('mirrorx',0)
    mirrory = kwargs.pop('mirrory',0)
    symmx = kwargs.pop('symmx',0)
    #cmap = kwargs.pop('cmap',cm.jet)
    isfilled = kwargs.pop('isfilled',False)
    xy = kwargs.pop('xy',0)
    xcoord = kwargs.pop("xcoord",None)
    ycoord = kwargs.pop("ycoord",None)
    lin = kwargs.pop('lin',0)
    xmax = kwargs.pop('xmax',10)
    ymax = kwargs.pop('ymax',5)
    cbxlabel = kwargs.pop('cbxla',None)
    cbylabel = kwargs.pop('cbyla',None)
    fntsize = kwargs.pop("fntsize",20)
    cbgoodticks = kwargs.pop("cbgoodticks",1)
    xlabel = kwargs.pop("xla",None)
    ylabel = kwargs.pop("yla",None)
    dobh = kwargs.pop("dobh",1)
    pretty = kwargs.pop("pretty",0)
    ax = kwargs.pop("ax",None)
    cbticks = kwargs.pop("cbticks",None)
    domathify = kwargs.pop("domathify",0)
    delta_phi = kwargs.pop("delta_phi",0.)
    if np.abs(xy)==1:
        if xcoord is None: xcoord = r * np.sin(h)
        if ycoord is None: ycoord = r * np.cos(h)
        if mirrory: ycoord *= -1
        if mirrorx: xcoord *= -1
    if xcoord is not None and ycoord is not None:
        xcoord = xcoord[:,:,None] if xcoord.ndim == 2 else xcoord[:,:,k:k+1]
        ycoord = ycoord[:,:,None] if ycoord.ndim == 2 else ycoord[:,:,k:k+1]
    if np.abs(xy)==1 and symmx:
        if myvar.ndim == 2:
            myvar = myvar[:,:,None] if myvar.ndim == 2 else myvar[:,:,k:k+1]
            myvar=np.concatenate((myvar[:,::-1],myvar),axis=1)
            xcoord=np.concatenate((-xcoord[:,::-1],xcoord),axis=1)
            ycoord=np.concatenate((ycoord[:,::-1],ycoord),axis=1)
        else:
            if myvar.shape[-1] > 1: 
                symmk = (k+nz/2)%nz 
            else: 
                symmk = k
            myvar=np.concatenate((myvar[:,ny-1:ny,k:k+1],myvar[:,::-1,symmk:symmk+1],myvar[:,:,k:k+1]),axis=1)
            xcoord=np.concatenate((xcoord[:,ny-1:ny,k:k+1],-xcoord[:,::-1],xcoord),axis=1)
            ycoord=np.concatenate((ycoord[:,ny-1:ny,k:k+1],ycoord[:,::-1],ycoord),axis=1)
    elif np.abs(xy) == 2 and symmx:
        #if fracphi == 0.5 done in a robust way
        if get_fracphi() < 0.75:
            r1 = np.concatenate((r,r,r[...,0:1]),axis=2)
            ph1 = np.concatenate((ph,ph+np.pi,ph[...,0:1]+2*np.pi),axis=2)
            myvar = np.concatenate((myvar,myvar,myvar[...,0:1]),axis=2)
        else:
            r1 = np.concatenate((r,r[...,0:1]),axis=2)
            ph1 = np.concatenate((ph,ph[...,0:1]+2*np.pi),axis=2)
            myvar = np.concatenate((myvar,myvar[...,0:1]),axis=2)
        if delta_phi != 0.:
            ph1 += delta_phi
        xcoord=(r1*cos(ph1))[:,ny/2,:,None]
        ycoord=(r1*sin(ph1))[:,ny/2,:,None]
        myvar = myvar[:,ny/2,:,None]
    else:
        myvar = myvar[:,:,None] if myvar.ndim == 2 else myvar[:,:,k:k+1]
    if lin:
        xcoord = r
        ycoord = h
    if ax is None:
        ax = plt.gca()
    if  xcoord is None or ycoord is None:
        if isfilled:
            res = ax.contourf(myvar[:,:,0].transpose(),nc,**kwargs)
        else:
            res = ax.contour(myvar[:,:,0].transpose(),nc,**kwargs)
    else:
        if isfilled:
            res = ax.contourf(xcoord[:,:,0],ycoord[:,:,0],myvar[:,:,0],nc,**kwargs)
        else:
            res = ax.contour(xcoord[:,:,0],ycoord[:,:,0],myvar[:,:,0],nc,**kwargs)
    if xy>0 and not symmx:
        ax.set_xlim(0,xmax)
        ax.set_ylim(-ymax,ymax)
    if xy> 0 and symmx:
        ax.set_xlim(-xmax,xmax)
        ax.set_ylim(-ymax,ymax)
    if xlabel is not None:
        ax.set_xlabel(xlabel,fontsize=fntsize)
    if ylabel is not None:
        ax.set_ylabel(ylabel,fontsize=fntsize)
    if pretty:
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontsize(fntsize)
            if domathify: mathify_axes_ticks(ax,fontsize=fntsize)
    if cb: #use color bar
        cb = plt.colorbar(res,ax=ax,pad=cbpad)
        if pretty and cbgoodticks and cbticks is None:
            vmin = cb.vmin
            vmax = cb.vmax
            #this returns incorrect ticks! so ignore it
            #ticks = cb.ax.get_yticks()
            #nticks = len(ticks)
            #if not too many ticks, then pretty them up
            rvmin = np.round(vmin)
            rvmax = np.round(vmax)
            if rvmin == vmin and rvmax == vmax and vmax-vmin <= 10:
                ticks = np.arange(rvmin,rvmax+1)
                cb.set_ticks(ticks)
                mathify_axes_ticks(cb.ax,fontsize=fntsize,yticks=ticks)
        if cbticks is not None:
            cb.set_ticks(cbticks)
            mathify_axes_ticks(cb.ax,fontsize=fntsize,yticks=cbticks)
        if cbxlabel is not None:
            cb.ax.set_xlabel(cbxlabel,fontsize=fntsize)
        if cbxlabel is not None:
            cb.ax.set_xlabel(cbxlabel,fontsize=fntsize)
        if cbylabel is not None:
            cb.ax.set_ylabel(cbylabel,fontsize=fntsize)
        if pretty:
            for label in cb.ax.get_yticklabels():
                label.set_fontsize(fntsize)
    if xy and dobh and "rhor" in globals(): 
        el = Ellipse((0,0), 2*rhor, 2*rhor, facecolor='k', alpha=1)
        art=ax.add_artist(el)
        art.set_zorder(20)
    if cb:
        return res, cb
    else:
        return res
    
    

def faraday():
    global fdd, fuu, omegaf1, omegaf2, omegaf1b, omegaf2b, rhoc
    if 'fdd' in globals():
        del fdd
    if 'fuu' in globals():
        del fuu
    if 'omegaf1' in globals():
        del omegaf1
    if 'omemaf2' in globals():
        del omegaf2
    # these are native values according to HARM
    fdd = np.zeros((4,4,nx,ny,nz),dtype=rho.dtype)
    #fdd[0,0]=0*gdet
    #fdd[1,1]=0*gdet
    #fdd[2,2]=0*gdet
    #fdd[3,3]=0*gdet
    fdd[0,1]=gdet*(uu[2]*bu[3]-uu[3]*bu[2]) # f_tr
    fdd[1,0]=-fdd[0,1]
    fdd[0,2]=gdet*(uu[3]*bu[1]-uu[1]*bu[3]) # f_th
    fdd[2,0]=-fdd[0,2]
    fdd[0,3]=gdet*(uu[1]*bu[2]-uu[2]*bu[1]) # f_tp
    fdd[3,0]=-fdd[0,3]
    fdd[1,3]=gdet*(uu[2]*bu[0]-uu[0]*bu[2]) # f_rp = gdet*B2
    fdd[3,1]=-fdd[1,3]
    fdd[2,3]=gdet*(uu[0]*bu[1]-uu[1]*bu[0]) # f_hp = gdet*B1
    fdd[3,2]=-fdd[2,3]
    fdd[1,2]=gdet*(uu[0]*bu[3]-uu[3]*bu[0]) # f_rh = gdet*B3
    fdd[2,1]=-fdd[1,2]
    #
    fuu = np.zeros((4,4,nx,ny,nz),dtype=rho.dtype)
    #fuu[0,0]=0*gdet
    #fuu[1,1]=0*gdet
    #fuu[2,2]=0*gdet
    #fuu[3,3]=0*gdet
    fuu[0,1]=-1/gdet*(ud[2]*bd[3]-ud[3]*bd[2]) # f^tr
    fuu[1,0]=-fuu[0,1]
    fuu[0,2]=-1/gdet*(ud[3]*bd[1]-ud[1]*bd[3]) # f^th
    fuu[2,0]=-fuu[0,2]
    fuu[0,3]=-1/gdet*(ud[1]*bd[2]-ud[2]*bd[1]) # f^tp
    fuu[3,0]=-fuu[0,3]
    fuu[1,3]=-1/gdet*(ud[2]*bd[0]-ud[0]*bd[2]) # f^rp
    fuu[3,1]=-fuu[1,3]
    fuu[2,3]=-1/gdet*(ud[0]*bd[1]-ud[1]*bd[0]) # f^hp
    fuu[3,2]=-fuu[2,3]
    fuu[1,2]=-1/gdet*(ud[0]*bd[3]-ud[3]*bd[0]) # f^rh
    fuu[2,1]=-fuu[1,2]
    #
    # these 2 are equal in degen electrodynamics when d/dt=d/dphi->0
    omegaf1=fdd[0,1]/fdd[1,3] # = ftr/frp
    omegaf2=fdd[0,2]/fdd[2,3] # = fth/fhp
    #
    # from jon branch, 04/10/2012
    #
    if 0:
        B1hat=B[1]*np.sqrt(gv3[1,1])
        B2hat=B[2]*np.sqrt(gv3[2,2])
        B3nonhat=B[3]
        v1hat=uu[1]*np.sqrt(gv3[1,1])/uu[0]
        v2hat=uu[2]*np.sqrt(gv3[2,2])/uu[0]
        v3nonhat=uu[3]/uu[0]
        #
        aB1hat=np.fabs(B1hat)
        aB2hat=np.fabs(B2hat)
        av1hat=np.fabs(v1hat)
        av2hat=np.fabs(v2hat)
        #
        vpol=np.sqrt(av1hat**2 + av2hat**2)
        Bpol=np.sqrt(aB1hat**2 + aB2hat**2)
        #
        #omegaf1b=(omegaf1*aB1hat+omegaf2*aB2hat)/(aB1hat+aB2hat)
        #E1hat=fdd[0,1]*np.sqrt(gn3[1,1])
        #E2hat=fdd[0,2]*np.sqrt(gn3[2,2])
        #Epabs=np.sqrt(E1hat**2+E2hat**2)
        #Bpabs=np.sqrt(aB1hat**2+aB2hat**2)+1E-15
        #omegaf2b=Epabs/Bpabs
        #
        # assume field swept back so omegaf is always larger than vphi (only true for outflow, so put in sign switch for inflow as relevant for disk near BH or even jet near BH)
        # GODMARK: These assume rotation about z-axis
        omegaf2b=np.fabs(v3nonhat) + np.sign(uu[1])*(vpol/Bpol)*np.fabs(B3nonhat)
        #
        omegaf1b=v3nonhat - B3nonhat*(v1hat*B1hat+v2hat*B2hat)/(B1hat**2+B2hat**2)
    #
    # charge
    #
    if 0:
        rhoc = np.zeros_like(rho)
        if nx>=2:
            rhoc[1:-1] += ((gdet*fuu[0,1])[2:]-(gdet*fuu[0,1])[:-2])/(2*_dx1)
        if ny>2:
            rhoc[:,1:-1] += ((gdet*fuu[0,2])[:,2:]-(gdet*fuu[0,2])[:,:-2])/(2*_dx2)
        if ny>=2 and nz > 1: #not sure if properly works for 2D XXX
            rhoc[:,0,:nz/2] += ((gdet*fuu[0,2])[:,1,:nz/2]+(gdet*fuu[0,2])[:,0,nz/2:])/(2*_dx2)
            rhoc[:,0,nz/2:] += ((gdet*fuu[0,2])[:,1,nz/2:]+(gdet*fuu[0,2])[:,0,:nz/2])/(2*_dx2)
        if nz>2:
            rhoc[:,:,1:-1] += ((gdet*fuu[0,3])[:,:,2:]-(gdet*fuu[0,3])[:,:,:-2])/(2*_dx3)
        if nz>=2:
            rhoc[:,:,0] += ((gdet*fuu[0,3])[:,:,1]-(gdet*fuu[0,3])[:,:,-1])/(2*_dx3)
            rhoc[:,:,-1] += ((gdet*fuu[0,3])[:,:,0]-(gdet*fuu[0,3])[:,:,-2])/(2*_dx3)
        rhoc /= gdet

def Tcalcud():
    global Tud, TudEM, TudMA
    global mu, sigma
    global enth
    global unb, isunbound
    pg = (gam-1)*ug
    w=rho+ug+pg
    eta=w+bsq
    if 'Tud' in globals():
        del Tud
    if 'TudMA' in globals():
        del TudMA
    if 'TudEM' in globals():
        del TudEM
    if 'mu' in globals():
        del mu
    if 'sigma' in globals():
        del sigma
    if 'unb' in globals():
        del unb
    if 'isunbound' in globals():
        del isunbound
    Tud = np.zeros((4,4,nx,ny,nz),dtype=np.float32,order='F')
    TudMA = np.zeros((4,4,nx,ny,nz),dtype=np.float32,order='F')
    TudEM = np.zeros((4,4,nx,ny,nz),dtype=np.float32,order='F')
    for kapa in np.arange(4):
        for nu in np.arange(4):
            if(kapa==nu): delta = 1
            else: delta = 0
            TudEM[kapa,nu] = bsq*uu[kapa]*ud[nu] + 0.5*bsq*delta - bu[kapa]*bd[nu]
            TudMA[kapa,nu] = w*uu[kapa]*ud[nu]+pg*delta
            #Tud[kapa,nu] = eta*uu[kapa]*ud[nu]+(pg+0.5*bsq)*delta-bu[kapa]*bd[nu]
            Tud[kapa,nu] = TudEM[kapa,nu] + TudMA[kapa,nu]
    mu = -Tud[1,0]/(rho*uu[1])
    sigma = TudEM[1,0]/TudMA[1,0]
    enth=1+ug*gam/rho
    unb=enth*ud[0]
    isunbound=(-unb>1.0)


def aux():
    faraday()
    Tcalcud()

if __name__ == "__main__":
    if False:
        #1D plot example
        plt.clf()
        rg("gdump")
        rd("dump000")
        plt.plot(r[:,ny/2,0],rho[:,ny/2,0])
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("r")
        plt.ylabel("rho")
    if False:
        #2D plot example
        plt.clf()
        rg("gdump")
        rd("dump000")
        #R-z plot of the logarithm of density distribution
        plc(r,np.log10(rho),cb=True,xy=1,xmax=100,ymax=50)

def ts(startn=0,endn=-1,whichi=0,whichj=0):
    flist1 = np.sort(glob.glob(os.path.join(dumps_path, "dump[0-9][0-9][0-9]") ) )
    flist2 = np.sort(glob.glob(os.path.join(dumps_path, "dump[0-9][0-9][0-9][0-9]") ) )
    flist1.sort()
    flist2.sort()
    flist = flist1.tolist() + flist2.tolist()
    rg("gdump")
    dic = {}
    dic["t"] = []
    dic["ktot"] = []
    dic["Ttot"] = []
    dic["Source"] = []
    dic["rho"] = []
    dic["qisosq"] = []
    dic["qisodotb"] = []
    dic["vx"] = []
    dic["whichi"] = whichi
    dic["whichj"] = whichi
    dic["vy"] = []
    dic["ug"] = []
    dic["ugscon"] = []
    dic["sel"] = []
    dic["kel4"] = []
    dic["kel5"] = []
    dic["keldis"] = []
    dic["uelvar"] = []
    dic["bsq"] = []

    if NEUTRON_STAR == True:
        dic["flrfrac"] = []

    for fname in flist:
        #find the index of the file
        findex = np.int(fname.split("p")[-1])
        if findex < startn:
            continue
        if endn>=0 and findex >= endn:
            break
        dumpname = "dump%03d" % findex
        print( "Reading from %s..." % dumpname )
        rd(dumpname)
        f = lambda x: x[:,:,0] #x[whichi,whichj,0]
        dic["ktot"].append( f(ktot) )
        dic["Ttot"].append( f(Ttot) )
        dic["Source"].append(f(qdot))
        dic["rho"].append(f(rho))
        dic["qisosq"].append(f(qisosq))
        dic["qisodotb"].append(f(qisodotb))
        dic["t"].append(t)
        dic["vx"].append(f(vu[1,:,:,:]))
        dic["vy"].append(f(vu[2,:,:,:]))
        dic["ug"].append(f(ug))
        dic["sel"].append(f(sel))
        dic["ugscon"].append(f(ugscon))
        dic["kel4"].append(f(kel4))
        dic["kel5"].append(f(kel5))
        dic["keldis"].append(f(keldis))
        dic["uelvar"].append(f(uelvar))
        dic["bsq"].append(mdot(f(bu),f(bd)))

        if NEUTRON_STAR == True:
            dic["flrfrac"].append( f(flrfrac) )

    return dic

def Tavg(startn=0,endn=-1,fname = "avg.npz"):
    rg('gdump')
    rd('dump000')
    divfac = (endn-startn+1)+0.0
    Tel4avg = Tel4-Tel4
    Teldisavg = Ttot -Ttot
    Tdisratavg = Ttot-Ttot
    qratavg = Ttot-Ttot
    rhoavg = rho-rho
    Tcondrat = rho-rho
    psiavg = rho-rho
    bsqavg = bsq-bsq
    Bavg = B-B
    for n in range(startn,endn+1):
        rd("dump%03d"%n)
        Tel4avg = Tel4avg + Tel4
        Teldisavg = Teldisavg + Teldis
        rhoavg = rho + rhoavg
        psiavg = psiavg + psicalc()
        qratavg = qratavg + qint/qdot
        Tcondrat = Tcondrat + Tel4/Teldis
        Tdisratavg = Tdisratavg + Teldis/Ttot
        bsqavg = bsqavg + bsq
        Bavg = Bavg + B
    
    
    Tel4avg = Tel4avg/divfac
    Teldisavg = Teldisavg/divfac
    rhoavg = rhoavg/divfac
    qratavg = qratavg/divfac
    psiavg = psiavg/divfac
    Tdisratavg = Tdisratavg/divfac
    Tcondrat = Tcondrat/divfac
    bsqavg = bsqavg/divfac
    Bavg = Bavg/divfac

    avg ={"rho":rhoavg,"psi":psiavg,"Tdis":Teldisavg,"Tel4":Tel4avg,"qrat":qratavg,"Tdisrat":Tdisratavg,"Tcrat":Tcondrat,"bsq":bsqavg,"B":Bavg }

    np.savez(fname,**avg)
    return avg

def bhole():

    ax = plt.gca()
    el = Ellipse((0,0), 2*rhor, 2*rhor, facecolor='k', alpha=1)
    art=ax.add_artist(el)
    art.set_zorder(20)
    plt.draw()


#
# 1D averaging
#
def mkts_wrapper(prefix = "qty"):
    if len(sys.argv[2:])>=2 and sys.argv[2].isdigit() and sys.argv[3].isdigit():
        whichi = int(sys.argv[2])
        whichn = int(sys.argv[3])
        endn = -1
        if len(sys.argv[2:])==3:
            if sys.argv[4].isdigit():
                endn = int(sys.argv[4])
        if whichi < whichn:
            postprocess1d(endn = endn, whichi = whichi, whichn = whichn, prefix = prefix)
        else:
            mrgnew(whichn, prefix = prefix)
    else:
        print("Syntax error")

def mrgnew(n=None,fin1=None,fin2=None,fout="qty.npz",**kwargs):
    #setup empty variables first (not necessary, here to remember what they are)
    v = {}
    if n is not None:
        #merge a numbered sequence of files
        v = mrgnew_n(n,v,**kwargs)
    elif fin1 is not None and fin2 is not None:
        #merge two named files
        v = mrgnew_f(fin1,v,**kwargs)
        v = mrgnew_f(fin2,v,**kwargs)
    print( "Saving into " + fout + " ..." )
    sys.stdout.flush()
    np.savez(fout, **v)
    print( "Done!" )

def mrgnew_n(n, v={},**kwargs):
    prefix = kwargs.pop("prefix","qty")
    if v is None:
        v = {}
    for i in np.arange(n):
        #load each file
        ft = "%s_%02d_%02d.npz" % (prefix, i, n)
        print( "Loading " + ft + " ..." )
        sys.stdout.flush()
        vt=np.load( ft )
        for key in vt.keys():
            if key not in ["ivals", "rvals"]:
                if key not in v:
                    v[key] = []
                v[key] += list(vt[key])
            else:
                #no time dependence for these, so same for all times and no need to append
                v[key] = list(vt[key])
        vt.close()    
    #now sort the resulting list
    imap = np.argsort(v["ind"])
    #sort keys that have time series
    for key in vt.keys():
        if key not in ["ivals", "rvals"]:
            if len(v[key]) > 0:
                v[key] = np.array(v[key])[imap,...]
    #hack to reintroduce the time if missing
    if "t" not in v or len(v["t"]) == 0:
        print( "Times are missing, assuming dumping period of 5 and setting t = 5*ind" )
        v["t"] = np.array(v["ind"])*5.
    return(v)

def mrgnew_f(ft, v={}, **kwargs):
    prefix = kwargs.pop("prefix","qty")
    #assumes that the dics in ft and v are sorted in time
    #load file
    print( "Loading " + ft + " ..." )
    sys.stdout.flush()
    vt=np.load( ft )
    #find last element in old array that does not exist in the new array
    if "ind" in v:
        ind = np.array(v["ind"])
        lasti = (ind<vt["ind"][0]).sum()
        FMbackup = np.copy(v["FM"])
    else:
        lasti = None
        ind = None
    for key in vt.keys():
        if key not in ["ivals", "rvals"]:
            if key not in v:
                if ind is not None:
                    #new field, set it to zeros for previoius time
                    v[key] = np.zeros_like(FMbackup)
                else:
                    #empty dictionary
                    v[key] = []
            #convert to list first
            v[key] = list(v[key])
            #has time dependence, so first discard repeated entries
            if lasti is not None: del v[key][lasti:]
            #and then append new entries
            v[key] += list(vt[key])
        else:
            v[key] = list(vt[key])
    vt.close()    
    #now sort the resulting list
    imap = np.argsort(v["ind"])
    #sort keys that have time series
    for key in vt.keys():
        if key not in ["ivals", "rvals"]:
            if len(v[key]) > 0:
                v[key] = np.array(v[key])[imap,...]
    #hack to reintroduce the time if missing
    if "t" not in v or len(v["t"]) == 0:
        print( "Times are missing, assuming dumping period of 5 and setting t = 5*ind" )
        v["t"] = np.array(v["ind"])*5.
    return(v)

def testfail(fldname = "dump000"):
    try: 
        rd(fldname)
    except IOError as e:
        print ("I/O error({0}): {1}".format(e.errno, e.strerror))


def get_sorted_file_list(prefix="dump"):
    flist0 = np.sort(glob.glob( os.path.join(dumps_path, "%s[0-9][0-9][0-9]"%prefix) ) )
    flist1 = np.sort(glob.glob( os.path.join(dumps_path, "%s[0-9][0-9][0-9][0-9]"%prefix) ) )
    flist2 = np.sort(glob.glob( os.path.join(dumps_path, "%s[0-9][0-9][0-9][0-9][0-9]"%prefix) ) )
    flist0.sort()
    flist1.sort()
    flist2.sort()
    flist = np.concatenate((flist0,flist1,flist2))
    return flist
    
def postprocess1d(startn=0,endn=-1,whichi=0,whichn=1,**kwargs):
    prefix = kwargs.pop("prefix","qty")
    rg("gdump")
    flist = get_sorted_file_list()
    v = {}
    #initialize some variables that are definitely going to be there
    v["t"] = []
    v["ind"] = []
    #the rest will be whatever compvals1d() returns
    fname = "%s_%02d_%02d.npz" % (prefix, whichi, whichn)
    if os.path.isfile( fname ):
        print( "File %s exists, loading from file..." % fname )
        sys.stdout.flush()
        vt=np.load( fname )
        for key in vt:
            v[key] = list(vt[key])
        vt.close()
    for fldname in flist:
        #find the index of the file
        ind = np.int(fldname.split("p")[-1])
        if ind < startn:
            continue
        if endn>=0 and ind >= endn:
            break
        if ind % whichn != whichi:
            #do every whichn'th snapshot starting with whichi'th snapshot
            continue
        if ind in v["ind"]:
            print( "File %04d already processed, skipping..." % ind )
            continue
        print( "Reading " + fldname + " ..." )
        sys.stdout.flush()
        try: 
            rd("../"+fldname)
        except IOError as e:
            print( "While reading %s:" % fldname )
            print( "I/O error({0}): {1}".format(e.errno, e.strerror) )
            print( "Skipping" )
            continue
        #cvellite()  #now unnecessary
        sys.stdout.flush()
        valdic = compvals1d()
        for key in valdic:
            if key == "ivals" or key == "rvals":
                #store only one copy of radial evaluation points
                v[key] = valdic[key]
            else:
                if key not in v:
                    v[key] = []
                v[key].append(valdic[key])
        v["ind"].append(ind)
        v["t"].append(t)
    np.savez("%s_%02d_%02d.npz" % (prefix, whichi, whichn), **v)

def compvals1d(di=5,hor=0.1,docool=1):
    dic = {}
    #rvals = np.array([rhor,5.,10.,20.,50.,100.])
    #ivals = np.int32(iofr(rvals))
    ihor = iofr(rhor)
    i0 = ihor%di
    #every di'th radial cell including the event horizon cell, so overall nx/di ~ 50 cells
    ivals = np.int32(ti[i0::di,ny/2,0]+0.5)
    rvals = r[i0::di,ny/2,0]
    #
    delta = lambda kapa,nu: (kapa==nu)
    fTudEM = lambda kapa,nu: bsq*uu[kapa]*ud[nu] + 0.5*bsq*delta(kapa,nu) - bu[kapa]*bd[nu]
    fTudMA = lambda kapa,nu: (rho+gam*ug)*uu[kapa]*ud[nu]+(gam-1)*ug*delta(kapa,nu)
    fTud = lambda kapa,nu: fTudEM(kapa,nu) + fTudMA(kapa,nu)
    fRud = lambda kapa,nu: 4./3.*Erf*uradu[kapa]*uradd[nu]+1./3.*Erf*delta(kapa,nu)
    #
    if "gdetF" in globals() and gdetF is not None:
        FM = -(gdetF[1,0]).sum(-1).sum(-1)*_dx2*_dx3
        dic["FM"] = FM[ivals]
        #this includes Mdot
        FE = (gdetF[1,1]).sum(-1).sum(-1)*_dx2*_dx3
        #subtract Mdot
        FE += FM
        dic["FE"] = FE[ivals]
    else:
        FM = -(gdet*rho*uu[1]).sum(-1).sum(-1)*_dx2*_dx3
        dic["FM"] = FM[ivals]
        #this includes Mdot
        FE = (gdet*fTud(1,0)).sum(-1).sum(-1)*_dx2*_dx3
        #no need to subtract Mdot here
        ##### no need for this ##### FE += FM
        dic["FE"] = FE[ivals]
    if docool:
        Ttarg = Tnt(a,r[:,ny/2:ny/2+1,0:1],hor=0.1)
        iscold = ((gam-1)*ug/rho) < np.exp(1)*Ttarg
        ishot = ~iscold
        dic["FMcold"] = -(rho*uu[1]*iscold)[ivals].sum(-1).sum(-1)*_dx2*_dx3
        dic["FMhot"] = -(rho*uu[1]*ishot)[ivals].sum(-1).sum(-1)*_dx2*_dx3
        dic["FMcoldin"] = -(rho*uu[1]*iscold*(rho*uu[1]<0))[ivals].sum(-1).sum(-1)*_dx2*_dx3
        dic["FMcoldout"] = -(rho*uu[1]*iscold*(rho*uu[1]>=0))[ivals].sum(-1).sum(-1)*_dx2*_dx3
        dic["FMhotin"] = -(rho*uu[1]*ishot*(rho*uu[1]<0))[ivals].sum(-1).sum(-1)*_dx2*_dx3
        dic["FMhotout"] = -(rho*uu[1]*ishot*(rho*uu[1]>=0))[ivals].sum(-1).sum(-1)*_dx2*_dx3
        dic["Sigcoldin"] = (gdet*rho*uu[0]*iscold*(rho*uu[1]<0))[ivals].mean(-1).sum(-1)*_dx2
        dic["Sighotin"] = (gdet*rho*uu[0]*ishot*(rho*uu[1]<0))[ivals].mean(-1).sum(-1)*_dx2
        dic["Sigcoldout"] = (gdet*rho*uu[0]*iscold*(rho*uu[1]>=0))[ivals].mean(-1).sum(-1)*_dx2
        dic["Sighotout"] = (gdet*rho*uu[0]*ishot*(rho*uu[1]>=0))[ivals].mean(-1).sum(-1)*_dx2
    FEM = (gdet*fTudEM(1,0)).sum(-1).sum(-1)*_dx2*_dx3
    dic["FEM"] = FEM[ivals]
    #radiation
    if "uradu" in globals():
        FR = (gdet*fRud(1,0)).sum(-1).sum(-1)*_dx2*_dx3
        FRj = (gdet*fRud(1,0)*(bsq/rho>10)).sum(-1).sum(-1)*_dx2*_dx3
        dic["FR"] = FR[ivals]
        dic["FRj"] = FRj[ivals]
        tau_dic = compute_taurad()
        tau1 = tau_dic["tau1"]
        tau2 = tau_dic["tau2"]
        dic["FRtau1=1"] = (gdet*fRud(1,0)*(tau1<1))[ivals].sum(-1).sum(-1)*_dx2*_dx3
        dic["FRtau1=2o3"] = (gdet*fRud(1,0)*(tau1<2./3.))[ivals].sum(-1).sum(-1)*_dx2*_dx3
        dic["FRtau2=1"] = (gdet*fRud(1,0)*(tau2<1))[ivals].sum(-1).sum(-1)*_dx2*_dx3
        dic["FRtau2=2o3"] = (gdet*fRud(1,0)*(tau2<2./3.))[ivals].sum(-1).sum(-1)*_dx2*_dx3
    #total absolute magnetic flux
    PhiBH = 0.5*np.abs(gdet*B[1]).sum(-1).sum(-1)*_dx2*_dx3
    PhiBHnorth = 0.5*(gdet*B[1])[:,:ny/2].sum(-1).sum(-1)*_dx2*_dx3
    PhiBHsouth = 0.5*(gdet*B[1])[:,ny/2:].sum(-1).sum(-1)*_dx2*_dx3
    dic["PhiBH"] = PhiBH[ivals]
    dic["PhiBHnorth"] = PhiBHnorth[ivals]
    dic["PhiBHsouth"] = PhiBHsouth[ivals]
    dic["rvals"] = rvals
    dic["ivals"] = ivals
    diskcondition = bsq/rho<10
    dic["hor"] = horcalc1d(which=diskcondition)[ivals]
    #jet power and mass outflow in the jets
    #SASMARK: UNBOUND CRITERION INCORRECT WITH RADIATION
    #Shouldn't compare radiation and gas energies in the same frame??! 
    #urad and ug are in different frames!
    isunb=(-(1+ug*gam/rho)*ud[0]>1.0)
    mu = -fTud(1,0)/(rho*uu[1])
    #EM-defined jet (minmu > ...)
    dic["FEMMAmu>2"] = (gdet*fTud(1,0)*(mu>2)*(isunb))[ivals].sum(-1).sum(-1)*_dx2*_dx3
    dic["FEMMAmu>1"] = (gdet*fTud(1,0)*(mu>1)*(isunb))[ivals].sum(-1).sum(-1)*_dx2*_dx3
    dic["FEMmu>2"] = (gdet*fTudEM(1,0)*(mu>2)*(isunb))[ivals].sum(-1).sum(-1)*_dx2*_dx3
    dic["FEMmu>1"] = (gdet*fTudEM(1,0)*(mu>1)*(isunb))[ivals].sum(-1).sum(-1)*_dx2*_dx3
    dic["FMAmu>2"] = (gdet*fTudMA(1,0)*(mu>2)*(isunb))[ivals].sum(-1).sum(-1)*_dx2*_dx3
    dic["FMAmu>1"] = (gdet*fTudMA(1,0)*(mu>1)*(isunb))[ivals].sum(-1).sum(-1)*_dx2*_dx3
    dic["FRMmu>2"] = (gdet*rho*uu[1]*(mu>2)*(isunb))[ivals].sum(-1).sum(-1)*_dx2*_dx3
    dic["FRMmu>1"] = (gdet*rho*uu[1]*(mu>1)*(isunb))[ivals].sum(-1).sum(-1)*_dx2*_dx3
    dic["Phimu>2"] = (gdet*B[1]*(mu>2)*(isunb))[ivals].sum(-1).sum(-1)*_dx2*_dx3
    dic["Phimu>1"] = (gdet*B[1]*(mu>1)*(isunb))[ivals].sum(-1).sum(-1)*_dx2*_dx3
    if "uradu" in globals():
        isunb=(-(1+(ug*gam+urad)/rho)*ud[0]>1.0)
        isbnd=~isunb
        murad = -(fTud(1,0)+fRud(1,0))/(rho*uu[1])
        # radiation fluxes
        dic["FRunb"]= (gdet*fRud(1,0)*(isunb))[ivals].sum(-1).sum(-1)*_dx2*_dx3
        dic["FRmurad>1"]= (gdet*fRud(1,0)*(murad>1)*(isunb))[ivals].sum(-1).sum(-1)*_dx2*_dx3
        dic["FRmurad>2"]= (gdet*fRud(1,0)*(murad>2)*(isunb))[ivals].sum(-1).sum(-1)*_dx2*_dx3
        dic["FRmu>1"]= (gdet*fRud(1,0)*(mu>1)*(isunb))[ivals].sum(-1).sum(-1)*_dx2*_dx3
        dic["FRmu>2"]= (gdet*fRud(1,0)*(mu>2)*(isunb))[ivals].sum(-1).sum(-1)*_dx2*_dx3
        # radiation-defined jet (minmurad > ...)
        dic["FEMMAmurad>2"] = (gdet*fTud(1,0)*(murad>2)*(isunb))[ivals].sum(-1).sum(-1)*_dx2*_dx3
        dic["FEMMAmurad>1"] = (gdet*fTud(1,0)*(murad>1)*(isunb))[ivals].sum(-1).sum(-1)*_dx2*_dx3
        dic["FEMmurad>2"] = (gdet*fTudEM(1,0)*(murad>2)*(isunb))[ivals].sum(-1).sum(-1)*_dx2*_dx3
        dic["FEMmurad>1"] = (gdet*fTudEM(1,0)*(murad>1)*(isunb))[ivals].sum(-1).sum(-1)*_dx2*_dx3
        dic["FMAmurad>2"] = (gdet*fTudMA(1,0)*(murad>2)*(isunb))[ivals].sum(-1).sum(-1)*_dx2*_dx3
        dic["FMAmurad>1"] = (gdet*fTudMA(1,0)*(murad>1)*(isunb))[ivals].sum(-1).sum(-1)*_dx2*_dx3
        dic["FRMmurad>2"] = (gdet*rho*uu[1]*(murad>2)*(isunb))[ivals].sum(-1).sum(-1)*_dx2*_dx3
        dic["FRMmurad>1"] = (gdet*rho*uu[1]*(murad>1)*(isunb))[ivals].sum(-1).sum(-1)*_dx2*_dx3
        dic["Phimurad>2"] = (gdet*B[1]*(murad>2)*(isunb))[ivals].sum(-1).sum(-1)*_dx2*_dx3
        dic["Phimurad>1"] = (gdet*B[1]*(murad>1)*(isunb))[ivals].sum(-1).sum(-1)*_dx2*_dx3
        #average energy transport velocity
        dic["vE"] =( (gdet*uu*(ug+urad)*isbnd).mean(-1).mean(-1)[:,ivals]
                   / (gdet*   (ug+urad)*isbnd).mean(-1).mean(-1)[ivals] )
    return( dic )

    #mu = -fTud(1,0)/(rho*uu[1])

    #left here for reference:
    # edtotvsr[:-1] = -0.5*(gdetF11[:-1]+gdetF11[1:])
    # mdtotvsr[:-1] = -0.5*(gdetF10[:-1]+gdetF10[1:])
    # edtotvsr -= mdtotvsr

def gaussf( x, x0, sigma ):
    """ Returns normalized Gaussian centered at x0 with stdev sigma """
    return( np.exp(-0.5*(x-x0)**2/sigma**2) / (np.sqrt(2*np.pi)*sigma) )

def timeavg( qty, ts, fti=None, ftf=None, step = 1, sigma = None ):
    cond = (ts<ftf)*(ts>=fti)
    if sigma is None:
        #use masked array to remove any stray NaN's
        qtycond = np.ma.masked_array(qty[cond],np.isnan(qty[cond]))
        qtycond = qtycond[::step]
        qtyavg = qtycond.mean(axis=0,dtype=np.float64)
    else:
        qtym = np.ma.masked_array(qty,np.isnan(qty))
        qtyavg=np.zeros_like(qtym)
        for (i0,t0) in enumerate(ts):
            mygauss_at_t0 = gaussf( ts, t0, sigma )
            #assumes uniform spacing in time
            if qtym.ndim == 2:
                mygauss_at_t0 = mygauss_at_t0[:,None]
            if t0 >= fti and t0 <= ftf:
                qtyavg[i0] += (qtym * mygauss_at_t0)[cond].sum(axis=0) / mygauss_at_t0[cond].sum()
        qtyavg[ts<fti] += qtyavg[ts>=fti][0]
        qtyavg[ts>ftf] += qtyavg[ts<=ftf][-1]
        #pdb.set_trace()

    return( qtyavg )


#
# 2D averaging
#

def mkavg2dnew(deltat = 100,endn = -1):
    if len(sys.argv[2:])>=2 and sys.argv[2].isdigit() and sys.argv[3].isdigit():
        whichi = int(sys.argv[2])
        whichn = int(sys.argv[3])
        if len(sys.argv[2:])==3:
            if sys.argv[4].isdigit():
                deltat = int(sys.argv[4])
        if sys.argv[1] == "mk2davg":
            postprocess2d(whichi = whichi, whichn = whichn, deltat = deltat, endn = endn)
        elif sys.argv[1] == "mk2dmrg":
            whichn1 = whichi
            whichn2 = whichn
            mrg2dnew(whichn1,whichn2,deltat = deltat)
    else:
        print("Syntax error")

def mrg2dnew(whichn1,whichn2,deltat=100):
    v = {}
    for n in range(whichn1,whichn2+1):
        f = "avg2d_%05d_%g.npz" % (n,deltat)
        print( "Reading %s..." % f )
        v = mrg2dnew_f(f,v)
    np.savez( "avg2d_%05d_%05d_%g.npz" % (whichn1, whichn2, deltat), **v )

def mrg2dnew_f(ft, v={}):
    vt = np.load(ft)
    for key in vt:
        #leave "t" for last so can properly average
        if key == "t":
            continue
        if key not in v:
            v[key]=vt[key]
        else:
            #properly average (taking into account the number of elements)
            lentot = len(v["t"])+len(vt["t"])
            w = 1.*len(v["t"])/lentot
            wt = 1.*len(vt["t"])/lentot
            v[key] = v[key]*w + vt[key]*wt
    #now deal with "t"
    if "t" not in v:
        v["t"] = []
    v["t"] = np.concatenate((v["t"],vt["t"]))
    vt.close()
    return v

def get_fieldline_time(fname):
    fin = open( dumps_path + fname, "rb" )
    header = fin.readline().split()
    #time of the dump
    t = myfloat(np.float64(header[0]))
    fin.close()
    return t
        
def postprocess2d(startn=0,endn=-1,whichi=0,whichn=1,**kwargs):
    rg("gdump")
    flist1 = np.sort(glob.glob( os.path.join(dumps_path, "dump[0-9][0-9][0-9]") ) )
    flist2 = np.sort(glob.glob( os.path.join(dumps_path, "dump[0-9][0-9][0-9][0-9]") ) )
    flist1.sort()
    flist2.sort()
    flist = np.concatenate((flist1,flist2))
    if len(flist) == 0:
        flist1 = np.sort(glob.glob( os.path.join(dumps_path, "dump[0-9][0-9][0-9]_0000") ) )
        flist2 = np.sort(glob.glob( os.path.join(dumps_path, "dump[0-9][0-9][0-9][0-9]_0000") ) )
        flist1.sort()
        flist2.sort()
        flist = np.concatenate((flist1,flist2))
    #default time interval for one averaging interval
    deltat = kwargs.pop("deltat",100)
    #number of intervals in a simulation
    try: 
        tlast = get_fieldline_time("../"+flist[-1])
    except IOError as e:
        print( "While reading %s:" % flist[-1] )
        print( "I/O error({0}): {1}".format(e.errno, e.strerror) )
        print( "Skipping" )
        return
    #round to smallest integer, i.e., ignore the last incomplete time averaging interval
    print("Last file name %s, last time %g" % (flist[-1], tlast))
    totintervals = np.int64(tlast/deltat) 
    #rg("gdump")
    print("Doing every %d interval out of %d total intervals..." % (whichn,totintervals))
    for numinterval in range(totintervals):
        #do every whichn'th averaging interval starting with whichi'th snapshot
        if numinterval < startn:
            continue
        if endn>=0 and numinterval >= endn:
            break
        if numinterval % whichn != whichi:
            continue
        #clear out the receptacle into which average will be stored
        v = {}
        v["t"] = []
        avgfname = "avg2d_%05d_%g.npz" % (numinterval, deltat)
        #start and end times of the interval
        tstartavg = numinterval     * deltat
        tendavg   = (numinterval+1) * deltat
        print( "Doing interval %d, t = [%g,%g)..." % (numinterval, tstartavg, tendavg) )
        if os.path.isfile( avgfname ):
            print( "File %s exists, skipping..." % avgfname )
            continue     
        for fldname in flist:
            #this only reads the first line in the file, so it is fast
            try: 
                t = get_fieldline_time("../"+fldname)
            except IOError as e:
                print( "While reading %s:" % fldname )
                print( "I/O error({0}): {1}".format(e.errno, e.strerror) )
                print( "Skipping" )
                continue
            #fieldline file falls outside of the averaging interval? skip it
            if t < tstartavg or t >= tendavg:
                #print("%s: t = %g falls outside interval, skipping" % (fldname, t))
                continue
            print("%s: t = %g falls inside interval [%g,%g) reading..." % (fldname, t, tstartavg, tendavg))
            sys.stdout.flush()
            try: 
                #remove the trailing _0000
                rd("../"+fldname.split("_")[0])
            except IOError as e:
                print( "While reading %s:" % fldname )
                print( "I/O error({0}): {1}".format(e.errno, e.strerror) )
                print( "Skipping" )
                continue
            cvellite()
            sys.stdout.flush()
            valdic = compvals2d()
            for key in valdic:
                #incorrect? should not assign an empty array?
                #if key not in v: v[key] = np.array([],dtype=np.float32())
                #collect all values of t averaged over
                if key == "t":
                    if key not in v:
                        v[key] = []
                    #append values of t that were averaged over
                    v[key].append(valdic[key])
                else:
                    if key not in v:
                        v[key] = np.array(valdic[key])
                    else:
                        v[key] += np.array(valdic[key])
        #obtain the time average by dividing by the number of files in the interval
        numfiles = np.float32(len(v["t"]))
        if numfiles != 0:
            #average all variables except time (times are stacked)
            for key in v:
                if key != "t":  v[key] /= numfiles
            print("Saving the average to %s..." % avgfname)
            np.savez(avgfname, **v)

def cvellite():
    global ud, bsq, bu, bd
    ud = mdot(gv3,uu)                  #g_mn u^n
    bu=np.empty_like(uu)              #allocate memory for bu
    #set component per component
    bu[0]=mdot(B[1:4], ud[1:4])             #B^i u_i
    bu[1:4]=(B[1:4] + bu[0]*uu[1:4])/uu[0]  #b^i = (B^i + b^t u^i)/u^t
    bd=mdot(gv3,bu)
    bsq=mdot(bu,bd)
    if "uradu" in globals():
        global uradd
        uradd = mdot(gv3,uradu)

def faradaylite():
    global omegaf1, omegaf2, omegaf1b, omegaf2b
    # these 2 are equal in degen electrodynamics when d/dt=d/dphi->0
    omegaf1=fFdd(0,1)/fFdd(1,3) # = ftr/frp
    omegaf2=fFdd(0,2)/fFdd(2,3) # = fth/fhp
    #
    # from jon branch, 04/10/2012
    #
    if 1:
        B1hat=B[1]*np.sqrt(gv3[1,1])
        B2hat=B[2]*np.sqrt(gv3[2,2])
        B3nonhat=B[3]
        v1hat=uu[1]*np.sqrt(gv3[1,1])/uu[0]
        v2hat=uu[2]*np.sqrt(gv3[2,2])/uu[0]
        v3nonhat=uu[3]/uu[0]
        #
        aB1hat=np.fabs(B1hat)
        aB2hat=np.fabs(B2hat)
        av1hat=np.fabs(v1hat)
        av2hat=np.fabs(v2hat)
        #
        vpol=np.sqrt(av1hat**2 + av2hat**2)
        Bpol=np.sqrt(aB1hat**2 + aB2hat**2)
        #
        # assume field swept back so omegaf is always larger than vphi (only true for outflow, so put in sign switch for inflow as relevant for disk near BH or even jet near BH)
        # GODMARK: These assume rotation about z-axis
        omegaf2b=np.fabs(v3nonhat) + np.sign(uu[1])*(vpol/Bpol)*np.fabs(B3nonhat)
        #
        omegaf1b=v3nonhat - B3nonhat*(v1hat*B1hat+v2hat*B2hat)/(B1hat**2+B2hat**2)

            
def fFdd(i,j):
    if i==0 and j==1:
        fdd =  gdet*(uu[2]*bu[3]-uu[3]*bu[2]) # f_tr
    elif i==1 and j==0:
        fdd = -gdet*(uu[2]*bu[3]-uu[3]*bu[2]) # -f_tr
    elif i==0 and j==2:
        fdd =  gdet*(uu[3]*bu[1]-uu[1]*bu[3]) # f_th
    elif i==2 and j==0:
        fdd = -gdet*(uu[3]*bu[1]-uu[1]*bu[3]) # -f_th
    elif i==0 and j==3:
        fdd =  gdet*(uu[1]*bu[2]-uu[2]*bu[1]) # f_tp
    elif i==3 and j==0:
        fdd = -gdet*(uu[1]*bu[2]-uu[2]*bu[1]) # -f_tp
    elif i==1 and j==3:
        fdd =  gdet*(uu[2]*bu[0]-uu[0]*bu[2]) # f_rp = gdet*B2
    elif i==3 and j==1:
        fdd = -gdet*(uu[2]*bu[0]-uu[0]*bu[2]) # -f_rp = gdet*B2
    elif i==2 and j==3:
        fdd =  gdet*(uu[0]*bu[1]-uu[1]*bu[0]) # f_hp = gdet*B1
    elif i==3 and j==2:
        fdd = -gdet*(uu[0]*bu[1]-uu[1]*bu[0]) # -f_hp = gdet*B1
    elif i==1 and j==2:
        fdd =  gdet*(uu[0]*bu[3]-uu[3]*bu[0]) # f_rh = gdet*B3
    elif i==2 and j==1:
        fdd = -gdet*(uu[0]*bu[3]-uu[3]*bu[0]) # -f_rh = gdet*B3
    else:
        fdd = np.zeros_like(uu[0])
    return fdd

delta = lambda kapa,nu: (kapa==nu)
fTudEM = lambda kapa,nu: bsq*uu[kapa]*ud[nu] + 0.5*bsq*delta(kapa,nu) - bu[kapa]*bd[nu]
fTudMA = lambda kapa,nu: (rho+gam*ug)*uu[kapa]*ud[nu]+(gam-1)*ug*delta(kapa,nu)
fTud = lambda kapa,nu: fTudEM(kapa,nu) + fTudMA(kapa,nu)
fRud = lambda kapa,nu: 4./3.*Erf*uradu[kapa]*uradd[nu]+1./3.*Erf*delta(kapa,nu)

def compvals2d():
    SMALL = 1e-20
    #cvellite()
    faradaylite()
    #returns angle-averaged values for the currently loaded dump
    dic = {}
    #
    #quantities
    dic["t"]=t
    dic["rho"]=rho.mean(-1)[:,:,None]
    dic["ug"]=ug.mean(-1)[:,:,None]
    dic["bsq"]=bsq.mean(-1)[:,:,None]
    #start betas
    dic["beta"]=(2*(gam-1)*ug/(bsq+SMALL)).mean(-1)[:,:,None]
    dic["ibeta"]=(bsq/(2*(gam-1)*ug)).mean(-1)[:,:,None]
    dic["rhobeta"]=(rho*2*(gam-1)*ug/(bsq+SMALL)).mean(-1)[:,:,None]
    dic["rhoibeta"]=(rho*bsq/(2*(gam-1)*ug)).mean(-1)[:,:,None]
    #end betas
    if "uradu" in globals():
        enth=1+(ug*gam+urad/3.)/rho
    else:
        enth=1+ug*gam/rho
    dic["unb"]=(enth*ud[0]).mean(-1)[:,:,None]
    dic["uu"]=uu.mean(-1)[:,:,:,None]
    dic["bu"]=bu.mean(-1)[:,:,:,None]
    dic["ud"]=ud.mean(-1)[:,:,:,None]
    dic["bd"]=bd.mean(-1)[:,:,:,None]
    #cell-centered magnetic field
    dic["B"]=B.mean(-1)[:,:,:,None]
    #face-centered magnetic field
    #
    dic["omegaf1"]=omegaf1.mean(-1)[:,:,None]
    dic["absomegaf1"]=np.abs(omegaf1).mean(-1)[:,:,None]
    dic["omegaf2"]=omegaf2.mean(-1)[:,:,None]
    dic["absomegaf2"]=np.abs(omegaf2).mean(-1)[:,:,None]
    dic["omegaf1b"]=omegaf1b.mean(-1)[:,:,None]
    dic["absomegaf1b"]=np.abs(omegaf1b).mean(-1)[:,:,None]
    dic["omegaf2b"]=omegaf2b.mean(-1)[:,:,None]
    dic["absomegaf2b"]=np.abs(omegaf2b).mean(-1)[:,:,None]
    dic["absbu"]=np.abs(bu).mean(-1)[:,:,:,None]
    dic["absbd"]=np.abs(bd).mean(-1)[:,:,:,None]
    dic["absuu"]=np.abs(uu).mean(-1)[:,:,:,None]
    dic["absud"]=np.abs(ud).mean(-1)[:,:,:,None]
    #
    dic["rhouu"]=(rho*uu).mean(-1)[:,:,:,None]
    dic["rhobu"]=(rho*bu).mean(-1)[:,:,:,None]
    dic["rhoud"]=(rho*ud).mean(-1)[:,:,:,None]
    dic["rhobd"]=(rho*bd).mean(-1)[:,:,:,None]
    dic["uguu"]=(ug*uu).mean(-1)[:,:,:,None]
    dic["ugud"]=(ug*ud).mean(-1)[:,:,:,None]
    #
    # electron-related diagnostic
    #
    #heating
    dic["qdot"] = qdot.mean(-1)[:,:,None]
    #heat fluxes
    dic["qisosq"] = qisosq.mean(-1)[:,:,None]
    dic["qisodotb"] = qisodotb.mean(-1)[:,:,None]
    dic["phi"] = phi.mean(-1)[:,:,None]
    #temperatures
    dic["Teldis"] = Teldis.mean(-1)[:,:,None]
    dic["Tel4"] = Tel4.mean(-1)[:,:,None]
    dic["Tel5"] = Tel5.mean(-1)[:,:,None]
    dic["Ttot"] = Ttot.mean(-1)[:,:,None]
    dic["rhoTeldis"] = (rho*Teldis).mean(-1)[:,:,None]
    dic["rhoTel4"] = (rho*Tel4).mean(-1)[:,:,None]
    dic["rhoTel5"] = (rho*Tel5).mean(-1)[:,:,None]
    dic["rhoTtot"] = (rho*Ttot).mean(-1)[:,:,None]
    #internal energies
    dic["ugeldis"] = ugeldis.mean(-1)[:,:,None]
    dic["ugel4"] = ugel4.mean(-1)[:,:,None]
    dic["ugel5"] = ugel5.mean(-1)[:,:,None]
    dic["terat"] = (ugel4/ugeldis-1.).mean(-1)[:,:,None]
    dic["teratinv"] = (ugeldis/ugel4-1.).mean(-1)[:,:,None]
    dic["feldis"] = (felcalc(ugeldis)).mean(-1)[:,:,None]
    dic["fel4"] = (felcalc(ugel4)).mean(-1)[:,:,None]
    dic["fel5"] = (felcalc(ugel5)).mean(-1)[:,:,None]
    dic["rhofeldis"] = (rho*felcalc(ugeldis)).mean(-1)[:,:,None]
    dic["rhofel4"] = (rho*felcalc(ugel4)).mean(-1)[:,:,None]
    dic["rhofel5"] = (rho*felcalc(ugel5)).mean(-1)[:,:,None]
    #properly compute average
    uuud=odot(uu,ud)[:,:,:,:,None]
    # part1: rho u^m u_l
    dic["rhouuud"]=(rho[:,:,None]*uuud).mean(-1)
    # part2: u u^m u_l
    dic["uguuud"]=(ug[:,:,None]*uuud).mean(-1)
    # part3: b^2 u^m u_l
    dic["bsquuud"]=(bsq[:,:,None]*uuud).mean(-1)
    # part6: b^m b_l
    dic["bubd"]=odot(bu,bd)[:,:,:,:,None].mean(-1)
    # u^m u_l
    dic["uuud"]=uuud.mean(-1)
    #
    #energy fluxes and faraday
    #EM/MA
    TudEMavgphi = np.zeros((4,4,nx,ny,1),dtype=np.float32)
    TudMAavgphi = np.zeros((4,4,nx,ny,1),dtype=np.float32)
    Tudavgphi = np.zeros((4,4,nx,ny,1),dtype=np.float32)
    Fddavgphi = np.zeros((4,4,nx,ny,1),dtype=np.float32)
    #
    # Radiation
    #
    if "uradu" in globals():
        dic["urad"] = urad.mean(-1)[:,:,None]
        dic["uradu"] = uradu.mean(-1)[:,:,:,None]
        Rudavgphi = np.zeros((4,4,nx,ny,1),dtype=np.float32)
        isunb=(-(1+(ug*gam+urad/3.)/rho)*ud[0]>1.0)
        isbnd=1-isunb
        #average energy transport velocity
        dic["uvE"] = (gdet*uu*(ug+urad)*isbnd).mean(-1)[...,None]
        dic["uE"]  = (gdet*   (ug+urad)*isbnd).mean(-1)[...,None]
    #to save memory use, average out each component in phi separately
    for i in range(4):
        for j in range(4):
            TudEMavgphi[i,j] = fTudEM(i,j).mean(-1)[...,None]
            TudMAavgphi[i,j] = fTudMA(i,j).mean(-1)[...,None]
            Tudavgphi[i,j] = fTud(i,j).mean(-1)[...,None]
            Fddavgphi[i,j] = fFdd(i,j).mean(-1)[...,None]
            if "uradu" in globals():
                Rudavgphi[i,j] = fRud(i,j).mean(-1)[...,None]
    dic["TudEM"]=TudEMavgphi
    dic["TudMA"]=TudMAavgphi
    dic["Tud"]=Tudavgphi
    if "uradu" in globals():
        dic["Rud"]=Rudavgphi
    dic["Tud"]=Tudavgphi
    dic["fdd"]=Fddavgphi
    #mu,sigma
    dic["mu"]= (-fTud(1,0)/(rho*uu[1])).mean(-1)[:,:,None]
    dic["sigma"]= (-fTudEM(1,0)/fTudMA(1,0)).mean(-1)[:,:,None]
    dic["bsqorho"]= (bsq/rho).mean(-1)[:,:,None]
    dic["absB"]= np.abs(B[1:4]).mean(-1)[:,:,:,None]
    aphi = psicalc()
    dic["aphisq"]= (aphi**2)[:,:,None]
    dic["bsquu"]= (bsq*uu).mean(-1)[:,:,:,None]
    dic["Bd3"]= (bd[3]*ud[0]-bd[0]*ud[3]).mean(-1)[:,:,None]
    dic["absBd3"]= np.abs(bd[3]*ud[0]-bd[0]*ud[3]).mean(-1)[:,:,None]
    return dic

def odot(a,b):
    """ Outer product of two vectors a^mu b_nu"""
    #the shape of the product is (4,4,nx,ny,max(a.nz,b.nz))
    outer_product = np.zeros(np.concatenate((np.array((4,4)),amax(a[0].shape,b[0].shape))),dtype=np.float32,order='F')
    for mu in np.arange(4):
        for nu in np.arange(4):
            outer_product[mu,nu] = a[mu]*b[nu]
    return(outer_product)


def amax(arg1,arg2):
    return(np.maximum(arg1,arg2))
    # arr1 = np.array(arg1)
    # arr2 = np.array(arg2)
    # #create storage array of size that's largest of arr1 and arr2
    # ret=np.zeros_like(arr1+arr2)
    # ret[arr1>=arr2]=arr1[arr1>=arr2]
    # ret[arr2>arr1]=arr2[arr2>arr1]
    # return(ret)
def amin(arg1,arg2):
    return(np.minimum(arg1,arg2))
    # arr1 = np.array(arg1)
    # arr2 = np.array(arg2)
    # #create storage array of size that's largest of arr1 and arr2
    # ret=np.zeros_like(arr1+arr2)
    # ret[arr1<=arr2]=arr1[arr1<=arr2]
    # ret[arr2<arr1]=arr2[arr2<arr1]
    # return(ret)

#############################
#
# Movie making
#
#############################

def mknewmov(**kwargs):
    startn = kwargs.pop("startn", 0)
    endn = kwargs.pop("endn",-1)
    dosavefig = kwargs.setdefault("dosavefig",1)
    dostreamlines = kwargs.setdefault("dostreamlines",1)
    ncont = kwargs.setdefault("ncont",200)
    maxaphi = kwargs.setdefault("maxaphi",1000)
    if len(sys.argv[2:])>=2 and sys.argv[2].isdigit() and sys.argv[3].isdigit():
        whichi = int(sys.argv[2])
        whichn = int(sys.argv[3])
        if len(sys.argv[2:])==3:
            dn = int(sys.argv[4])
        else:
            dn = 1
    else:
        print( "Usage: %s %s <whichi> <whichn> [<dn>]" % (sys.argv[0], sys.argv[1]) )
        return
    rg("gdump")
    v = np.load("qty.npz")
    flist = get_sorted_file_list()
    for fldname in flist:
        #find the index of the file
        fldindex = np.int(fldname.split("p")[-1])
        if fldindex < startn:
            continue
        if endn>=0 and fldindex >= endn:
            break
        if fldindex % dn: continue
        if fldindex/dn % (whichn) != whichi:
            #do every whichn'th snapshot starting with whichi'th snapshot
            continue
        if dosavefig:
            fname = "frame%04d.png"%fldindex
            if os.path.isfile( fname ):
                print("File %s exists, skipping..." % fname)
                continue
        print( "Reading " + fldname + " ..." )
        sys.stdout.flush()
        rd("../"+fldname)
        sys.stdout.flush()
        plt.clf()
        mkmfnew(v,findex=fldindex,doreload=0,**kwargs)

    
#new movie frame
def mkmfnew(v,**kwargs): #vmin=-6,vmax=0.5625):
    iti = kwargs.get("iti",None)
    itf = kwargs.get("itf",None)
    fti = kwargs.get("fti",None)
    ftf = kwargs.get("ftf",None)
    dic = get_titf(iti,itf,fti,ftf)
    #copy extra parameters over
    for key in dic: kwargs[key] = dic[key]
    findex = kwargs.get("findex",10000)
    maxsBphi = kwargs.get("maxsBphi",3)
    vmin = kwargs.get("vmin",-9)
    vmax = kwargs.get("vmax",-3)
    dostreamlines = kwargs.get("dostreamlines",1)
    ncont = kwargs.get("ncont",100)
    maxaphi = kwargs.get("maxaphi",1000)
    dovarylw = kwargs.get("dovarylw",2)
    domakeframes = kwargs.get("domakeframes",1)
    doreload = kwargs.get("doreload",1)
    plotlen = kwargs.get("plotlen",25)
    fntsize = kwargs.get("fntsize",16)
    dosavefig = kwargs.get("dosavefig",1)
    bbox_props = dict(boxstyle="round,pad=0.1", fc="w", ec="w", alpha=0.)
    plt.figure(0, figsize=(12,9), dpi=100)
    plt.clf()
    plot_ts_panels(v,**kwargs)
    #Rz xy
    gs1 = GridSpec(1, 1)
    gs1.update(left=0.07, right=0.45, top=0.995, bottom=0.5, wspace=0.05)
    #gs1.update(left=0.05, right=0.45, top=0.99, bottom=0.45, wspace=0.05)
    ax1 = plt.subplot(gs1[:, -1])
    if domakeframes:
        if doreload: rd("dump%03d" % findex)
        mkframe("lrho%04d_Rz%g" % (findex,plotlen),len=plotlen,ax=ax1,cb=False,pt=False,maxsBphi=maxsBphi,whichr=1.5,domask=0.5,vmin=vmin,vmax=vmax,dostreamlines=dostreamlines,ncont=ncont,maxaphi=maxaphi,dovarylw=dovarylw) #domask = 0.5 is important so that magnetic field lines extend down all the way to BH
    ax1.set_ylabel(r'$z\ [r_g]$',fontsize=16,ha='center')
    ax1.set_xlabel(r'$x\ [r_g]$',fontsize=16)
    for label in ax1.get_xticklabels() + ax1.get_yticklabels():
        label.set_fontsize(fntsize)
    placeletter(ax1,"$(\mathrm{a})$",va="center",bbox=bbox_props)
    gs2 = GridSpec(1, 1)
    gs2.update(left=0.5, right=1, top=0.995, bottom=0.5, wspace=0.05)
    ax2 = plt.subplot(gs2[:, -1])
    if domakeframes:
        mkframexy("lrho%04d_xy%g" % (findex,plotlen),len=plotlen,ax=ax2,cb=True,pt=False,dostreamlines=True,dovarylw=dovarylw,domask=0.5,vmin=vmin,vmax=vmax) #,label=r"$\log\rho$",fontsize=20) #domask = 0.5 is important so that magnetic field lines extend down all the way to BH
    ax2.set_ylabel(r'$y\ [r_g]$',fontsize=16,ha='center',labelpad=0)
    ax2.set_xlabel(r'$x\ [r_g]$',fontsize=16)
    for label in ax2.get_xticklabels() + ax2.get_yticklabels():
        label.set_fontsize(fntsize)
    placeletter(ax2,"$(\mathrm{b})$",va="center",bbox=bbox_props)
    if dosavefig:
        plt.savefig("frame%04d.png" % findex,dpi=120)

def get_titf(iti=None, itf=None, fti=None, ftf=None):
    if iti is None or itf is None or fti is None or ftf is None or os.path.isfile(os.path.join("titf.txt")):
        dic = {}
        gd1 = np.loadtxt( "titf.txt",
                          dtype=np.float64, 
                          unpack = True )
        dic["iti"] = gd1[0]
        dic["itf"] = gd1[1]
        dic["fti"] = gd1[2]
        dic["ftf"] = gd1[3]
        sys.stdout.write( "Found titf.txt: iti = %g, itf = %g, fti = %g, ftf = %g"
                   % (dic["iti"],dic["itf"],dic["fti"],dic["ftf"]) )
        if len(gd1)>=6:
            dic["vmin"] = gd1[4]
            dic["vmax"] = gd1[5]
            sys.stdout.write( ", vmin = %g, vmax = %g" % (dic["vmin"], dic["vmax"]) )
        if len(gd1)>=9:
            dic["mdotmax"] = gd1[6]
            dic["phibhmax"] = gd1[7]
            dic["etabhmax"] = gd1[8]
            sys.stdout.write( ", mdotmax = %g, phibhmax = %g, etabhmax = %g" % (dic["mdotmax"], dic["phibhmax"], dic["etabhmax"]) )
        if len(gd1)>=11:
            dic["ncont"] = gd1[9]
            dic["maxaphi"] = gd1[10]
            sys.stdout.write( ", ncont = %g, aphimax = %g" % (dic["ncont"], dic["maxaphi"]) )
        if len(gd1)>=12:
            dic["plotlen"] = gd1[11]
            sys.stdout.write( ", plotlen = %g" % (dic["plotlen"]) )
        if len(gd1)>=14:
            dic["thetaemin"] = gd1[12]
            dic["thetaemax"] = gd1[13]
            sys.stdout.write( ", thetaemin = %g, thetaemax = %g" % (dic["thetaemin"],dic["thetaemax"]) )
        sys.stdout.write("\n")
    else:
        dic["iti"] = 3500
        dic["itf"] = 9500
        dic["fti"] = 3500
        dic["ftf"] = 9500
        print( "Warning: titf.txt not found: using default numbers for averaging: iti = %g, itf = %g, fti = %g, ftf = %g" % (dic["iti"],dic["itf"],dic["fti"],dic["ftf"]) )
    return dic
    

def plot_ts_panels(v,findex=10000,
            sigma=1500,sigma1=None,prefactor=100,domakeframes=1,maxsBphi=3,
            doreload=1,dosavefig = 1,fntsize=16,myi=None,vmin=-9,vmax=-3,
            dostreamlines=1,dovarylw=2,**kwargs):
    global FMavg, t
    iti = kwargs.pop("iti",2000)
    itf = kwargs.pop("itf",5000)
    fti = kwargs.pop("fti",2000)
    ftf = kwargs.pop("ftf",5000)
    mdotmax = kwargs.pop("mdotmax",None)
    phibhmax = kwargs.pop("phibhmax",None)
    etabhmax = kwargs.pop("etabhmax",None)
    ncont = kwargs.pop("ncont",100)
    maxaphi = kwargs.pop("maxaphi",1000)
    plotlen = kwargs.pop("plotlen",25)
    top = kwargs.pop("top",0.42)
    left = kwargs.pop("left",0.1)
    right = kwargs.pop("right",0.94)
    bottom = kwargs.pop("bottom",0.06)
    wspace = kwargs.pop("wspace",0.01)
    hspace = kwargs.pop("hspace",0.04)
    xlabelpad =kwargs.pop("xlabelpad",0)
    doannotate=kwargs.pop("doannotate",1)
    #plt.clf()
    if myi is None:
        myi = np.sum(v["rvals"]<5)
    hori = np.sum(v["rvals"]<=rhor)
    tmin=max(v["t"][0],0)
    tmax=min(v["t"][-1],max(fti,ftf))
    print( "tmin = %g, tmax = %g" % (tmin, tmax) )
    print( "Radius at which Mdot and powers are evaluated = %g" % v["rvals"][myi] )
    print( "Radius at which PhiBH is evaluated = %g = %g*rhor" % (v["rvals"][hori], v["rvals"][hori]/rhor) )
    #mdot,pjet,pjet/mdot plots
    gs3 = GridSpec(3, 1)
    gs3.update(left=left, right=right, top=top, bottom=bottom, wspace=wspace, hspace=hspace)
    #
    #mdot plot
    #
    ax31 = plt.subplot(gs3[-3,:])
    ax31.set_ylabel(r'$\dot Mc^2$',fontsize=fntsize,ha="left") #,labelpad=9)
    plt.setp( ax31.get_xticklabels(), visible=False)
    which = (v["t"] <= ftf)
    #start plotting
    ind = np.where(v["t"]==t)
    tt = v["t"]
    ax31.plot(v["t"][ind],v["FM"][ind,myi],'o',mfc='r')
    if 'FMavg' not in globals() and sigma is not None:
        FMavg = 0*v["t"]
        FMavg[which] = timeavg(v["FM"][which,myi],v["t"][which],fti=iti,ftf=ftf,sigma=sigma)
        FMavg[~which]  = FMavg[~which] + FMavg[which][-1]      
    elif sigma is None:
        FMavg1 = 0*v["t"]
        FMavg2 = 0*v["t"]
        FMavg1 = timeavg(v["FM"][which,myi],v["t"][which],fti=iti,ftf=itf,sigma=sigma)+0*v["t"][which]
        FMavg2 = timeavg(v["FM"][which,myi],v["t"][which],fti=fti,ftf=ftf,sigma=sigma)+0*v["t"][which]
        FMavg = np.copy(FMavg1)
        FMavg[(iti<=t)*(t<=itf)] = FMavg1[(iti<=t)*(t<=itf)]
        FMavg[(fti<=t)*(t<=ftf)] = FMavg2[(fti<=t)*(t<=ftf)]
    ax31.plot(v["t"][which],v["FM"][which,myi],"k")
    #pdb.set_trace()
    #ensure of the same shape as the rest
    l, = ax31.plot(tt[which],FMavg[which],"k")
    l.set_dashes([10,5])
    #end plotting
    ax31.set_xlim(tmin,tmax)
    ymax=ax31.get_ylim()[1]
    ymed=np.median(v["FM"][which,myi])
    if ymax > 1:
        ymax=2*(np.floor(np.floor(ymax+1.5)/2))
    # if ymax > 5*ymed:
    #     ymax = 2*ymed
    ymax = float("{0:.2g}".format(ymax))
    if mdotmax is not None:
        ymax = np.abs(mdotmax)
    print( "max(FM) = %g" % ymax )
    # ymax = 26.
    if mdotmax is not None:
        if mdotmax > 0:
            ax31.set_yticks((ymax/2.,ymax))
        else:
            ax31.set_ylim((ymax/100.,ymax))
            ax31.set_yscale("log")
    ax31.grid(True)
    bbox_props = dict(boxstyle="round,pad=0.1", fc="w", ec="w", alpha=0.)
    #placeletter(ax31,"$(\mathrm{c})$",fx=0.15,fy=0.5,bbox=bbox_props)
    if doannotate:
        ax31.annotate("$(\mathrm{c})$",(0,1),xytext=(5,-5),xycoords="axes fraction",textcoords="offset points",fontsize=20,rotation=0,rotation_mode="anchor",ha="left",va="top",color="black")
    ax31r = ax31.twinx()
    if mdotmax is not None:
        if mdotmax > 0:
            ax31r.set_yticks((ymax/2.,ymax))
            ax31r.set_ylim(0,ymax)
            ax31.set_ylim(0,ymax)
        else:
            #ax31r.set_yticks((ymax/2.,ymax))
            ax31r.set_ylim(ymax/100.,ymax)
            ax31r.set_yscale("log")
            ax31.set_ylim(ymax/100.,ymax)
    for label in ax31.get_xticklabels() + ax31.get_yticklabels() + ax31r.get_yticklabels():
        label.set_fontsize(fntsize)
    mathify_axes_ticks(ax31,fontsize=fntsize)
    mathify_axes_ticks(ax31r,fontsize=fntsize)
    #
    #\phi plot
    #
    ax35 = plt.subplot(gs3[-2,:])
    #start plotting
    PhiBHcgs = v["PhiBH"][:,hori]*(4*np.pi)**0.5
    phibh = PhiBHcgs/(get_fracphi()*FMavg)**0.5
    ax35.plot(v["t"][which],phibh[which],"k")
    ax35.plot(v["t"][ind],phibh[ind],'o',mfc='r')
    phiavg1 = timeavg(phibh[:],v["t"],fti=iti,ftf=itf,sigma=sigma1)+0*v["t"]
    phiavg2 = timeavg(phibh[:],v["t"],fti=fti,ftf=ftf,sigma=sigma1)+0*v["t"]
    l, = ax35.plot(v["t"][(iti<=tt)*(tt<=itf)],phiavg1[(iti<=tt)*(tt<=itf)],"k")
    l.set_dashes([10,5])
    l, = ax35.plot(v["t"][(fti<=tt)*(tt<=ftf)],phiavg2[(fti<=tt)*(tt<=ftf)],"k")
    l.set_dashes([10,5])
    #end plotting
    ax35.set_xlim(tmin,tmax)
    ymax=ax35.get_ylim()[1]
    if phibhmax is not None:
        ymax = phibhmax
        ax35.set_ylim(0,ymax)
    if 1 < ymax and ymax < 2: 
        #ymax = 2
        tck=(1,)
        ax35.set_yticks(tck)
        #ax35.set_yticklabels(('','1','2'))
    elif ymax < 1: 
        ymax = 1
        tck=(0.5,1)
        ax35.set_yticks(tck)
        ax35.set_yticklabels(('','1'))
    else:
        ymax=np.floor(ymax)+1
        if ymax >= 60:
            tck=np.arange(1,ymax/30.)*30.
        elif ymax >= 20:
            tck=np.arange(1,ymax/10.)*10.
        elif ymax >= 10:
            tck=np.arange(1,ymax/5.)*5.
        else:
            tck=(np.int32(ymax)/2,2*(np.int32(ymax)/2))
        ax35.set_yticks(tck)
    print( "max(phi) = %g" % ymax )
    ax35.grid(True)
    plt.setp( ax35.get_xticklabels(), visible=False)
    #placeletter(ax35,"$(\mathrm{d})$",fx=0.15,fy=0.1,bbox=bbox_props)
    if doannotate:
        ax35.annotate("$(\mathrm{d})$",(0,1),xytext=(5,-5),xycoords="axes fraction",textcoords="offset points",fontsize=20,rotation=0,rotation_mode="anchor",ha="left",va="top",color="black")
    ax35.set_ylabel(r"$\phi$",size=fntsize) #,labelpad=22) #labelpad=25
    ax35.grid(True)
    ax35r = ax35.twinx()
    ax35r.set_ylim(ax35.get_ylim())
    #if phibhmax is None:
    ax35r.set_yticks(tck)
    for label in ax35.get_xticklabels() + ax35.get_yticklabels() + ax35r.get_yticklabels():
        label.set_fontsize(fntsize)
    mathify_axes_ticks(ax35,fontsize=fntsize)
    mathify_axes_ticks(ax35r,fontsize=fntsize)
    #
    #pjet/<mdot>
    #
    ax34 = plt.subplot(gs3[-1,:])
    #start plotting
    etabh = (v["FM"]-v["FE"])[:,myi]/FMavg
    #print( "FMavg = %g" % FMavg )
    ax34.plot(v["t"][which],etabh[which]*prefactor,"k")
    ax34.plot(v["t"][ind],etabh[ind]*prefactor,'o',mfc='r')
    etabhavg1 = timeavg(etabh[:],v["t"],fti=iti,ftf=itf,sigma=sigma1) + 0*v["t"]
    etabhavg2 = timeavg(etabh[:],v["t"],fti=fti,ftf=ftf,sigma=sigma1) + 0*v["t"]
    l, = ax34.plot(v["t"][(iti<=tt)*(tt<=itf)],etabhavg1[(iti<=tt)*(tt<=itf)]*prefactor,"k")
    l.set_dashes([10,5])
    l, = ax34.plot(v["t"][(fti<=tt)*(tt<=ftf)],etabhavg2[(fti<=tt)*(tt<=ftf)]*prefactor,"k")
    l.set_dashes([10,5])
    #end plotting
    ax34.set_xlim(tmin,tmax)    
    #ymax=ax34.get_ylim()[1]
    # if ymax == 0:
    #     print("Got max(etabh) = 0, recomputing...")
    #     ymax = etabh[v["t"]<ftf].nanmax()
    maxval = np.nanmax(ma.filled(etabh[v["t"]<ftf]))
    medval = np.median(ma.filled(etabh[v["t"]<ftf]))
    # if maxval > 5*medval:
    #     maxval = 2*medval
    if maxval > 0.5:
        ymax=np.floor(maxval+1)*prefactor
    else:
        ymax = maxval*prefactor #ax34.get_ylim()[1]
    ymax = float("{0:.2g}".format(ymax))    
    ax34.set_ylim(0,ymax)
    if etabhmax is not None:
        ymax = etabhmax
    print( "max(etabh) = %g" % ymax )
    if ymax > 50 and ymax < 200: 
        tck = np.arange(0,ymax,50)
    elif ymax >= 200:
        tck = np.arange(0,ymax,100)
        ax34.set_yticks(tck)
    else:
        tck = np.array([0,ymax/2,ymax])
    ax34.set_yticks(tck)
    #placeletter(ax34,"$(\mathrm{e})$",fx=0.15,fy=0.35,bbox=bbox_props)
    if doannotate:
        ax34.annotate("$(\mathrm{e})$",(0,1),xytext=(5,-5),xycoords="axes fraction",textcoords="offset points",fontsize=20,rotation=0,rotation_mode="anchor",ha="left",va="top",color="black")
    #reset lower limit to 0
    ax34.set_xlabel(r'$t\ [r_g/c]$',fontsize=fntsize,labelpad=xlabelpad)
    ax34.set_ylim(0,ymax)
    ax34.set_ylabel(r"$p$",size=fntsize) #,labelpad=12)
    ax34.grid(True)
    ax34r = ax34.twinx()
    ax34r.set_ylim(ax34.get_ylim())
    ax34r.set_yticks(tck)
    for label in ax34.get_xticklabels() + ax34.get_yticklabels() + ax34r.get_yticklabels():
        label.set_fontsize(fntsize)
    mathify_axes_ticks(ax34,fontsize=fntsize)
    mathify_axes_ticks(ax34r,fontsize=fntsize)
    
        
def mkframe(fname,ax=None,cb=True,vmin=None,vmax=None,len=20,ncell=800,pt=True,shrink=1,dostreamlines=True,downsample=4,density=2,dodiskfield=False,minlendiskfield=0.2,minlenbhfield=0.2,dovarylw=True,dobhfield=True,dsval=0.01,color='k',dorandomcolor=False,doarrows=True,lw=None,skipblankint=False,detectLoops=True,minindent=1,minlengthdefault=0.2,startatmidplane=True,showjet=False,arrowsize=1,startxabs=None,startyabs=None,populatestreamlines=True,useblankdiskfield=True,dnarrow=2,whichr=1.0,ncont=100,maxaphi=100,aspect=1.0,isnstar=False,kval=0,kvalvar=0,onlyeta=True,maxsBphi=None,domirror=True,nanout=False,whichvar="lrho",avgbsqorho=None,fntsize=None,aphiaccent=None,cmap=None,domask=1,**kwargs):
    #######################
    #setup PYTHONPATH because needed for streamlines: fstreamplot() function
    import os
    import sys
    sys.path.append(os.path.join(os.environ["HOME"],"py"))
    from streamlines import fstreamplot
    #######################
    extent=(-len,len,-len/aspect,len/aspect)
    if cmap is None:
        palette=cm.jet
    else:
        palette=cmap
    palette.set_bad('k', 1.0)
    #palette.set_over('r', 1.0)
    #palette.set_under('g', 1.0)
    if isnstar:
        domask = Rin
    else:
        domask = domask
    if avgbsqorho is None:
        avgbsqorho = lambda: rho
    if not isnstar:
        rhor=1+(1-a**2)**0.5
        ihor = iofr(rhor)
    else:
        rhor=1
        ihor = 0
        #a=1
    if 'rho' in globals():
        ilrho = reinterp(np.log10(rho),extent,ncell,domask=1.0,rhor=rhor,kval=kval)
    else:
        ilrho = None
    if True:
        aphi = psicalc()
        iaphi = reinterp(aphi,extent,ncell,domask=0,rhor=rhor,kval=kval)
        #maxabsiaphi=np.max(np.abs(iaphi))
        maxabsiaphi = maxaphi #50
        #ncont = 100 #30
        levs=np.linspace(-maxabsiaphi,maxabsiaphi,ncont)
        Br = dxdxp[1,1]*B[1]+dxdxp[1,2]*B[2]
        Bh = dxdxp[2,1]*B[1]+dxdxp[2,2]*B[2]
        #note toroidal field located at faces
        #Bp = gdetB[3]/gdet*dxdxp[3,3]
        if "gdetB" in globals():
            if 'Bstag' not in globals():
                print("Bstag is same as B, so will use gdetB/gdet to show perpendicular field component")
                Bp = gdetB[3]/gdet*dxdxp[3,3]
            #elif Bstag is B:
            else:
                Bp = Bstag[3]*dxdxp[3,3]
        else:
            Bp = B[3]*dxdxp[3,3]
        #Bp[(h<0)+(h>np.pi)] *= -1
        #Bp = Bstag[2]
        #
        Brnorm=Br
        Bhnorm=Bh*np.abs(r)
        Bpnorm=Bp*np.abs(r*np.sin(h))
        #
        Bznorm=Brnorm*np.cos(h)-Bhnorm*np.sin(h)
        BRnorm=Brnorm*np.sin(h)+Bhnorm*np.cos(h)
        if nanout:
            if not domirror:
                Bznorm[:,-1]*=NaN
                Bznorm[:,0]*=NaN
                BRnorm[:,-1]*=NaN
                BRnorm[:,0]*=NaN
                Bpnorm[:,-1]*=NaN
                Bpnorm[:,0]*=NaN
            else:
                NBND=4
                Bznorm[:,:NBND]*=NaN
                Bznorm[:,-NBND+1:]*=NaN
                BRnorm[:,:NBND]*=NaN
                BRnorm[:,-NBND+1:]*=NaN
                Bpnorm[:,:NBND]*=NaN
                Bpnorm[:,-NBND+1:]*=NaN
        #
        iBz = reinterp(Bznorm,extent,ncell,isasymmetric=False,domask=domask,rhor=rhor,kval=kval-0.5,domirror=domirror)
        iBR = reinterp(BRnorm,extent,ncell,isasymmetric=True,domask=domask,rhor=rhor,kval=kval-0.5,domirror=domirror) #isasymmetric = True tells to flip the sign across polar axis
        iBp = reinterp(Bpnorm,extent,ncell,isasymmetric=True,domask=domask,rhor=rhor,kval=kval,domirror=domirror) #isasymmetric = True tells to flip the sign         #
        if whichvar is not None:
            irho = reinterp(rho,extent,ncell,isasymmetric=False,domask=domask,rhor=rhor,kval=kval-0.5,domirror=domirror)
            iug  = reinterp(ug,extent,ncell,isasymmetric=False,domask=domask,rhor=rhor,kval=kval-0.5,domirror=domirror)
            ibsq = reinterp(bsq,extent,ncell,isasymmetric=False,domask=domask,rhor=rhor,kval=kval-0.5,domirror=domirror)
            iavgbsqorho = reinterp(avgbsqorho(),extent,ncell,isasymmetric=False,domask=domask,rhor=rhor,kval=kvalvar-0.5,domirror=domirror)
            if not isinstance(whichvar,basestring):
                var = whichvar()
                ivar = reinterp(var,extent,ncell,isasymmetric=False,domask=domask,rhor=rhor,kval=kvalvar-0.5,domirror=domirror)
        if 0 and dorandomcolor:
            Ba=np.copy(B)
            cond = (B[1]<0)
            Ba[2,cond]*=-1
            Ba[3,cond]*=-1
            Ba[1,cond]*=-1
            Bar = dxdxp[1,1]*Ba[1]+dxdxp[1,2]*Ba[2]
            Bah = dxdxp[2,1]*Ba[1]+dxdxp[2,2]*Ba[2]
            Bap = Ba[3]*dxdxp[3,3]
            #
            Barnorm=Bar
            Bahnorm=Bah*np.abs(r)
            Bapnorm=Bap*np.abs(r*np.sin(h))
            #
            Baznorm=Barnorm*np.cos(h)-Bahnorm*np.sin(h)
            BaRnorm=Barnorm*np.sin(h)+Bahnorm*np.cos(h)
            #
            iBaz = reinterp(Baznorm,extent,ncell,domask=0.8,rhor=rhor,kval=kval)
            iBaR = reinterp(BaRnorm,extent,ncell,isasymmetric=True,domask=0.8,rhor=rhor,kval=kval) #isasymmetric = True tells to flip the sign across polar axis
        else:
            iBaz = None
            iBaR = None
        if showjet:
            imu = reinterp(mu,extent,ncell,domask=0.8,rhor=rhor,kval=kval-0.5)
        #
        if dovarylw:
            if "urad" in globals():
                iibeta = reinterp(0.5*bsq/(gam-1)/(ug+urad),extent,ncell,domask=0,rhor=rhor,kval=kval-0.5)
            else:
                iibeta = reinterp(0.5*bsq/(gam-1)/(ug),extent,ncell,domask=0,rhor=rhor,kval=kval-0.5)
            ibsqorho = reinterp(bsq/rho,extent,ncell,domask=0,rhor=rhor,kval=kval-0.5)
            ibsqo2rho = 0.5 * ibsqorho
        xi = np.linspace(extent[0], extent[1], ncell)
        yi = np.linspace(extent[2], extent[3], ncell)
        #myspeed=np.sqrt(iBR**2+iBz**2)
    #
    #myslines=streamplot(ti[:,0,0],tj[0,:,0],avg_gdetB[0,:,:,0].transpose(),avg_gdetB[1,:,:,0].transpose(),density=2,downsample=4,linewidth=1)
    #for c in cset2.collections:
    #    c.set_linestyle('solid')
    #CS = plt.contourf(xi,yi,zi,15,cmap=palette)
    if ax == None:
        ax = plt.gca()
        # CS = plt.imshow(ilrho, extent=extent, cmap = palette, norm = colors.Normalize(clip = False),origin='lower',vmin=vmin,vmax=vmax)
        # if not dostreamlines:
        #     cset2 = plt.contour(iaphi,linewidths=0.5,colors='k', extent=extent,hold='on',origin='lower',levels=levs)
        # else:
        #     lw = 0.5+1*ftr(np.log10(amax(ibsqo2rho,1e-6+0*ibsqorho)),np.log10(1),np.log10(2))
        #     lw += 1*ftr(np.log10(amax(iibeta,1e-6+0*ibsqorho)),np.log10(1),np.log10(2))
        #     lw *= ftr(np.log10(amax(iibeta,1e-6+0*iibeta)),-3.5,-3.4)
        #     # if t < 1500:
        #     #     lw *= ftr(ilrho,-2.,-1.9)
        #     fstreamplot(yi,xi,iBR,iBz,density=2,downsample=4,linewidth=lw,ax=ax,detectLoops=True,dodiskfield=False,dobhfield=True,startatmidplane=True,a=a)
        #     #streamplot(yi,xi,iBR,iBz,density=3,linewidth=1,ax=ax)
        # plt.xlim(extent[0],extent[1])
        # plt.ylim(extent[2],extent[3])
    if whichvar == "lrho":
        CS = ax.imshow(ilrho, extent=extent, cmap = palette, norm = colors.Normalize(clip = False),origin='lower',vmin=vmin,vmax=vmax)
    if not isinstance(whichvar,basestring):
        CS = ax.imshow(ivar, extent=extent, cmap = palette, norm = colors.Normalize(clip = False),origin='lower',vmin=vmin,vmax=vmax)
    if whichvar == "Bphi":
        siBp = np.sqrt(np.abs(iBp))
        if maxsBphi is None:
            maxsiBp = np.max(siBp)
        else:
            maxsiBp = maxsBphi
        print( "Max(Sqrt(Bout)) = %g" % maxsiBp )
        CS = ax.imshow(np.sign(iBp)*siBp, extent=extent, cmap = palette, norm = colors.Normalize(clip = False),origin='lower', vmin=-0.3*maxsiBp,vmax=0.3*maxsiBp)
    if whichvar == 'bsqorho':
        CS = ax.imshow(ibsq/irho, extent=extent, cmap = palette, norm = colors.Normalize(clip = False),origin='lower', vmin=0,vmax=100)
    if whichvar == 'avgbsqorho':
        CS = ax.imshow(iavgbsqorho, extent=extent, cmap = palette, norm = colors.Normalize(clip = False),origin='lower', vmin=vmin,vmax=vmax)
    if whichvar == 'rho':
        CS = ax.imshow(irho, extent=extent, cmap = palette, norm = colors.Normalize(clip = False),origin='lower')
    if showjet:
        ax.contour(imu,linewidths=0.5,colors='g', extent=extent,hold='on',origin='lower',levels=(2,))
        ax.contour(iaphi,linewidths=0.5,colors='b', extent=extent,hold='on',origin='lower',levels=(aphi[ihor,ny/2,0],))
    if not dostreamlines:
        cset2 = ax.contour(iaphi,linewidths=0.5,colors='k', extent=extent,hold='on',origin='lower',levels=levs)
        # if aphiaccent is not None:
        #     ax.contour(iaphi,linewidths=2,colors='k', extent=extent,hold='on',origin='lower',levels=(aphiaccent,))
        traj = None
    elif dostreamlines == 1:
        if dovarylw:
            if False:
                #old way
                lw = 0.5+1*ftr(np.log10(amax(ibsqo2rho,1e-6+0*ibsqorho)),np.log10(1),np.log10(2))
                lw += 1*ftr(np.log10(amax(iibeta,1e-6+0*ibsqorho)),np.log10(1),np.log10(2))
                lw *= ftr(np.log10(amax(iibeta,1e-6+0*iibeta)),-3.5,-3.4)
                # if t < 1500:
                lw *= ftr(iaphi,0.001,0.002)
            elif dovarylw==1:
                #new way, to avoid glitches in u_g in jet region to affect field line thickness
                lw1 = 2*ftr(np.log10(amax(ibsqo2rho,1e-6+0*ibsqorho)),np.log10(1),np.log10(2))
                lw2 = ftr(np.log10(amax(iibeta,1e-6+0*ibsqorho)),np.log10(1),np.log10(2))
                lw = 0.5 + amax(lw1,lw2)
                #lw *= ftr(np.log10(amax(iibeta,1e-6+0*iibeta)),-3.5,-3.4)
                # if t < 1500:
                lw *= ftr(iaphi,0.001,0.002)
            elif dovarylw==2:
                #new way, to avoid glitches in u_g in jet region to affect field line thickness
                lw1 = 2*ftr(np.log10(amax(ibsqo2rho,1e-6+0*ibsqorho)),np.log10(1),np.log10(2))
                lw2 = ftr(np.log10(amax(iibeta,1e-6+0*ibsqorho)),np.log10(1),np.log10(2))
                lw = 0.5 + amax(lw1,lw2)
                #lw *= ftr(np.log10(amax(iibeta,1e-6+0*iibeta)),-3.5,-3.4)
                # if t < 1500:
                #lw *= ftr(iaphi,0.001,0.002)
            elif dovarylw==3:
                #new way, to avoid glitches in u_g in jet region to affect field line thickness
                lw1 = 2*ftr(np.log10(amax(ibsqo2rho,1e-6+0*ibsqorho)),np.log10(1),np.log10(2))
                lw2 = ftr(np.log10(amax(iibeta,1e-6+0*ibsqorho)),np.log10(1),np.log10(2))
                lw = 0.5 + amax(lw1,lw2)
            elif dovarylw==4:
                #new way, to avoid glitches in u_g in jet region to affect field line thickness
                lw1 = 2*ftr(np.log10(amax(ibsqorho,1e-6+0*ibsqorho)),np.log10(0.9),np.log10(1))
                val = 1. #/(6.*(gam-1.))
                lw2 = ftr(np.log10(amax(iibeta,1e-6+0*ibsqorho)),np.log10(0.9*val),np.log10(val))
                lw = 0.5 + amax(lw1,lw2)
        #pdb.set_trace()
        traj = fstreamplot(yi,xi,iBR,iBz,ua=iBaR,va=iBaz,density=density,downsample=downsample,linewidth=lw,ax=ax,detectLoops=detectLoops,dodiskfield=dodiskfield,dobhfield=dobhfield,startatmidplane=startatmidplane,a=a,minlendiskfield=minlendiskfield,minlenbhfield=minlenbhfield,dsval=dsval,color=color,doarrows=doarrows,dorandomcolor=dorandomcolor,skipblankint=skipblankint,minindent=minindent,minlengthdefault=minlengthdefault,arrowsize=arrowsize,startxabs=startxabs,startyabs=startyabs,populatestreamlines=populatestreamlines,useblankdiskfield=useblankdiskfield,dnarrow=dnarrow,whichr=whichr)
    elif dostreamlines == 2:
        quiver(yi, xi, iBR, iBz)
        #streamplot(yi,xi,iBR,iBz,density=3,linewidth=1,ax=ax)
    ax.set_xlim(extent[0],extent[1])
    ax.set_ylim(extent[2],extent[3])
    #CS.cmap=cm.jet
    #CS.set_axis_bgcolor("#bdb76b")
    if True == cb:
        cbar=plt.colorbar(CS,ax=ax,shrink=shrink) # draw colorbar
        if fntsize is not None:
            #set font size of colorbar tick labels
            cl = plt.getp(cbar.ax, 'ymajorticklabels')
            plt.setp(cl, fontsize=fntsize)
    #plt.title(r'$\log_{10}\rho$ at $t = %4.0f$' % t)
    if True == pt:
        plt.title('log rho at t = %4.0f' % t)
    #if None != fname:
    #    plt.savefig( fname + '.png' )
    if dostreamlines == 1:
        return(traj)

def mkframexy(fname,ax=None,cb=True,vmin=None,vmax=None,len=20,ncell=800,pt=True,shrink=1,dostreamlines=True,arrowsize=1,isnstar=False,avgbsqorho=None,whichvar=None,fntsize=None,aphiaccent=None,cmap=None,domask=1,**kwargs):
    #######################
    #setup PYTHONPATH because needed for streamlines: fstreamplot() function
    import os
    import sys
    sys.path.append(os.path.join(os.environ["HOME"],"py"))
    from streamlines import fstreamplot
    #######################
    extent=(-len,len,-len,len)
    if cmap is None:
        palette=cm.jet
    else:
        palette=cmap
    if isnstar:
        domask = Rin
    else:
        domask = domask
    if avgbsqorho is None:
        avgbsqorho = lambda: rho
    if not isnstar:
        rhor=1+(1-a**2)**0.5
        ihor = iofr(rhor)
    else:
        rhor=1
        ihor = 0
        #a=1
    label=kwargs.pop('label',None)
    fontsize=kwargs.pop('fontsize',None)
    lw=kwargs.pop('lw',None)
    dovarylw=kwargs.pop('dovarylw',0)
    dobhfield=kwargs.pop('dobhfield',True)
    density=kwargs.pop('density',1)
    downsample=kwargs.pop('downsample',1)
    thetarot=kwargs.pop('thetarot',0)
    dnarrow=kwargs.pop('dnarrow',0)
    kval=kwargs.pop('kval',None)
    minlengthdefault=kwargs.pop('minlengthdefault',0.2)
    if kval is not None:
        thetarot = ph[0,0,kval]
    palette.set_bad('k', 1.0)
    #palette.set_over('r', 1.0)
    #palette.set_under('g', 1.0)
    ilrho = reinterpxy(np.log10(rho),extent,ncell,domask=1.0,rhor=rhor,thetarot=thetarot)
    aphi = psicalc()[:,:,None]+rho*0
    iaphi = reinterpxy(aphi,extent,ncell,rhor=rhor,thetarot=thetarot)
    if whichvar is not None:
        if not isinstance(whichvar,basestring):
            var = whichvar()
            ivar = reinterpxy(var,extent,ncell,domask=domask,rhor=rhor,thetarot=thetarot)
        elif whichvar == "lrho":
            ivar = ilrho
        else:
            ivar = reinterpxy(avgbsqorho(),extent,ncell,domask=domask,rhor=rhor)
    #maxabsiaphi=np.max(np.abs(iaphi))
    #maxabsiaphi = 100 #50
    #ncont = 100 #30
    #levs=np.linspace(-maxabsiaphi,maxabsiaphi,ncont)
    #cset2 = plt.contour(iaphi,linewidths=0.5,colors='k', extent=extent,hold='on',origin='lower',levels=levs)
    #for c in cset2.collections:
    #    c.set_linestyle('solid')
    #CS = plt.contourf(xi,yi,zi,15,cmap=palette)
    if dostreamlines:
        Br = dxdxp[1,1]*B[1]+dxdxp[1,2]*B[2]
        Bh = dxdxp[2,1]*B[1]+dxdxp[2,2]*B[2]
        Bp = B[3]*dxdxp[3,3]
        #
        Brnorm=Br
        Bhnorm=Bh*np.abs(r)
        Bpnorm=Bp*np.abs(r*np.sin(h))
        #
        Bznorm=Brnorm*np.cos(h)-Bhnorm*np.sin(h)
        BRnorm=Brnorm*np.sin(h)+Bhnorm*np.cos(h)
        Bxnorm=BRnorm*np.cos(ph-thetarot)-Bpnorm*np.sin(ph-thetarot)
        Bynorm=BRnorm*np.sin(ph-thetarot)+Bpnorm*np.cos(ph-thetarot)
        #
        iBx = reinterpxy(Bxnorm,extent,ncell,domask=1,mirrorfactor=-1.,rhor=rhor,thetarot=thetarot)
        iBy = reinterpxy(Bynorm,extent,ncell,domask=1,mirrorfactor=-1.,rhor=rhor,thetarot=thetarot)
        if "urad" in globals():
            iibeta = reinterpxy(0.5*bsq/(gam-1)/(ug+urad),extent,ncell,domask=0,rhor=rhor,thetarot=thetarot)
        else:
            iibeta = reinterpxy(0.5*bsq/(gam-1)/(ug),extent,ncell,domask=0,rhor=rhor,thetarot=thetarot)
        ibsqorho = reinterpxy(bsq/rho,extent,ncell,domask=0,rhor=rhor,thetarot=thetarot)
        ibsqo2rho = 0.5 * ibsqorho
        xi = np.linspace(extent[0], extent[1], ncell)
        yi = np.linspace(extent[2], extent[3], ncell)
    if ax == None:
        if whichvar is not None:
            CS = plt.imshow(ivar, extent=extent, cmap = palette, norm = colors.Normalize(clip = False),origin='lower', vmin=vmin,vmax=vmax)
        else:
            CS = plt.imshow(ilrho, extent=extent, cmap = palette, norm = colors.Normalize(clip = False),origin='lower',vmin=vmin,vmax=vmax)
        plt.xlim(extent[0],extent[1])
        plt.ylim(extent[2],extent[3])
    else:
        if whichvar is not None:
            CS = ax.imshow(ivar, extent=extent, cmap = palette, norm = colors.Normalize(clip = False),origin='lower', vmin=vmin,vmax=vmax)
        else:
            CS = ax.imshow(ilrho, extent=extent, cmap = palette, norm = colors.Normalize(clip = False),origin='lower',vmin=vmin,vmax=vmax)
        if dostreamlines:
            if dovarylw==1:
                lw = 0.5+1*ftr(np.log10(amax(ibsqo2rho,1e-6+0*ibsqorho)),np.log10(1),np.log10(2))
                lw += 1*ftr(np.log10(amax(iibeta,1e-6+0*ibsqorho)),np.log10(1),np.log10(2))
                lw *= ftr(np.log10(amax(iibeta,1e-6+0*iibeta)),-3.5,-3.4)
                # if t < 1500:
                #     lw *= ftr(ilrho,-2.,-1.9)
                lw *= ftr(iaphi,0.001,0.002)
            elif dovarylw==2:
                #same as in mkframe
                lw1 = 2*ftr(np.log10(amax(ibsqo2rho,1e-6+0*ibsqorho)),np.log10(1),np.log10(2))
                lw2 = ftr(np.log10(amax(iibeta,1e-6+0*ibsqorho)),np.log10(1),np.log10(2))
                lw = 0.5 + amax(lw1,lw2)
                # lw = 0.5+1*ftr(np.log10(amax(ibsqo2rho,1e-6+0*ibsqorho)),np.log10(1),np.log10(2))
                # lw += 1*ftr(np.log10(amax(iibeta,1e-6+0*ibsqorho)),np.log10(1),np.log10(2))
                # lw *= ftr(np.log10(amax(iibeta,1e-6+0*iibeta)),-3.5,-3.4)
            elif dovarylw==3:
                lw1 = 2*ftr(np.log10(amax(ibsqo2rho,1e-6+0*ibsqorho)),np.log10(1),np.log10(2))
                lw2 = ftr(np.log10(amax(iibeta,1e-6+0*ibsqorho)),np.log10(1),np.log10(2))
                lw = 0.5 + amax(lw1,lw2)
            elif dovarylw==4:
                #new way, to avoid glitches in u_g in jet region to affect field line thickness
                lw1 = 2*ftr(np.log10(amax(ibsqorho,1e-6+0*ibsqorho)),np.log10(0.9),np.log10(1))
                val = 1. #/(6.*(gam-1.))
                lw2 = ftr(np.log10(amax(iibeta,1e-6+0*ibsqorho)),np.log10(0.9*val),np.log10(val))
                lw = 0.5 + amax(lw1,lw2)
            fstreamplot(yi,xi,iBx,iBy,density=density,downsample=downsample,linewidth=lw,detectLoops=True,dodiskfield=False,dobhfield=dobhfield,startatmidplane=False,a=a,arrowsize=arrowsize,dnarrow=dnarrow,minlengthdefault=minlengthdefault)
        ax.set_xlim(extent[0],extent[1])
        ax.set_ylim(extent[2],extent[3])
    #CS.cmap=cm.jet
    #CS.set_axis_bgcolor("#bdb76b")
    if True == cb:
        cbar = plt.colorbar(CS,ax=ax,shrink=shrink) # draw colorbar
        if fntsize is not None:
            #set font size of colorbar tick labels
            cl = plt.getp(cbar.ax, 'ymajorticklabels')
            plt.setp(cl, fontsize=fntsize)
        if label is not None: ax.set_xlabel(label,fontsize=fntsize)
    #plt.title(r'$\log_{10}\rho$ at $t = %4.0f$' % t)
    if True == pt:
        plt.title('log rho at t = %4.0f' % t)
    #if None != fname:
    #    plt.savefig( fname + '.png' )
        
def placeletter(ax1,lab,size=16,fx=0.07,fy=0.07,ha="center",va="top",color='k',bbox=None):
    ax1.text(
        ax1.get_xlim()[0]+(ax1.get_xlim()[1]-ax1.get_xlim()[0])*fx,
        ax1.get_ylim()[0]+(ax1.get_ylim()[1]-ax1.get_ylim()[0])*(1-fy), 
        r"%s" % lab, size=size,
        rotation=0., ha=ha, va=va,
        color=color,weight='regular',bbox=bbox )

def reinterp(vartointerp,extent,ncell,ncelly=None,domirrory=0,domask=1,isasymmetric=False,isasymmetricy=False,rhor=None,kval=0,domirror=True,dolimitr=True,method='linear'):
    global xi,yi,zi
    #grid3d("gdump")
    #rfd("fieldline0250.bin")
    if rhor is None:
        rhor = (1+np.sqrt(1-a**2))
    if len(vartointerp.shape) == 2:
        vartointerp = vartointerp[:,:,None]
    if kval >= vartointerp.shape[2]:
        kval = 0
    if ncelly is None:
        ncellx = ncell
        ncelly = ncell
    else:
        ncellx = ncell
        ncelly = ncelly
    maxr = 2*np.max(np.abs(np.array(extent)))
    xraw=r*np.sin(h)
    yraw=r*np.cos(h)
    x=xraw[:,:,int(kval):int(kval+1.5)].mean(2).view().reshape(-1)
    y=yraw[:,:,int(kval):int(kval+1.5)].mean(2).view().reshape(-1)
    var=vartointerp[:,:,int(kval):int(kval+1.5)].mean(2).view().reshape(-1)
    if dolimitr:
        myr=r[:,:,kval].view().reshape(-1)
        x = x[myr<maxr]
        y = y[myr<maxr]
        var = var[myr<maxr]
    #mirror
    if domirror:
        x=np.concatenate((-x,x))
        y=np.concatenate((y,y))
        kvalmirror=(kval+nz/2) % (vartointerp.shape[2])
        varmirror = np.copy(vartointerp[:,:,kvalmirror].view().reshape(-1))
        if dolimitr:
            varmirror = varmirror[myr<maxr]
        if isasymmetric==True:
            varmirror *= -1.
        var=np.concatenate((varmirror,var))
    if domirrory:
        x=np.concatenate((x,x))
        y=np.concatenate((y,-y))
        varmirror = np.copy(var)
        if isasymmetricy:
            varmirror *= -1
        var=np.concatenate((var,varmirror))
    #else do not do asymmetric part
    # define grid.
    xi = np.linspace(extent[0], extent[1], ncelly)
    yi = np.linspace(extent[2], extent[3], ncellx)
    # grid the data.
    zi = griddata((x, y), var, (xi[None,:], yi[:,None]), method=method)
    #zi[interior] = np.ma.masked
    if domask!=0:
        interior = np.sqrt((xi[None,:]**2) + (yi[:,None]**2)) < rhor*domask
        varinterpolated = ma.masked_where(interior, zi)
    else:
        varinterpolated = zi
    return(varinterpolated)

def reinterpxy(vartointerp,extent,ncell,domask=1,mirrorfactor=1,rhor=None,thetarot=0):
    global xi,yi,zi
    #grid3d("gdump")
    #rfd("fieldline0250.bin")
    if rhor is None:
        rhor = (1+np.sqrt(1-a**2))
    xraw=r*np.sin(h)*np.cos(ph-thetarot)
    yraw=r*np.sin(h)*np.sin(ph-thetarot)
    #2 cells below the midplane
    x=xraw[:,ny/2+1,:].view().reshape(-1)
    y=yraw[:,ny/2+1,:].view().reshape(-1)
    var=vartointerp[:,ny/2-1,:].view().reshape(-1)
    #mirror
    if nz*_dx3*dxdxp[3,3,0,0,0] < 0.99 * 2 * np.pi:
        x=np.concatenate((-x,x))
        y=np.concatenate((-y,y))
        var=np.concatenate((var*mirrorfactor,var))
    # define grid.
    xi = np.linspace(extent[0], extent[1], ncell)
    yi = np.linspace(extent[2], extent[3], ncell)
    # grid the data.
    zi = griddata((x, y), var, (xi[None,:], yi[:,None]), method='cubic')
    #zi[interior] = np.ma.masked
    if domask!=0:
        interior = np.sqrt((xi[None,:]**2) + (yi[:,None]**2)) < rhor*domask
        varinterpolated = ma.masked_where(interior, zi)
    else:
        varinterpolated = zi
    return(varinterpolated)
    
def ftr(x,xb,xf):
    return( amax(0.0*x,amin(1.0+0.0*x,1.0*(x-xb)/(xf-xb))) )

#############################
#
# End of movie making
#
#############################
        
    
if __name__ == "__main__":
    if len(sys.argv)>1:
        if sys.argv[1] == "mk2davg" or sys.argv[1] == "mk2dmrg":
            mkavg2dnew()
        elif sys.argv[1].startswith("mkfrm"):
            mkmov_wrapper(which=sys.argv[1])
        elif sys.argv[1].startswith("mkts"):
            mkts_wrapper()
        elif sys.argv[1].startswith("mkjethead"):
            mkjethead_wrapper()
        elif sys.argv[1].startswith("convertdump"):
            convert_wrapper(which=sys.argv[1])
        elif sys.argv[1] == "mknewmov":
            mknewmov()
        elif sys.argv[1] == "mknewmovaphi":
            mknewmov(dostreamlines=0,ncont=200,maxaphi=1000)
        elif sys.argv[1] == "mknewmovaphibondi":
            mknewmov(dostreamlines=0,ncont=200,maxaphi=200)
        elif sys.argv[1] == "mknewmovaphisane":
            mknewmov(dostreamlines=0,ncont=200,maxaphi=100,sigma=None)
#        else:
#            print( "Unknown command %s" % sys.argv[1] )



#mov_ns_accr(0, 20, 1)


def readd(n):
    if(n>=100.):
        return "dump"+str(n)
    if(n>=10.):
        return "dump0"+str(n)
    else:
        return "dump00"+str(n)

if (len(sys.argv) > 2 and sys.argv[2] == "inits"):

    #Initial state plots
    rg("gdump")
    rd("dump000")
    print("t0 = ", t)

    rd(readd(0))


    plt.clf()
    fig = plt.figure(figsize=(25, 5))

    plt.subplot(1, 3, 1)
    plt.plot(np.log10(r[:, ny//2, 0]), np.log10(bsq[:, ny//2, 0]/rho[:, ny//2, 0]))
    plt.xlabel("log $r/r_g$")
    plt.ylabel("log $(B^2/\\rho)$")


    plt.subplot(1, 3, 2)
    plt.plot(np.log10(r[:, ny//2, 0]), np.log10(rho[:, ny//2, 0]))
    plt.title("t = {:.4g} $r_g/c$".format(t), fontsize = 20)
    plt.xlabel("log $r/r_g$")
    plt.ylabel("log $\\rho$")

    plt.subplot(1, 3, 3)
    plt.plot(np.log10(r[:, ny//2, 0]), np.log(pg[:, ny//2, 0]/bsq[:, ny//2, 0]))
    plt.xlabel("log $r/r_g$")
    plt.ylabel("log $(p_g/B^2)$")

    plt.show()

#1D
#plt.clf()
#plt.loglog(r[:,ny//2,0],rho[:,ny//2,0])
#plt.xlabel("r @ equator and $\phi = 0$")
#plt.ylabel("rho")
#plt.show()

#2D plot example
#plt.clf()
#R-z plot of the logarithm of density distribution
#plt.title("density distribution")
#plc(np.log10(rho),cb=True,xy=1,xmax=100,ymax=50)
#plt.show()

#levsb = np.arange(0.01, 10.0, 0.01)
#levsrho = np.arange(-6.0, 1.0, 0.1)

#plt.clf()
#fig=plt.figure()
#plt.title("density and psicalc")
#plc(np.log10(rho),isfilled=True,cb=True,xy=1,cmap='jet',levels=levsrho)
#plc(psicalc(),isfilled=False,cb=False,xy=1,colors='black',levels=levsb)
#plt.xlim(-200,200)
#plt.ylim(-200,200)
#fig.set_size_inches(12, 10)
#plt.show()
#plt.savefig("aaa_plot.png")


def transform_up_to_polar (vec_up,isNorm=1):
    vec_up_r = dxdxp[1,1]*vec_up[1]+dxdxp[1,2]*vec_up[2]
    vec_up_h = dxdxp[2,1]*vec_up[1]+dxdxp[2,2]*vec_up[2]
    vec_up_p = vec_up[3]*dxdxp[3,3]
    #vec_up_r = dxdxp[1,1]*vec_up[1]+dxdxp[1,2]*vec_up[2]+dxdxp[1,3]*vec_up[3]
    #vec_up_h = dxdxp[2,1]*vec_up[1]+dxdxp[2,2]*vec_up[2]+dxdxp[2,3]*vec_up[3]
    #vec_up_p = vec_up[3]*dxdxp[3,3]+vec_up[2]*dxdxp[3,2]+vec_up[1]*dxdxp[3,1]
    #
    vec_rnorm=vec_up_r
    vec_hnorm=vec_up_h*np.abs(r)
    vec_pnorm=vec_up_p*np.abs(r*np.sin(h))
    #
    vec_znorm=vec_rnorm*np.cos(h)-vec_hnorm*np.sin(h)
    vec_Rnorm=vec_rnorm*np.sin(h)+vec_hnorm*np.cos(h)
    #
    vec_sq = vec_rnorm**2 + vec_hnorm**2 + vec_pnorm**2

    if (isNorm):
        return(np.array([vec_up[0],vec_rnorm,vec_hnorm,vec_pnorm]))
    else:
        return(np.array([vec_up[0],vec_up_r,vec_up_h,vec_up_p]))
