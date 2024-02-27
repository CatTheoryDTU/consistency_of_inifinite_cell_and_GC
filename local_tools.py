import sys,os
import pickle as  pkl
from general_tools import lin_fun
from ase.io import read, write
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from scipy.optimize import curve_fit
from scripts.parse_tools import *


def plot_charge_on_Na():
    alldata=parse_meanpot(sizes,Nacharges=True)
    for thissizes in [sizes,sizes[1:]]:
        plt.plot(np.nan,np.nan,'o-',label='Mean potential',color='k')
        plt.plot(np.nan,np.nan,'d-',label='Grand canonical',markerfacecolor='None',color='k')

        for isize,size in enumerate(sizes):
            results=[[],[]]
            if size  not in alldata.keys(): continue
            if 'interface_ion'  not in alldata[size].keys(): continue
            for imet,method in enumerate(alldata[size]['interface_ion'].keys()):
                if not alldata[size]['interface_ion'][method]['FS']: continue
                dat=alldata[size]['interface_ion'][method]['FS']
                for charge in sorted(dat.keys()):
                    results[imet].append([charge,dat[charge]])
            gcresults=np.array(results[1])
            mpresults=np.array(results[0])
            if len(gcresults):
                plt.plot(gcresults[:,0],gcresults[:,1],'d-',markerfacecolor='None',color=colors[isize])
            if len(mpresults):
                plt.plot(mpresults[:,0],mpresults[:,1],'o-',color=colors[isize])
            plt.plot(np.nan,np.nan,'s',color=colors[isize],label=r'$\theta$=%1.2f'%coverages[size])
        plt.xlabel('Applied charge [e]')
        plt.ylabel('Bader charge on Na [e]')
        plt.legend()
        plt.tight_layout()
        if thissizes[0] == '2x2':
            plt.savefig('Na_charge_vs_applied_charge_all.pdf')
        else:
            plt.savefig('Na_charge_vs_applied_charge_no2x2.pdf')
        plt.close()

def plot_MP_dpot(alldata):
    for isize,size in enumerate(alldata.keys()):
        results=[]
        if 'mean_pot'  not  in alldata[size]:  continue
        mpdata=alldata[size]['mean_pot']
        if ('FS' not in mpdata.keys() or
            'slab'not in mpdata.keys()): continue
        for chFS in mpdata['FS'].keys():
            for chsl in mpdata['slab'].keys():
                if np.around(chsl-chFS,3): continue
                pot_diff = (mpdata['FS'][chsl][1]-
                            mpdata['slab'][chsl][1])
                results.append([chsl,pot_diff])
        results=np.array(results)
        if len(results) == 0:  continue
        plt.plot(results[:,0],results[:,1],'o:',color=colors[isize],linewidth=0.5)
        plt.plot(np.nan,np.nan,'s',color=colors[isize],label=r'$\theta$=%1.2f'%coverages[size])
    plt.xlabel('Applied charge [e]')
    plt.ylabel('Potential difference [V]')
    plt.legend()
    plt.tight_layout()
    plt.savefig('MP_potdiff_vs_potential.pdf')
    plt.close()
    #plt.show()

if __name__ == '__main__':
    main()
