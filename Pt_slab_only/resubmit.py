#!/usr/bin/python
import numpy as np
from structure_setup_tools import *
import shutil

scriptdir='/home/cat/geokast/with_Surdarshan/Cations_on_Pt_method_comparison/relaxed/scripts'

""" Create cation structures """
import os,sys
from copy import deepcopy
from ase.io import write,read

home=os.getcwd()
if __name__ == '__main__':
    # charges you want to run at
    ads = home.split('/')[-1].split('_')[0]
    charges = [-0.2, -0.1, 0.0, 0.1, 0.2]
    # cell sizes
    cell_sizes = [['sqrt3','sqrt3'], ['2','2'], ['sqrt7','sqrt7'], ['3','3'], ['sqrt13', 'sqrt13'], ['4','4']]

    for cell in cell_sizes:
        homedir=cell[0]+'x'+cell[1]
        if homedir not in os.listdir(home):
            os.mkdir(home+'/'+homedir)

        for scheme in ['GC','mean_pot']:
            schemedir=home+'/'+homedir+'/'+scheme
            if scheme not in os.listdir(homedir):
                os.mkdir(schemedir)
            for state in ['state_slab']:
                statedir=schemedir+'/'+state
                if state not in os.listdir(schemedir):
                    print(statedir+' doesnt exist')
                    continue
                trajfiles=[i for i in os.listdir(statedir)
                        if (i.split('.')[-1] == 'traj' and
                            i[:2] == 'sp')]

                if 'finished' in os.listdir(statedir):
                    print(statedir+' is finished')
                    continue

                if scheme == 'GC':
                    relaxfile='pot_4.40.traj'
                else:
                    relaxfile='q_-0.00.traj'

                if relaxfile in os.listdir(statedir):
                    if os.path.getsize(statedir+'/'+relaxfile):
                                print('inputfile has been updated for '+statedir)
                                shutil.copy(statedir+'/'+relaxfile,statedir+'/'+relaxfile+'.bakk')
                                shutil.move(statedir+'/init.traj',statedir+'/init_orig.traj')
                                shutil.move(statedir+'/'+relaxfile,statedir+'/init.traj')
                shutil.copy(scriptdir+'/run.py',statedir)

        continue
        for scheme in ['CE']:
            schemedir=home+'/'+homedir+'/'+scheme
            if scheme not in os.listdir(homedir):
                os.mkdir(schemedir)
            repetitions=[]
            for multix in [1,2]:
                for multiy in [1,2]:
                    if multiy*multix in repetitions: continue
                    repetitions.append(multix*multiy)

                    multdir=schemedir+'/rep_%s_%s'%(multix,multiy)
                    if 'rep_%s_%s'%(multix,multiy) not in os.listdir(schemedir):
                        print(multdir+'does not exist')

                    for state in ['state_FS','state_slab']:

                        statedir=multdir+'/'+state
                        if state not in os.listdir(multdir):
                            os.mkdir(statedir)
                        relaxfile='out.traj'
                        if relaxfile in os.listdir(statedir):
                            if os.path.getsize(statedir+'/'+relaxfile):
                                        shutil.copy(statedir+'/'+relaxfile,statedir+'/'+relaxfile+'.bakk')
                                        shutil.move(statedir+'/init.traj',statedir+'/init_orig.traj')
                                        shutil.move(statedir+'/'+relaxfile,statedir+'/init.traj')
                                        print(statedir+' has been updated')
                        shutil.copy(scriptdir+'/run.py',statedir)
os.chdir(home)
#os.system('bash %s/submit.sh'%scriptdir)

