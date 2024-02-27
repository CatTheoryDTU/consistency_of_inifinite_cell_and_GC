import os,sys
import  numpy as np
from ase.optimize import BFGS
from ase.visualize import view
from ase.parallel import world
from gpaw.utilities import h2gpts
from gpaw import FermiDirac
from ase.io import write, read
from ase.units import Pascal, m, Bohr

# Import solvation modules
from ase.data.vdw import vdw_radii
from gpaw.solvation import (LinearDielectric, GradientSurface,
                            SurfaceInteraction, EffectivePotentialCavity)

from gpaw.solvation.sjm import SJM, SJMPower12Potential


"""
Solvent parameters from JCP 141, 174108 (2014):
"""
u0 = 0.180                          # eV
epsinf = 78.36                      # dielectric constant of water at 298K
gamma = 18.4 * 1e-3 * Pascal * m    # Surface tension
T = 298.15                          # Temperature
vdw_radii = vdw_radii.copy()
atomic_radii = lambda atoms: [vdw_radii[n] for n in atoms.numbers]

kpts_base={'sqrt3':8,'2':6,'3':4,'4':4,'sqrt7': 6,'sqrt13':4}
grid_base={'sqrt3': 24,
               '2': 32,
               'sqrt7': 40,
               '3': 40,
               'sqrt13': 48,
               '4': 56}

home=os.getcwd()
pathinfo=home.split('/')
if pathinfo[-2] == 'GC':
    gcen=True
elif (pathinfo[-2] in ['mean_pot'] or
        pathinfo[-3] == 'CE'):
    gcen=False
else:
    print('Paths are weird!')
    sys.exit()

if pathinfo[-2] in ['GC','mean_pot']:
    size=pathinfo[-3].split('x')
    kpts=(kpts_base[size[0]],kpts_base[size[1]],1)
    grid=(grid_base[size[0]],grid_base[size[1]],120)
else:
    size=pathinfo[-4].split('x')
    rep=[int(i) for i in pathinfo[-2].split('_')[1:]]
    kpts=(kpts_base[size[0]]//rep[0],kpts_base[size[1]]//rep[1],1)
    grid=(grid_base[size[0]]*rep[0],grid_base[size[1]]*rep[1],120)

atoms=read('init.traj')

nes = [0.0,-0.2,-0.1,0.1,0.2]
if pathinfo[-3] == 'CE':
    nes=[0.0]

calc = SJM(
               doublelayer={'start': 17,'upper_limit':19.5},
               poissonsolver={'dipolelayer': 'xy'},
               gpts = grid,#h2gpts(0.18,atoms.get_cell(),idiv=8),
               kpts=kpts,
               xc='BEEF-vdW',
               potential_equilibration_mode='sim',
               write_grandcanonical_energy=gcen,
               max_pot_deviation=0.005,
               spinpol=False,
               occupations=FermiDirac(0.1),
               maxiter=400,
               cavity=EffectivePotentialCavity(
                   effective_potential=SJMPower12Potential(atomic_radii,
                                                           u0,
                                                           H2O_layer=False),
                   temperature=T,
                   surface_calculator=GradientSurface()),
               dielectric=LinearDielectric(epsinf=epsinf),
               interactions=[SurfaceInteraction(surface_tension=gamma)])

atoms.set_calculator(calc)

for ne in nes:
    #atoms.calc.ne=ne
    neglect=False
    if pathinfo[-2] == 'mean_pot':
        for prev_file in os.listdir():
            if prev_file[:2] != 'sp': continue
            if np.abs(float(prev_file.split('_')[2]) - ne) < 0.01:
                neglect=True

    elif pathinfo[-2] == 'GC':
        pot=ne*5+4.40
        for prev_file in os.listdir():
            if prev_file[:2] != 'sp': continue
            if np.abs(float(prev_file.split('_')[4].rstrip('.traj')) - pot) < 0.1:
                neglect=True
    else:
        for prev_file in os.listdir():
            if prev_file[:2] != 'sp': continue
            if np.abs(float(prev_file.split('_')[2]) - ne) < 0.01:
                neglect=True

    if neglect: continue
    if pathinfo[-2] == 'mean_pot':
        outname='q_%1.2f'%(-ne)
        calc.set(txt=outname+'.txt',
                ne=ne)
    elif pathinfo[-2] == 'GC':
        outname='pot_%1.2f'%(pot)
        calc.set(txt=outname+'.txt',
                potential=pot)
    else:
        outname='out'
        calc.set(txt=outname+'.txt')


    if np.around(ne,2) == 0:
        opt=BFGS(atoms,trajectory=outname+'.traj',logfile=outname+'.log')
        opt.run(fmax=0.05)
    else:
        atoms.get_potential_energy()


    write('sp_q_%1.3f_pot_%1.3f.traj'%(-calc.results['ne'],calc.results['electrode_potential']),atoms)

if world.rank == 0:
    os.system('touch finished')
