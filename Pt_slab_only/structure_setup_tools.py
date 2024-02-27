import numpy as np
from ase.constraints import FixAtoms
from ase import Atoms
from ase.build import fcc111_root, add_adsorbate, fcc111
def create_slab(cell):
        print('Creating cell: {l}'.format(l=cell))
        if 'sqrt' in cell[0]:
            # Create root surfaces
            slab = fcc111_root('Pt', root=int(cell[0].replace('sqrt','')), a=4.00, size=[1, 1, 3], vacuum=10)
        else:
            slab = fcc111('Pt', a=4.00, size=[int(a) for a in cell]+[3], vacuum=10)

        slab.translate([0,0,-3])
        slab.cell[2,2]= 21.
        slabz=np.unique(slab.positions[:,2])[:2]
        c = FixAtoms([atom.index for atom in slab if atom.position[2] in slabz])
        slab.set_constraint(c)

        return slab

def add_cation(FS,ads,cell):
        cation = Atoms(ads)
        if ads == 'H':
            ads_dist=1.5
        else:
            ads_dist=2.

        if cell[0] == 'sqrt13':
            cation.positions[0] = [5.034,1.696,11.619+ads_dist]
            FS+=cation
        elif cell[0] == 'sqrt7':
            cation.positions[0] = [3.653,1.697,11.619+ads_dist]
            FS+=cation
        elif cell[0] == 'sqrt3':
            cation.positions[0] = [4.491,3.536,11.619+ads_dist]
            FS+=cation
        else:
            add_adsorbate(FS, cation,  ads_dist, 'fcc')
