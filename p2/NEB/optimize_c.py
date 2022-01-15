from ase.io import read, write
from gpaw import *
from ase.eos import EquationOfState
import numpy as np

import matplotlib
matplotlib.use('Agg')

# User inputs 
############################################################ 
kpoints = [3, 4, 4]
gpoints = [64, 48, 64]
xc = 'BEEF-vdW'
############################################################ 

v = []
e = []

def genstr(scale):
    atoms = read('str.cif')

    # Scale the cell
    cell = atoms.get_cell()
    cell[2][2] = scale * cell[2][2]
    atoms.set_cell(cell, scale_atoms=False)

    # Translate Li, F atoms vertically
    for j in range(len(atoms)):
        if atoms[j].symbol != 'C':
            atoms[j].position[2] += (scale-1) * 11.08/2
    return atoms


# Main
for i, scale in enumerate(np.linspace(0.85, 1.0, num=5)):
    atoms = genstr(scale)
    atoms.calc = GPAW(xc=xc, kpts=kpoints, gpts=gpoints, txt=str(i)+'.txt')
    e.append(atoms.get_potential_energy())
    v.append(atoms.get_volume())

eos = EquationOfState(v, e, eos='birchmurnaghan')
v0, e0, B = eos.fit()
eos.plot('eos.png', show=False)
opt_scale = v0/read('str.cif').get_volume()

write('opt_str.cif', genstr(opt_scale))
