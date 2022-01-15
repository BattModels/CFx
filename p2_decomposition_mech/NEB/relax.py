from gpaw import *
from ase.io import read
from ase.optimize import BFGS
from ase.constraints import FixAtoms


# User inputs
#################################################
kpoints = [3, 4, 4]
gpoints = [64, 48, 72]
xc = 'BEEF-vdW'
#################################################

input_structure = "str.cif"
fmaxx = 0.03
maxstep = 0.1


atoms = read(input_structure)
indices = []
for i in range(len(atoms)):
    if atoms[i].symbol == 'C':
        indices.append(i)
    elif atoms[i].scaled_position[0] <= 0.25:
        indices.append(i)
atoms.set_constraint(FixAtoms(indices))
calc = GPAW(xc=xc, gpts=gpoints, kpts=kpoints, txt='bfgs.txt')
atoms.set_calculator(calc)
dyn=BFGS(atoms=atoms, trajectory='traj.traj', logfile = 'qn.log', maxstep=maxstep)
dyn.run(fmax=fmaxx)
