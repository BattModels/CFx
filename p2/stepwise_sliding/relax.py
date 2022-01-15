from gpaw import GPAW
from ase.io import read
from ase.optimize import MDMin
from ase.constraints import FixAtoms

# User inputs
#################################################
kpoints = [1, 12, 3]
gpoints = [232, 16, 64]
xc = 'BEEF-vdW'
#################################################

vacuum = 8 # Angstrom
a = 2.53959*2 # Angstrom
input_structure = "str.cif"
fmaxx = 0.03
maxstep = 0.05


atoms = read(input_structure)
atoms.set_calculator(GPAW(xc=xc, kpts = kpoints, gpts=gpoints, txt='bfgs.txt'))
cellpar = atoms.cell.cellpar()
scaled_pos = atoms.get_scaled_positions()

# Constraints- fix atoms in the first first unit cell from origin and the last two unit cells from origin
indices = []
for i in range(len(atoms)):
    if atoms[i].symbol == 'C':
        indices.append(i)

    elif scaled_pos[i][0] <= (3*a/2 + vacuum)/cellpar[0]:
        indices.append(i)

    elif scaled_pos[i][0] >= 1 - (vacuum + 5*a/2)/cellpar[0]:
        indices.append(i)

atoms.set_constraint(FixAtoms(indices))
dyn=MDMin(atoms=atoms, trajectory='traj.traj', logfile = 'qn.log')
dyn.run(fmax=fmaxx)
