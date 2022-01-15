from ase.optimize import BFGS
from gpaw import *
from ase.io import read
from ase.constraints import FixAtoms, FixedPlane

# User inputs
input_structure = "str.cif"
atoms = read(input_structure)
maxstep = 0.1
fmaxx = 0.01

kpoints = [6, 12, 3]
gpoints = [32, 16, 72]


index = [1,1,0]



calc = GPAW(xc='BEEF-vdW', gpts=gpoints, kpts=kpoints)  #symmetry={'point_group': False})
atoms.set_calculator(calc)

# Set constraints- Fix carbon atoms, move Li,F atoms only in plane normal to displacement direction given by index
constraints = []
constraints.append(FixAtoms(indices=[atom.index for atom in atoms if atom.symbol == 'C']))
sym = atoms.get_chemical_symbols()
for i in range(len(sym)):
    if sym[i] != 'C':
        constraints.append(FixedPlane(i,index))

atoms.set_constraint(constraints)
atoms.calc.set(txt='bfgs.txt')
dyn=BFGS(atoms=atoms, trajectory='traj.traj', logfile = 'qn.log', maxstep=maxstep)
dyn.run(fmax=fmaxx)
