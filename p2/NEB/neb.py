from ase.io import read
from gpaw import *
from ase.neb import NEB
from ase.optimize import BFGS
from ase.parallel import parprint
import glob


#########################
n_interpol = 3
fmaxx = 0.03
kpoints = [3, 4, 4]
gpoints = [64, 48, 72]
xc = 'BEEF-vdW'
#########################


initial = read(glob.glob('initial/*.traj')[0])
final = read(glob.glob('final/*.traj')[0])
parprint('Initial energy: ',initial.get_potential_energy())
parprint('Final energy: ',final.get_potential_energy())
parprint('Initial - Final energy: ', initial.get_potential_energy()-final.get_potential_energy())


images = [initial]
images += [initial.copy() for i in range(n_interpol)]
images += [final]
neb = NEB(images)
neb.interpolate('IDPP')


for i,image in enumerate(images):
    calc = GPAW(xc=xc, kpts=kpoints, gpts=gpoints, txt=str(i+1)+".txt", parallel = {'sl_auto':True})
    image.set_calculator(calc)


dyn = BFGS(neb, trajectory = 'neb.traj', logfile = 'qn.log')
dyn.run(fmax=fmaxx)
