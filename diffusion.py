from ase.io import read
from gpaw import *
from ase.neb import NEB
from ase.optimize import BFGS
from gpaw.mpi import rank, size
from ase.parallel import parprint

#########################
n_interpol = 3
fmaxx = 0.03
kpoints = [4, 4, 3]
gpoints = [48, 48, 88]
#########################


initial = read('initial/C36F36Li.traj')
final = read('final/C36F36Li.traj')
parprint('Initial energy: ',initial.get_potential_energy())
parprint('Final energy: ',final.get_potential_energy())


images = [initial]
images += [initial.copy() for i in range(n_interpol)]
images += [final]
neb = NEB(images)
neb.interpolate('idpp')

for i,image in enumerate(images):
    calc = GPAW(xc='PBE', kpts=kpoints, gpts=gpoints, txt=str(i+1)+".txt", parallel = {'sl_auto':True})
    image.set_calculator(calc)

dyn = BFGS(neb, trajectory = 'neb.traj', logfile = 'qn.log')
dyn.run(fmax=fmaxx)












"""
n = size // n_interpol     # number of cpu's per image
j = 1 + rank // n  # my image number
assert 3 * n == size

images = [initial]

for i in range(n_interpol):
    ranks = range(i * n, (i + 1) * n)
    image = initial.copy()

    if rank in ranks:

        calc = GPAW(gpts=gpoints,
                    kpts=kpoints,
                    txt='neb{}.txt'.format(j),
                    communicator=ranks,
                    parallel = {'sl_auto':True})

        image.set_calculator(calc)
    images.append(image)

images.append(final)

neb = NEB(images, parallel=True, climb=True)
neb.interpolate()

qn = BFGS(neb, logfile='qn.log', trajectory='neb.traj')
qn.run(fmax=fmaxx)
"""
