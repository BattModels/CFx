from ase.io import read,write
from ase.io.trajectory import TrajectoryWriter
import os 

##################
# Tunable parameters to have: vacuum, number of moving units, number of fixed units, number of intermediate images (nstep)
# Vaccum is always to the left i.e. all atoms displaced along a by a length = vacuum for setup
# \
#  \
#   \
#  b \
#     \_______________
#             a

vacuum = 8 # Angstroms per side
num_moving_units = 4
num_fixed_units = 3        # Using converged value
nstep = 10
#a = 2.53959
##################

supercell = [num_moving_units+num_fixed_units, 1, 1]
writer= TrajectoryWriter('initial.traj','w')

for i in range(nstep):
    os.system('mkdir '+str(i))
    atoms = read('str.cif')

    atoms = atoms.repeat(supercell)
    scaled_pos = atoms.get_scaled_positions()

    for j in range(len(atoms)): # Cycle through all atoms; shift atoms that need to be shifted
        if atoms[j].symbol!='C':
            if scaled_pos[j][0] < num_moving_units/(num_moving_units+num_fixed_units):
                scaled_pos[j][1] += i/nstep
    atoms.set_scaled_positions(scaled_pos)

    # Add vacuum to the left by expanding the cell and shifting all atoms to the right
    cell = atoms.get_cell()
    cell[0][0] += 2*vacuum
    atoms.set_cell(cell, scale_atoms=False)
    pos = atoms.get_positions()
    for j in range(len(pos)):
        pos[j][0]+= vacuum
    atoms.set_positions(pos)

    writer.write(atoms)

    # Write structure to ./i/str.cif
    write(str(i)+'/str.cif',atoms)
    os.system('cp relax.py '+str(i))
    os.system('cp job.sh '+str(i))
    os.system('sed -i \'s/job_name/disl_'+str(i)+'/g\' '+str(i)+'/job.sh ')

writer.close()
