from ase.io import *
from ase.io.trajectory import TrajectoryWriter
import os

filename = "final.cif"
index = [1,1,0]
nstep = 10

path = str(index[0])+str(index[1])+str(index[2])
command = "mkdir "+path
os.system(command)

writer = TrajectoryWriter(path+'/initial.traj','a')
latprm = read(filename).get_cell_lengths_and_angles()

for i in range(nstep+1):
    atoms = read(filename)
    scaled_pos = atoms.get_scaled_positions()
    sym = atoms.get_chemical_symbols()
    for j in range(len(sym)):
        if sym[j]!='C':
            scaled_pos[j][0] = scaled_pos[j][0]+i*index[0]/nstep
            scaled_pos[j][1] = scaled_pos[j][1]+i*index[1]/nstep
            scaled_pos[j][2] = scaled_pos[j][2]+i*index[2]/nstep
    atoms.set_scaled_positions(scaled_pos)
    writer.write(atoms)
    command = "mkdir "+path+"/"+str(i)
    os.system(command)
    write(path+'/'+str(i)+'/str.cif',atoms)

    command = "cp job.sh relax.py "+path+"/"+str(i)+"/."
    os.system(command)

    # Edit job.sh job name
    newjobname = path+'_'+str(i)
    command = 'sed -i \'s/jobname/'+newjobname+'/g\' ./'+path+'/'+str(i)+'/job.sh '
    os.system(command)

writer.close()
