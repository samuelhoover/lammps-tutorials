import os
import re

import matplotlib.pyplot as plt
import numpy as np

from tqdm import trange

from neutral-flex import exercise


def extract_trajectories(lammps_file, data_file, dump_file, save_file):
    """
    Extract atom trajectories from LAMMPS dump file and store in 3d-array [frames, atoms, 3].

    Args:
        lammps_file ([string]): LAMMPS input file
        data_file ([string]): LAMMPS data file
        dump_file ([string]): LAMMPS trajectory dump file
        save_file ([string]): NumPy save file
    """

    # Get number of atoms
    for line in open(data_file, mode='r'):
        if re.search('atoms$', line):
            num_atoms = int(line.split()[0])
            break
    
    # Get number of frames
    time_steps = 0
    for line in open(lammps_file, mode='r'):
        if re.search(dump_file.split('/')[-1], line):
            del_time = int(line.split()[4])

        if re.search('^run  ', line):
            time_steps = int(line.split()[-1])

    num_frames = int(time_steps / del_time) + 1  # include first frame

    # Extract atom trajectories from dump file
    dump_array = np.zeros((num_frames, num_atoms, 3))
    atom_traj = np.zeros((num_frames, 3))
    for atom in trange(num_atoms):
        t_step = 0
        for line in open(dump_file, mode='r'):
            if re.search('^{} '.format(atom+1), line):
                x, y, z = line.split()[2:5]
                atom_traj[t_step, :] = float(x), float(y), float(z)
                t_step += 1
            
        dump_array[:, atom, :] = atom_traj

    np.save(save_file, dump_array)


def calculate_cm(trajectories):
    """
    Calculate center of mass from atom trajectories data. 
    Assuming data is in num_time_steps x num_atoms x 3 shape.

    Args:
        trajectories ([numpy.ndarray]): trajectory file

    Returns:
        cm: center of mass trajectory
    """

    return np.mean(trajectories, axis=1)


cwd = os.getcwd()


def main():
    if cwd.split('/')[-2].startswith('neutral'):
        if cwd.split('/')[-1].startswith('exercise-1'):
            exercise.exercise-1_1()
            
        if cwd.split('/')[-1].startswith('exercise-2'):
            if cwd.split('/')[-1].endswith('1'):

            if cwd.split('/')[-1].endswith('2'):

            if cwd.split('/')[-1].endswith('3'):

            if cwd.split('/')[-1].endswith('4'):

        if cwd.split('/')[-1].startswith('exercise-3'):
            if cwd.split('/')[-1].endswith('1'):

            if cwd.split('/')[-1].endswith('2'):

            if cwd.split('/')[-1].endswith('3'):

            if cwd.split('/')[-1].endswith('4'):

    if cwd.split('/')[-2].startswith('charged'):


if __name__ == '__main__':
    main()