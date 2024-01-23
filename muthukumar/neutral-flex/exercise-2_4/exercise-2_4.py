import matplotlib
matplotlib.use('Agg')

import os

import matplotlib.pyplot as plt
import numpy as np

from tqdm import tqdm


def extract_traj(num_header, num_atoms, dump_file, save_file):
    """
    Extract atom trajectories from LAMMPS dump file.
    Store and save in 3d-array [num_atoms, 3, num_frames].

    Args:
        num_header ([int]): number of header lines in dump file
        num_atoms ([int]): number of atoms (assuming constant number)
        dump_file ([string]): LAMMPS trajectory dump file
        save_file ([string]): NumPy save file
    """
    with open(dump_file, mode='r') as file:
        lines = file.readlines()
        num_block = num_header + num_atoms
        num_frames = int(len(lines) / num_block)
        traj = np.zeros((num_atoms, 3, num_frames))

        for frame in range(num_frames):
            # extract atom positions block for frame
            block = lines[frame * num_block:(frame+1) * num_block][-num_atoms:]

            # convert each line from string to float, store as array
            frame_pos = np.array(
                [np.array(x.split()) for x in block],
            )[:, [0, 3, 4, 5]].astype(np.float32)

            # sort frame by atom #, slice into traj array
            traj[:, :, frame] = frame_pos[frame_pos[:, 0].argsort()][:, -3:]

    file.close()
    
    np.save(save_file, traj)


def calculate_cm(traj):
    """
    Calculate center of mass from atom trajectories data. 
    Assuming data is in [num_atoms, 3, num_frames] shape.

    Args:
        traj ([numpy.ndarray]): trajectory array

    Returns:
        center of mass trajectory
    """
    return np.mean(traj, axis=0)


def calculate_msd(traj, frame_step, prod_start, idx=None):
    """
    Calculate the Mean Square Displacement (MSD) of a single particle.

    MSD_i(tau) = <[R_i(tau + t) - R_i(t)]^2>

    Args:
        traj ([numpy.ndarray]): trajectory array
        frame_step ([int]): time step between each frame
        prod_start ([int]): time step when production stage starts
        idx ([int]): the label of the particle of interest [default: None]

    Returns:
        msd_i ([numpy.ndarray]): MSD of particle idx
        taus ([numpy.ndarray]): time delays
        idx ([int]): the label of the particle of interest
    """
    if idx is None:
        idx = np.random.randint(low=1, high=traj.shape[0] + 1, dtype=int)

    # remove duplicate frames due to restart
    traj = check_restart(traj)

    # extract monomer trajectory
    if traj.ndim == 2:
        traj_i = traj  # assuming we're calculating center of mass

    if traj.ndim == 3:
        traj_i = traj[idx-1, :, :]

    # remove equilibration stage
    prod_idx = int(prod_start / frame_step)
    traj_i_prod = traj_i[:, prod_idx:]

    # calculate position vector at each frame
    r_i = np.linalg.norm(traj_i_prod, axis=0)

    t = len(r_i)
    msd_i = np.zeros(t - 1)
    taus = np.arange(1, t)

    # calculate MSD
    for i, tau in enumerate(taus):
        diff = t - tau
        msd_i[i] = np.square(r_i[-diff:] - r_i[:diff]).mean()

    return msd_i, taus, idx


def calculate_msd_V2(traj, frame_step, prod_start, idx=None):
    """
    Calculate the Mean Square Displacement (MSD) of a single particle.

    MSD_i(tau) = <[R_i(tau) - R_i(0)]^2>

    Args:
        traj ([numpy.ndarray]): trajectory array
        frame_step ([int]): time step between each frame
        prod_start ([int]): time step when production stage starts
        idx ([int]): the label of the particle of interest [default: None]

    Returns:
        msd_i ([numpy.ndarray]): MSD of particle idx
        taus ([numpy.ndarray]): time delays
        idx ([int]): the label of the particle of interest
    """
    if idx is None:
        idx = np.random.randint(low=1, high=traj.shape[0] + 1, dtype=int)
    
    # remove duplicate frames due to restarts
    traj = check_restart(traj)

    # extract monomer trajectory
    if traj.ndim == 2:
        traj_i = traj  # assuming we're calculating center of mass
        
    if traj.ndim == 3:
        traj_i = traj[idx-1, :, :]

    # remove equilibration stage
    prod_idx = int(prod_start / frame_step)
    traj_i_prod = traj_i[:, prod_idx:]

    # calculate position vector at each frame
    r_i = np.linalg.norm(traj_i_prod, axis=0)
    taus = np.arange(0, len(r_i))

    # calculate MSD
    msd_i = np.square(r_i - r_i[0])

    return msd_i, taus, idx


def check_restart(traj):
    """
    Remove repeating frames in a trajectory array that are a result of
    starting a simulation from a restart file.

    Args:
        traj ([numpy.ndarray]): trajectory array

    Returns:
        traj ([numpy.ndarray]): trajectory array with repeats removed
    """
    repeats = []

    if traj.ndim == 2:
        for i in range(traj.shape[-1] - 1):
            if np.array_equal(traj[:, i], traj[:, i+1]):
                repeats.append(i)
    
    if traj.ndim == 3:
        for i in range(traj.shape[-1] - 1):
            if np.array_equal(traj[:, :, i], traj[:, :, i+1]):
                repeats.append(i)
    
    traj = np.delete(traj, repeats, axis=-1)
    
    return traj 


def plot_Rg_evol():
    """Plot time evolution of radius of gyration."""
    filenames = ['{}/radius_of_gyration.dat'.format(dir) for dir in os.listdir() 
      if dir.startswith('RUN_')]
    data = np.array([np.loadtxt(fname, skiprows=1) for fname in filenames]).mean(axis=0)

    dt = 1e-3
    Rg = np.sqrt(data[:, 1])
    t = data[:, 0] * dt

    plt.figure()
    for fname in filenames:
        data = np.loadtxt(fname, skiprows=1)
        plt.plot(data[:, 0] * dt, np.sqrt(data[:, 1]), 'b', alpha=0.01, linewidth=0.75)

    plt.plot(t , Rg, 'b', linewidth=0.75)
    plt.title('N = 120 ({} independent runs)'.format(len(filenames)))
    plt.xlabel('t')
    plt.ylabel('$R_g$')
    plt.tight_layout()
    plt.savefig('figures/Rg_vs_t.png', dpi=600)


# Exercise 2
#######################################################################
#for i in tqdm([dir for dir in os.listdir() if dir.startswith('RUN_')]):
#    extract_traj(
#        num_header=9,
#        num_atoms=120,
#        dump_file='{}/dump.lammpstrj'.format(i),
#        save_file='{}/atom_traj.npy'.format(i),
#    )
#
#plot_Rg_evol()
#
# Part 4
########
frame_step = 5000
prod_start = 3000000
filenames = ['{}/atom_traj.npy'.format(dir) for dir in os.listdir() if dir.startswith('RUN_')]

# <[R_i(tau + t) - R_i(t)]^2>
print('<[R_i(tau + t) - R_i(t)]^2>')
for i in tqdm([20, 40, 60, 80, 100, 'CM']):
    if type(i) == int:
        msd_i, TAUS, _ = np.array(
            [calculate_msd(np.load(fname), frame_step, prod_start, idx=i) \
            for fname in filenames],
        ).mean(axis=0)
        np.savetxt('msd-files/msd_{}.txt'.format(i), np.stack((TAUS, msd_i)).T, header='tau MSD')

    if type(i) == str:
        msdCM, TAUS, _ = np.array(
            [calculate_msd(calculate_cm(np.load(fname)), frame_step, prod_start, idx=0) \
            for fname in filenames],
        ).mean(axis=0)
        np.savetxt('msd-files/msd_CM.txt', np.stack((TAUS, msdCM)).T, header='tau MSD')

# <[R_i(tau) - R_i(0)]^2>
print('<[R_i(tau) - R_i(0)]^2>')
for i in tqdm([20, 40, 60, 80, 100, 'CM']):
    if type(i) == int:    
        msd_i, TAUS, _ = np.array(
            [calculate_msd_V2(np.load(fname), frame_step, prod_start, idx=i) \
            for fname in filenames],
        ).mean(axis=0)
        np.savetxt(
            'msd-files/msd_{}_V2.txt'.format(i), np.stack((TAUS, msd_i)).T, 
            header='tau MSD',
        )

    if type(i) == str:
        msdCM, TAUS, _ = np.array(
            [calculate_msd_V2(calculate_cm(np.load(fname)), frame_step, prod_start, idx=0) \
            for fname in filenames],
        ).mean(axis=0)
        np.savetxt('msd-files/msd_CM_V2.txt', np.stack((TAUS, msdCM)).T, header='tau MSD')

