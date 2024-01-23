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


def extract_traj_V2(num_header, num_atoms, dump_file, save_dir):
    """
    Extract atom trajectories from LAMMPS dump file.
    Store and save in 3d-array [num_atoms, 3, num_frames].
    *Store additional atom type array.*

    Args:
        num_header ([int]): number of header lines in dump file
        num_atoms ([int]): number of atoms (assuming constant number)
        dump_file ([string]): LAMMPS trajectory dump file
        save_dir ([string]): directory to save files to
    """
    with open(dump_file, mode='r') as file:
        lines = file.readlines()
        num_block = num_header + num_atoms
        num_frames = int(len(lines) / num_block)
        traj = np.zeros((num_atoms, 3, num_frames))
        types = np.zeros((num_atoms, 2))

        for frame in range(num_frames):
            # extract atom positions block for frame
            block = lines[frame * num_block:(frame+1) * num_block][-num_atoms:]

            # convert each line from string to float, store as array
            frame_pos = np.array(
                [np.array(x.split()) for x in block],
            )[:, [0, 3, 4, 5]].astype(np.float32)
            
            # store atom type information using first frame
            if frame == 0:
                types = np.array(
                    [np.array(x.split()) for x in block],
                )[:, [0, 2]].astype(int)

            # sort frame by atom #, slice into traj array
            traj[:, :, frame] = frame_pos[frame_pos[:, 0].argsort()][:, -3:]
            types = types[types[:, 0].argsort()]

    file.close()
    
    np.save('{}/atom_traj.npy'.format(save_dir), traj)
    np.savetxt('{}/atom_types.txt'.format(save_dir), types, fmt='%i',
        header='atom type')


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


def plot_Rg_evol(path, dt_eq=1e-5, dt_pr=1e-3, prod_start=None):
    """Plot time evolution of radius of gyration.

    Args:
        path [(string)]: path to data, either a directory or file
        dt_eq [(float)]: timestep size for equilibration [default: 1e-5]
        dt_pr [(float)]: timestep size for production [default: 1e-3]
        prod_start [(int)]: timestep where production run begins [default: None]
    """
    path = os.path.abspath(path)

    if path.endswith('.dat'):
        data = np.loadtxt(path, skiprows=1)
        box_str = int(path.split('/')[-4].split('-')[-1])
        sigma_str = float(path.split('/')[-3].split('-')[-1].replace('_', '.'))  # UGLY WOW

        t = data[:, 0]
        Rg = data[:, 1]

    else:
        if path.split('/')[-1].startswith('sigma'):
            box_str = int(path.split('/')[-2].split('-')[-1])
            sigma_str = float(path.split('/')[-1].split('-')[-1].replace('_', '.'))
            # plot single curve
            filenames = [f'{dir}/radius_of_gyration.dat' for
              dir in os.listdir(path) if dir.startswith('RUN_')]
            data = np.array(
                [np.loadtxt(fname, skiprows=1) for fname in filenames],
            ).mean(axis=0)
            t = data[:, 0]
            Rg = data[:, 1]

        if path.split('/')[-1].startswith('box-1'):
            # plot curve for each charge separation distance (sigma)
            sigmas = ['0.25 nm', '0.7 nm', '1.0 nm']
            data = np.array(
                [[np.loadtxt(f'{pdir}/{cdir}/radius_of_gyration.dat',
              skiprows=1) for cdir in os.listdir(os.path.join(path, pdir)) if 
              cdir.startswith('RUN_')] for pdir in os.listdir(path) if 
              pdir.startswith('sigma')
            ]).mean(axis=1)
            t = data[0, :, 0]
            Rg = data[:, :, 1]

        else:  # assuming only other option is "RUN_*" directories
            box_str = int(path.split('/')[-2].split('-')[-1])
            sigma_str = float(path.split('/')[-1].split('-')[-1].replace('_', '.'))
            filenames = [f'{dir}/radius_of_gyration.dat' for dir
              in os.listdir(path) if dir.startswith('RUN_')]
            data = np.array(
                [np.loadtxt(fname, skiprows=1) for fname in filenames],
            ).mean(axis=0)
            t = data[:, 0]
            Rg = data[:, 1]

    if prod_start is not None:
        t = np.hstack((t[t <= prod_start] * dt_eq, t[t > prod_start] * dt_pr))
        
        fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
        if Rg.ndim == 2:  # if plotting multiple curves (for multiple charge distances)
            for i, x in enumerate(Rg):
                ax1.plot(t, x, linewidth=0.75, label=sigmas[i])
                ax2.plot(t, x, linewidth=0.75, label=sigmas[i])
        
        if Rg.ndim == 1:
            ax1.plot(t, Rg, linewidth=0.75)
            ax2.plot(t, Rg, linewidth=0.75)

        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
        plt.grid(False)
        try:
            fig.suptitle(f'Polyelectrolyte [N = 60, $\sigma$ = {sigma_str} nm, box = {box_str} LJ units]')
        except:
            pass
        
        plt.xlabel('$t$ [LJ units]')
        ax1.set_ylabel('$R_g$')
        ax1.set_xlim(0, prod_start * dt_eq)
        ax2.set_xlim(prod_start * dt_pr, t[-1])

        ax1.spines.right.set_visible(False)
        ax2.spines.left.set_visible(False)
        ax1.yaxis.tick_left()
        ax1.tick_params(labelright=False)
        ax2.yaxis.tick_right()

        d = 0.015
        kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
        ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)
        ax1.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

        kwargs.update(transform=ax2.transAxes)
        ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)
        ax2.plot((-d, +d), (-d, +d), **kwargs)

    else:
        t = t * dt_eq
        sigma = path.split('/')[-1].split('-')[-1].replace('_', '.')

        fig, ax = plt.subplots(1, 1)
        ax.plot(t, Rg, 'b', linewidth=0.75)
        try:
            ax.set_title(f'Polyelectrolyte [N = 60, $\sigma$ = {sigma_str} nm, box = {box_str} LJ units]')
        except:
            pass
        
        ax.set_xlabel('$t$ [LJ units]')
        ax.set_ylabel('$R_g$')
        
    #ax2.legend()
    fig.tight_layout()
    fig.savefig('figures/Rg_vs_t.png', dpi=600)
   

def plot_MSD(path, dt=1e-3, frame_step=5000):
    """Plot time evolution of mean square displacement.
    <[R_i(t + tau) - R_i(t)]^2>

    Args:
        path [(string)]: path to data, either a directory or file
        dt [(scalar)]: timestep size for production [default: 1e-3]
        frame_step([int]): steps between each frame [default: 5000]
    """
    path = os.path.abspath(path)

    if path.endswith('.txt'):
        data = np.loadtxt(path, skiprows=1)

        t = data[:, 0]
        msd = data[:, 1]

    else:
        filenames = [f'msd-files/{file}' for file in 
          os.listdir(os.path.join(path, 'msd-files')) if '_V3' not in file]
        data = np.array(
            [np.loadtxt(fname, skiprows=1) for fname in filenames],
        )

    idx = [20, 60, 100, 'CM']
    fig, ax = plt.subplots()

    for i, a in enumerate(data):
        ax.plot(a[:, 0] * frame_step * dt, a[:, 1], linewidth=0.75, label=f'i = {idx[i]}')

    ax.set_xscale('log')
    ax.set_yscale('log') 
    ax.set_xlabel('$\\tau$')
    ax.set_ylabel(r'$ \langle [R_{i}(t + \tau) - R_{i}(t)]^2 \rangle $')
    ax.set_title('N = 120')
    ax.legend()
    
    fig.tight_layout()
    fig.savefig('figures/msd_vs_tau.png', dpi=600)

def plot_MSD_V2(path, dt=1e-3, frame_step=5000):
    """Plot time evolution of mean square displacement.
    <[R_i(tau) - R_i(0)]^2>

    Args:
        path [(string)]: path to data, either a directory or file
        dt [(scalar)]: timestep size for production [default: 1e-3]
        frame_step([int]): steps between each frame [default: 5000]
    """
    path = os.path.abspath(path)

    if path.endswith('.txt'):
        data = np.loadtxt(path, skiprows=1)

        t = data[:, 0]
        msd = data[:, 1]

    else:
        filenames = [f'msd-files/{file}' for file in 
          os.listdir(os.path.join(path, 'msd-files')) if '_V2' in file]
        data = np.array(
            [np.loadtxt(fname, skiprows=1) for fname in filenames],
        )

    idx = [20, 60, 100, 'CM']
    fig, ax = plt.subplots()

    for i, a in enumerate(data):
        ax.plot(a[:, 0] * frame_step * dt, a[:, 1], linewidth=0.75, label=f'i = {idx[i]}')

    ax.set_xscale('log')
    ax.set_yscale('log') 
    ax.set_xlabel('$\\tau$')
    ax.set_ylabel(r'$ \langle [R_{i}(\tau) - R_{i}(0)]^2 \rangle $')
    ax.set_title('N = 120')
    ax.legend()
    
    fig.tight_layout()
    fig.savefig('figures/msd_vs_tau_V2.png', dpi=600)


def plot_gamma_Rg(path):
    """
    Args:
        path ([string])
    
    Returns:

    """
    path = os.path.abspath(path)

    data = np.array(
        [[np.loadtxt('{}/{}/radius_of_gyration.dat'.format(pdir, cdir),
      skiprows=1) for cdir in os.listdir(os.path.join(path, pdir)) if 
      cdir.startswith('RUN_')] for pdir in os.listdir(path) if 
      pdir.startswith('sigma')
    ]).mean(axis=1)
    t = data[0, :, 0]
    Rg = data[:, -500:, 1].mean(axis=1)
    gammas = 1 / np.array([0.25, 0.7, 1.0])

    fig, ax = plt.subplots()
    ax.scatter(gammas, Rg)
    ax.set_xlabel('t')
    ax.set_ylabel('$R_g$')
     
    fig.tight_layout
    fig.savefig('figures/gamma_vs_Rg.png', dpi=600)


def calculate_adsorbed_ions(traj, types, chain_type, ion_type, cutoff):
    """
    Calculate the fraction of adsorbed ions on a polymer chain backbone.
    
    Args:
        traj [(numpy.ndarray)]: trajectory array 
        types [(numpy.ndarray)]: atom types array
        chain_type [(int)]: atom type for chain
        ion_type [(int)]: atom type for ions
        cutoff [(scalar)]: cutoff criterion for adsorption

    Returns:
        frac_ads [(numpy.ndarray)]: time evolution of fraction of adsorbed ions
    """
    # commenting out for now, i will need to update this func for handling multiple types
    #if type(chain_type) is not list:
    #    chain_type = [chain_type]

    #if type(ion_type) is not list:
    #    ion_type = [ion_type]

    window = 100
    chain_traj = traj[types[:, 1] == chain_type, :, -window:]
    ions_traj = traj[types[:, 1] == ion_type, :, -window:]
    print(chain_traj.shape, ions_traj.shape)
    #frac_ads = np.zeros(window)
    #for i in range(window):
    #    count = 0
    #    for ion in ions_traj[:, :, i]:
    #        dist = np.abs(np.linalg.norm(ion - chain_traj[:, :, i]).min())
    #        if i == 99:
    #            print(dist)
    #        if dist <= cutoff:
    #            count += 1

    #    frac_ads[i] = count / ions_traj.shape[0]
    #
    #return frac_ads

