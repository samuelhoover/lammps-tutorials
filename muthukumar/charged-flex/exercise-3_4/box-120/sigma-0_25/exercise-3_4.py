import matplotlib
matplotlib.use('Agg')

import operator
import os

import matplotlib.pyplot as plt
import numpy as np

from functools import reduce
from tqdm import tqdm

from utils.utils import extract_traj_V2, check_restart, \
                        plot_Rg_evol, calculate_adsorbed_ions, \
                        plot_gamma_Rg


# Exercise 3
#######################################################################
path = os.path.abspath('.')
#dir_list = [
#  [os.path.join(pdir, cdir) for cdir in os.listdir(os.path.join(path, pdir)) \
#  if cdir.startswith('RUN_')] for pdir in os.listdir(path) if \
#  pdir.startswith('sigma')
#]
#dirs = reduce(operator.add, dir_list)
#for i in tqdm(dirs):
#    extract_traj_V2(
#        num_header=9,
#        num_atoms=120,
#        dump_file=f'{i}/dump.lammpstrj',
#        save_dir=f'{i}',
#    )

plot_Rg_evol(path='.', dt_pr=5e-3, prod_start=6e6)

#for d in tqdm([x for x in os.listdir(path) if x.startswith('RUN_')]):
#    extract_traj_V2(
#        num_header=9,
#        num_atoms=120,
#        dump_file=f'{d}/dump.lammpstrj',
#        save_dir=f'{d}',
#    )
#
#rundir = 'RUN_01'
#traj = np.load(f'{rundir}/atom_traj.npy')
#types = np.loadtxt(f'{rundir}/atom_types.txt', skiprows=1)
#cutoff = 1.12
#frac = calculate_adsorbed_ions(traj, types, 1, 2, cutoff)
##print(frac, frac.var())

#plot_gamma_Rg(path='.')

