#!/bin/bash

DATAPATH=/home/samhoover/projects/lammps-tutorial/tutorial-1/neutral-flex/exercise-2_4-NO-RECENTER

rsync -avzP samhoover@172.30.64.26:$DATAPATH/msd-files/msd*.txt ./msd_files/
rsync -avzP samhoover@172.30.64.26:$DATAPATH/figures/*.png ./figures/
