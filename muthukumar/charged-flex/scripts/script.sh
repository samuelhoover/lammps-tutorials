#!/bin/sh
#!/bin/bash
#SBATCH -p long
#SBATCH -J SinM75-0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=14-00
#SBATCH --output=energy.out
#SBATCH --cpus-per-task=8
export OMP_NUM_THREADS=1

mpirun -np 8 ./lmp -in in.complexation 


