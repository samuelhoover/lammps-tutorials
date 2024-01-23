#!/bin/sh
#!/bin/bash
#SBATCH -p short
#SBATCH -J SinM75-0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=1-00
#SBATCH --output=energy.out
#SBATCH --cpus-per-task=8
export OMP_NUM_THREADS=1

mpirun -np 1 ./lmp -in in.complexation 


