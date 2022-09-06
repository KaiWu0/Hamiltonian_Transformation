#!/bin/sh
#SBATCH -J QE
#SBATCH -o out
#SBATCH -e err
#SBATCH -N 1 -n 48
#SBATCH -p normal
#scontrol update job $SLURM_JOB_ID JobName $(pwd|awk -F '/' '{printf("%s/%s",$(NF-1),$NF)}')

mpirun -np 12 pw.x -npool 2 < scf1.in > scf1.out
mpirun -np 12 pw.x -npool 2 < scf2.in > scf2.out
#mpirun -np 2 pw.x -npool 2 < band.in > band.out
mpirun -np 1 ht.x < band.in >band.out
