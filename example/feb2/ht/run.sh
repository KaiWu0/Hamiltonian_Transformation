#!/bin/sh
#SBATCH -J QE
#SBATCH -o out
#SBATCH -e err
#SBATCH -N 1 -n 48
#SBATCH -p normal
#scontrol update job $SLURM_JOB_ID JobName $(pwd|awk -F '/' '{printf("%s/%s",$(NF-1),$NF)}')

mpirun -np 24 pw.of -npool 4 < scf.in >scf.out
#mpirun -np 24 pw.x -npool 2 < band.in >band.out
ht.o < band.in >band.out
