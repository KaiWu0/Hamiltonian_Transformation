#!/bin/sh
#SBATCH -J QE
#SBATCH -o out
#SBATCH -e err
#SBATCH -N 1 -n 48
#SBATCH -p short
scontrol update job $SLURM_JOB_ID JobName $(pwd|awk -F '/' '{printf("%s/%s",$(NF-1),$NF)}')

mpirun -np $SLURM_JOB_CPUS_PER_NODE pw.x -npool 4 < scf.in >scf.out
mpirun -np $SLURM_JOB_CPUS_PER_NODE pw.x -npool 4 < band.in >band.out
#mpirun -np 1 ht.x < band.in >band.out
