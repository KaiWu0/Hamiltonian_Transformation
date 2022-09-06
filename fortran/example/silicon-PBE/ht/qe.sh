#!/bin/sh
#SBATCH -J QE
#SBATCH -o out
#SBATCH -e err
#SBATCH -N 1 -n 48
#SBATCH -p short
scontrol update job $SLURM_JOB_ID JobName $(pwd|awk -F '/' '{printf("%s/%s",$(NF-1),$NF)}')

mpirun -np $SLURM_JOB_CPUS_PER_NODE pw.o -npool 4 < scf.in >scf.out
ht.of < band.in >band.out
