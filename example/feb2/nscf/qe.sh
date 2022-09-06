#!/bin/sh
#SBATCH -J QE
#SBATCH -o out
#SBATCH -e err
#SBATCH -N 1 -n 48
#SBATCH -p normal
if [ -z "$SLURM_JOB_CPUS_PER_NODE" ];then
    SLURM_JOB_CPUS_PER_NODE=24
else
    scontrol update job $SLURM_JOB_ID JobName $(pwd|awk -F '/' '{printf("%s/%s",$(NF-1),$NF)}')
fi

mpirun -np $SLURM_JOB_CPUS_PER_NODE pw.of -npool 2 < scf.in >scf.out
mpirun -np $SLURM_JOB_CPUS_PER_NODE pw.of -npool 2 < band.in >band.out
