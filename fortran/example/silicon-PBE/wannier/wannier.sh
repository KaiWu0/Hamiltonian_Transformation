#!/bin/bash
#SBATCH -J QE
#SBATCH -o out
#SBATCH -e err
#SBATCH -N 1 -n 48
#SBATCH -p normal
#scontrol update job $SLURM_JOB_ID JobName $(pwd|awk -F '/' '{printf("%s/%s",$(NF-1),$NF)}')

seedname=$(find . -name "*.win"|cut -d '/' -f2|cut -d '.' -f1)
wannier_kpoints.py
echo scf
mpirun -np $SLURM_JOB_CPUS_PER_NODE pw.x -npool 4 <scf.in >scf.out
echo prepare
mpirun -np $SLURM_JOB_CPUS_PER_NODE wannier90.x -pp ${seedname}
echo pw2wan
mpirun -np 12 pw2wannier90.x < pw2wan.in > pw2wan.out
echo wannier90
mpirun -np $SLURM_JOB_CPUS_PER_NODE wannier90.x ${seedname}
