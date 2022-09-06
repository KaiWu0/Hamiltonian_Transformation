#!/bin/bash
#SBATCH -J QE
#SBATCH -o out
#SBATCH -e err
#SBATCH -N 1 -n 48
#SBATCH -p long
scontrol update job $SLURM_JOB_ID JobName $(pwd|awk -F '/' '{printf("%s/%s",$(NF-1),$NF)}')

seedname=$(find . -name "*.win"|cut -d '/' -f2|cut -d '.' -f1)
#echo scf
#cp -r scf_for_wan/bi2se3.save .
echo prepare
mpirun -np 48 wannier90.x -pp ${seedname}
echo pw2wan
mpirun -np 4 pw2wannier90.x < pw2wan.in > pw2wan.out
echo wannier90
mpirun -np 48 wannier90.x ${seedname}
