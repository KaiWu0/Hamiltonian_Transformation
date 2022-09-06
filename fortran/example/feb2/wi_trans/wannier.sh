#!/bin/bash
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

seedname=$(find . -name "*.win"|cut -d '/' -f2|cut -d '.' -f1)
echo scf
mpirun -np 24 pw.of -npool 4 < scf.in >scf.out
echo prepare
wannier_kpoints.py
mpirun -np 12 wannier90.o -pp ${seedname}
echo pw2wan
mpirun -np 4 /public/home/jlyang/wuk/bin/qe-7.0/bin/pw2wannier90.x < pw2wan.in > pw2wan.out
echo wannier90
./transform.py feb2 forward
mpirun -np 12 wannier90.o ${seedname}
./transform.py feb2 backward
