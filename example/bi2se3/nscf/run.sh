mpirun -np 24 pw.x -npool 2 < scf.in >scf.out
mpirun -np 24 pw.x -npool 2 < band.in >band.out
#mpirun -np 1 ht.x < band.in >band.out
