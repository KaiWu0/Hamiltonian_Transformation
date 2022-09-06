# Hamiltonian_Transformation
Introduction of Hamiltonian transformation(HT).
1. It is a first-principles band structure calculation method, which is fast, accurate, parameter-free and functional independent.
2. The codes in "fortran" directory is HT, which is implemented in Quantum ESPRESSO (QE).
3. The codes in "python" directory is eigenvalue transformation, which is used for Wannier interpolation.
4. HT supports OpenMP, does not support MPI now. I recommend to compile two versions of QE, one is MPI version and runs pw.x, another is OpenMP version and runs ht.x, since pw.x is quite slow when OpenMP is open.

Install HT (fortran):
1. Copy the QE directory to a new place. The HT is written based on QE-6.7, you may need to modify some codes of HT for newest QE version (see Compile Errors section).
2. Copy “ht.f90” and “ht.sh” to the new directory.
3. Run “bash ht.f90”, it will output “success”.
4. run configure. Add OpenMP flags in you configure parameters, i.e. “-enable-openmp”. This step is optional, ht.x will run slower without OpenMP.
   example: "./configure MPIF90=mpiifort CC=icc F77=ifort FC=ifort -enable-openmp"
5. Compile QE by “make clean; make pw pp -j”

HT Compile Errors:
1. If you see "No rule to make target `@spin_orb@', needed by `ht.o'.  Stop." That is because the new version of QE changed module name.
   Open ht.f90 file, search and replace "USE spin_orb" with "USE noncollin_module".

Run HT examples (fortran/example):
1. The examples contain silicon-PBE, bi2se3 and feb2 scripts.
2. You should change some parameters for your own environment.
3. The nscf directory contains script for non-scf calculation, ht directory contains script for Hamiltonian transformation, wannier contains script for wannier interpolation.
4. qe.sh or wannier.sh contains commands about how to run the script.
5. The eigenvalues calculated by HT are saved in band.txt as a matrix.
6. There are python scripts to plot band structures in silicon-PBE directory. Run plot_band.py to plot bands. You may need to modify some parameters.

Run eigenvalue transformation script (python):
1. Run "python transform.py [prefix of QE] forward" before the Wannier90 optimization.
2. Run "python transform.py [prefix of QE] backward" after the Wannier90 optimization.

HT Input Parameters:
1. The necessary parameters are prefix, outdir and K_POINTS block. You can obtain satisfactory results without other parameters.
2. funtype:
    Type of transform functions.
    0 is simple transform function in paper.
    1 is deprecated.
    2 is complex transform function in paper (default).
3. power, gap:
    power is n, gap is a of f_{a,n}(x) or f_n(x) in paper.
    Default: power = 3, gap = -1 (less than 0 to make program calculate automatically)
4. delete_top_bands:
    Whether delete top bands automatically, default is .false.
5. qr_eps:
    Controls the rank of QRCP. Assume the largest singular value to 1, and discard singular values smaller than qr_eps. Default 4e-3. Sometimes it can be reduced to 1e-2.
6. rank:
    Assign the rank to QRCP to a particular value. Default -1, less than 0 to make program calculate automatically.
7. num_threads:
    Control the threads of OpenMP, default 8. Set it larger will not make the program faster.
8. window_min, search_start, search_end:
    Discard bands lower than window_min. In fact, the inner logic is complicated, the program can search and set window_min automatically. You do not need to modify them.
9. read_eig:
    Bool value, default is .false.. Whether read in SCF eigenvalues from the specified file instead of QE default file, which is used to interpolate GW bands.
10. eig_file:
    The name of file to read in SCF eigenvalues. Only used when read_eig is .true..
