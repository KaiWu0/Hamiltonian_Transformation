# Hamiltonian_Transformation
Introduction of Hamiltonian transformation(HT).
1. It is a first-principles band structure calculation method, which is fast, accurate, parameter-free and functional independent.
2. The codes in "fortran" directory is HT, which is implemented in Quantum ESPRESSO (QE).
3. The codes in "python" directory is eigenvalue transformation, which is used for Wannier interpolation.
4. HT supports OpenMP, does not support MPI now. We recommend to compile two versions of QE, one is MPI version and runs pw.x, another is OpenMP version and runs ht.x, since pw.x is quite slow when OpenMP is open.

Install HT (fortran):
1. Download QE-7.1 or higer version from official website (*do not use the QE codes in github*).
2. Copy “ht.f90” and “ht.sh” to the main directory of QE.
3. Run `bash ht.sh`, it will copy files to the correct position and change Makefiles, then output “success”.
4. run configure. Add OpenMP flags in you configure parameters, i.e. `-enable-openmp`. This step is optional, ht.x will run slower without OpenMP. Examples:
   - Intel compiler: `./configure MPIF90=mpiifort CC=icc F77=ifort FC=ifort --with-scalapack=intel -enable-openmp`.
   - Gcc compiler: `./configure MPIF90=mpif90 CC=gcc F77=gfortran FC=gfortran --with-scalapack=intel -enable-openmp`. If you use intel scalapack and gcc compiler, search and replace "lmkl_blacs_intelmpi_lp64" with "lmkl_blacs_openmpi_lp64" in make.inc.
7. Compile QE by `make clean; make pw pp -j`

Run HT examples (fortran/example):
1. The examples contain silicon-PBE, bi2se3 and feb2 scripts.
2. You should change some parameters for your own environment.
3. The nscf directory contains script for non-scf calculation, ht directory contains script for Hamiltonian transformation, wannier contains script for wannier interpolation.
4. qe.sh or wannier.sh contains commands about how to run the script. The `pw.x` and `ht.x` commands are renamed (`pw.*` and `ht.*`) just to distinguish different versions.
5. The eigenvalues calculated by HT are saved in band.txt as a matrix.
6. There are python scripts to plot band structures in silicon-PBE directory. Run plot_band.py to plot bands. You may need to modify some parameters.
7. Silicon-GW only contains ht directory. For Wannier interpolating GW band structures, see the silicon example in BerkelyGW, you can use eigenvalue transformation script in the calculation.

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
    Controls the rank of QRCP. Assume the largest singular value to 1, and discard singular values smaller than qr_eps. Default 1e-3. Sometimes it can be reduced to 1e-2.
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
11. max_rank:
    This parameter is used in randomized QRCP. If max_rank is set to 0 (default) in input script, the code will set max_rank=20*number_of_bands. If max_rank <0, max_rank=abs(max_rank)*number_of_bands. If ht.x throw an error with max_rank is too small, set max_rank=-50 is usually enough.
