# Introduction
Code developed by Prof. Jeferson W. D. Fernandes

OMEGA - Object-oriented tool for Multi-physics Engineering Analyses

[![Run Unit Tests](https://github.com/jefersondossa/FSArl/actions/workflows/runUnitTests.yml/badge.svg)](https://github.com/jefersondossa/FSArl/actions/workflows/runUnitTests.yml)

# Usage
Perform Fluid-Structure interactions using the Arlequin method.

Once all the packages described in Presentation are provided, the software can be used as follows:

- To compile use the command 

	<CODE> make f </CODE>

- To run use

	<CODE> make run4(run8,run16,etc.) </CODE>

- To clear all files created running the program use

	<CODE> make clear </CODE>

- To visualize Eigen matrices and vectors while debugging with gdb, please do the following:

https://gitlab.com/libeigen/eigen/-/tree/master/debug/gdb

- PETSc installation for debug:
<CODE> sudo ./configure PETSC_ARCH=arch-linux2-c-debug --with-64-bit-indices --with-mpi-dir=/home/jeferson/MPICH/install/ --with-cxx-dialect=C++11 --with-debugging=1 --with-X=1 --download-mumps --download-scalapack --download-ptscotch --download-fblaslapack  --download-suitesparse </CODE>


- PETSc installation for Release:

<CODE> sudo ./configure PETSC_ARCH=arch-linux2-c-opt --with-mpi-dir=/home/jeferson/MPICH/install/ --with-cxx-dialect=C++11 --with-debugging=0 --with-X=1 COPTFLAGS='-O3 -march=native -mtune=native' CXXOPTFLAGS='-O3 -march=native -mtune=native' FOPTFLAGS='-O3 -march=native -mtune=native' --download-mumps --download-scalapack --download-ptscotch --download-fblaslapack --with-64-bit-indices --download-suitesparse   </CODE>

-Graphic Interface:

Provide the following libraries:

<CODE> sudo apt-get install libglfw3-dev libfltk1.3-dev libgl1-mesa-dev</CODE>