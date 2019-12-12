MANSEC           = KSP
CLEANFILES       = rhs.vtk solution.vtk
NP               = 1
CSOURCES 				 = $(wildcard *.cpp)
FCOMPILER        = gfortran -O2
CXXFLAGS        += -w

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test

f: modules.o porticomb.o $(CSOURCES:.cpp=.o) 
	@-${CLINKER} -o $@ $^ ${PETSC_KSP_LIB} -lboost_system -lgfortran -std=c++0x

debug: $(CSOURCES:.cpp=.o)
	@-${CLINKER} -o $@ $^ ${PETSC_KSP_LIB} -g
	@gdb debug

modules.o: modules.for
	@ ${FCOMPILER} -c modules.for

porticomb.o: porticomb.for
	@ ${FCOMPILER} -c porticomb.for

clear:
	@$ rm *.o *~ f *.vtu mirror* domain* *.mod *.dat ma26* tensao* esforc* saida omega.txt *.msh *.geo

run1:
	@$ mpirun -np 1 ./f

run2:
	@$ mpirun -np 2 ./f

run3:
	@$ mpirun -np 3 ./f

run4:
	@$ export OMP_NUM_THREADS=1 && mpirun -n 4 ./f

run5:
	@$ mpirun -np 5 ./f -info -logsummary

run6:
	@$ export OMP_NUM_THREADS=1 && mpirun -n 6 ./f

run7:
	@$ mpirun -np 7 ./f

run8:
	@$ mpirun -np 8 ./f -pc_factor_nonzeros_along_diagonal 1.e-8

run16:
	@$ mpirun -np 16 ./f -pc_type jacobi -ksp_type gmres -ksp_monitor_singular_value -ksp_gmres_restart 1000
