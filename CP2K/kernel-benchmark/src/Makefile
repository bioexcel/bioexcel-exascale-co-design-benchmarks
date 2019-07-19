all: build_dependencies qmmm_gpw_forces.f90 kernel_benchmark.f90
	$(FC) $(FCFLAGS) -I dependencies -I dependencies/base -I dependencies/subsys -c qmmm_gpw_forces.f90
	$(FC) $(FCFLAGS) -I dependencies -I dependencies/base -I dependencies/subsys -c kernel_benchmark.f90
	$(FC) $(LDFLAGS) kernel_benchmark.o qmmm_gpw_forces.o dependencies/subsys/cell_types.o dependencies/base/kinds.o dependencies/qmmm_util.o dependencies/common/math*.o

build_dependencies:
	( cd dependencies; make )

clean:
	rm -f *.o *.mod
	( cd dependencies; make clean )

