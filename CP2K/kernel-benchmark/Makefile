export FC = gfortran
export FCFLAGS = -g #-O3
export LDFLAGS = $(FCFLAGS)

all:
	( cd src; make )

clean:
	rm -f kernel_benchmark
	( cd src; make clean )
