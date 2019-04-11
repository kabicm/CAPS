# The following assumes that the module system has been used to load the GNU compiler
# mpi c++ compiler, with whatever flag is required for openmp
CC=CC -fopenmp
# whatever options are necessary to link against the desired version of BLAS
LIBS=-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -lgomp -lpthread -lm -ldl /users/kabicm/local/mpip/lib/libmpiP.so -lm -liberty

CFLAGS= -O3 -Wall

# set to the number of cores to be used per process
THREADS= -DNUM_THREADS=${OMP_NUM_THREADS}

TARGETS=bench-rect bench-rect-nc

all:	$(TARGETS)

bench-rect:	bench-rect.cpp communication.o library.o rectsizes.o rectmm.o blas.h communication.h rectmm.h library.h rectsizes.h sizes.o
	$(CC) $(CFLAGS) bench-rect.cpp -o bench-rect communication.o library.o rectsizes.o rectmm.o sizes.o $(LIBS)

bench-rect-nc:	bench-rect.cpp communication.o library.o rectsizes.o rectmm-nc.o blas.h communication.h rectmm.h library.h rectsizes.h sizes.o
	$(CC) $(CFLAGS) bench-rect.cpp -o bench-rect-nc communication.o library.o rectsizes.o rectmm-nc.o sizes.o $(LIBS)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $<

sequential:	sequential.cpp
	$(CC) sequential.cpp -o sequential $(LIBS)

clean:
	rm $(TARGETS) *.o
