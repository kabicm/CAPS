# The following assumes that the module system has been used to load the GNU compiler
# mpi c++ compiler, with whatever flag is required for openmp
CC=CC -fopenmp
# whatever options are necessary to link against the desired version of BLAS
LIBS=-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -lgomp -lpthread -lm -ldl

CFLAGS= -O3 -Wall

# set to the number of cores to be used per process
THREADS= -DNUM_THREADS=${OMP_NUM_THREADS}

# uncomment for sanity checks to abort on invalid input; may slightly affect performance
#SANITY= -DSANITY_CHECKS 
# uncomment to use daxpy instead of openmp for additions; the best choice seems to be machine dependent
#DAXPY= -DDAXPY

DFLAGS = $(THREADS) $(SANITY) $(DAXPY)

OFLAGS = $(CFLAGS) $(DFLAGS)
BFLAGS = $(CFLAGS) $(DFLAGS) $(LIBS)
