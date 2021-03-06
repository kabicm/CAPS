#include "outer.h"
#include "testlib.h"
#include "local-multiply.h"
#include "types.h"
#include "generate.h"
#include "library.h"
#include <mpi.h>

int main( int argc, char **argv ) {
  MPI_Init( &argc, &argv );
  int rank, P;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &P );

  int s = read_int( argc, argv, "-s", 10 );
  int d = read_int( argc, argv, "-d", 1 );
  int ds = read_int( argc, argv, "-ds", 0 );
  int_d n = ((int_d) 1) << s;
  int v = read_int( argc, argv, "-v", 0 );
  int_d colsPerProc = (n+P-1)/P;
  int_d colsThisProc = min(colsPerProc,n-rank*colsPerProc);
  double density = 1.*d/n/(1<<ds);

  BlockIndexRowMajor bi = BlockIndexRowMajor(colsPerProc*rank, 0, colsThisProc, n );
  BlockIndexColMajor ai = BlockIndexColMajor(0, colsPerProc*rank, n, colsThisProc );

  Matrix *A = generateMatrix( &ai, density, time(0)+100*rank );
  Matrix *B = generateMatrix( &bi, density, time(0)+100*rank+1 );
  Matrix *C;
  vector<double> times;
  MPI_Barrier( MPI_COMM_WORLD );
  double startTime = read_timer();
  if( v ) {
    C = outerProduct( A, B, n, P, &times );    
  } else {
    C = outerProduct( A, B, n, P );
  }
  MPI_Barrier( MPI_COMM_WORLD );
  double stopTime = read_timer();
  if( rank == 0 ) {
    printf("n %ld d %d/(2^%d) P %d time %f matrix sizes this proc %lu %lu %lu\n", n, d, ds, P, stopTime-startTime, A->size(), B->size(), C->size());
  }
  if( v ) {
    double totalTimes[times.size()];
    double minTimes[times.size()];
    double maxTimes[times.size()];
    MPI_Reduce( times.data(), totalTimes, times.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( times.data(), minTimes, times.size(), MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
    MPI_Reduce( times.data(), maxTimes, times.size(), MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    if( rank == 0 ) {
      for( unsigned i = 0; i < times.size(); i++ )
	printf("%f (%f-%f)\n", totalTimes[i]/P,minTimes[i],maxTimes[i]);
    }
  }

  MPI_Finalize();
}
