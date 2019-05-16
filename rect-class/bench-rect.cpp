#include "communication.h"
#include "rectmm.h"
#include "rectsizes.h"
#include "library.h"
#include <cmath>
#include <vector>
#include <tuple>
#include <limits>
#include <chrono>
#include <iostream>
#include <algorithm>
#include <sstream>

void fillInt( double *A, int n ) {
  for( int i = 0; i < n; i++ )
    A[i] = double((int)( 10*drand48()));
}

int round_up(int n) {
    return (int) std::pow(2, std::ceil(std::log2(n)));
}

int round_down(int n) {
    return (int) std::pow(2, std::floor(std::log2(n)));
}

long long compute_memory(int m, int n, int P) {
    return m * n / P;
}

long long total_memory(std::vector<long long>& memory) {
    return memory[0] + memory[1] + memory[2];
}

// index of the dimension to divide
long long new_memory(std::vector<long long> current_memory, std::vector<int> dims, int index) {
    long long total_current_memory = total_memory(current_memory);

    dims[index] /= 2;

    if (index == 0) {
        total_current_memory += 2 * current_memory[1];
    } else if (index == 1) {
        total_current_memory += 2 * current_memory[0];
    } else {
        total_current_memory += 2 * current_memory[2];
    }

    return total_current_memory;
}

// assumes m, n, k, P are powers of 2
std::tuple<std::string, std::vector<int>, int, int> get_carma_pattern(int m, int n, int k, 
        int P, long long memory_limit = std::numeric_limits<long long>::max()) {
    int r = 0;
    std::vector<int> dims = {m, n, k};
    std::string pattern;
    std::vector<int> divPattern;

    std::vector<long long> memory = {compute_memory(m, k, P),
        compute_memory(k, n, P),
        compute_memory(m, n, P)};

    // if some of the dimensions reaches 1 before number of processors reached P,
    // then we will use less processors
    int P_actual = 1;
    while (P_actual < P) {
        auto ptr_to_max = std::max_element(dims.begin(), dims.end());
        if (*ptr_to_max > 1) {
            int index = std::distance(dims.begin(), ptr_to_max);

            auto required_memory = new_memory(memory, dims, index);

            if (required_memory <= memory_limit) {
                // divide the largest dimension by 2
                dims[index] /= 2;
                // add BFS stcommep to the pattern
                pattern += "b";
                // the maximum dimension is divided by 2, others by 1
                for (int j = 0; j < 3; ++j) {
                    divPattern.push_back((index==j ? 2 : 1));
                }
                P_actual *= 2;
                r++;
            } else {
                // divide the largest dimension by 2
                dims[index] /= 2;
                // add DFS step to the pattern
                pattern += "d";
                // the maximum dimension is divided by 2, others by 1
                for (int j = 0; j < 3; ++j) {
                    divPattern.push_back((index==j ? 2 : 1));
                }
                memory[index] += 2 * memory[index];
                r++;
            }
        } else {
            return std::make_tuple(pattern, divPattern, P_actual, r);
        }
    }
    return std::make_tuple(pattern, divPattern, P, r);
}

int get_n_iter() {
    const char* value = std::getenv("n_iter");
    std::stringstream strValue;
    strValue << value;

    unsigned int intValue;
    strValue >> intValue;

    if (intValue<0 || intValue > 100) {
        std::cout << "Number of iteration must be in the interval [1, 100]" << std::endl;
        std::cout << "Setting it to 1 iteration instead" << std::endl;
        return 1;
    }

    return intValue;
}

int main( int argc, char **argv ) {
  initCommunication(&argc, &argv);

  int rank = getRank();

  double *A, *B, *C;

  int m = read_int( argc, argv, "-m", 128 );
  int n = read_int( argc, argv, "-n", 128 );
  int k = read_int( argc, argv, "-k", 128 );
  m = round_up(m);
  n = round_up(n);
  k = round_up(k);

  long long memory_limit = read_long_long( argc, argv, "-L", -1);
  if (memory_limit < 0) {
      memory_limit = std::numeric_limits<long long>::max();
  }

  int n_iter = get_n_iter();

  int r;
  std::string pattern;
  std::vector<int> divPattern;
  int P_actual;

  int P_orig;
  MPI_Comm_size( MPI_COMM_WORLD, &P_orig );
  int P = round_down(P_orig);

  std::tie(pattern, divPattern, P_actual, r) = get_carma_pattern(m, n, k, P, memory_limit);

  if( getRank() == 0 ) {
    printf("Benchmarking %dx%dx%d multiplication using %d processes\n", m,n,k,P);
    printf("Division pattern: ");
    for( int i = 0; i < r; i++ ) {
      for( int j = 0; j < 3; j++ )
        printf("%d,",divPattern[3*i+j]); 
      printf(" ");
    }
    printf("execution pattern: %s\n", pattern);
  }

  if (rank < P) {
      initSizesRect( m, n, k, P, r, divPattern.data() );

      // allocate the initial matrices
      A = (double*) malloc( sizeof(double)*getSizeRect(m,1,k,P) );
      B = (double*) malloc( sizeof(double)*getSizeRect(1,n,k,P) );
      C = (double*) malloc( sizeof(double)*getSizeRect(m,n,1,P) );

      // fill the matrices with random data
      srand48(getRank());
  }

  std::vector<long long> times;

  for (int i = 0; i < n_iter; ++i) {
      if (rank < P) {
          fillInt( A, getSizeRect(m,1,k,P) );
          fillInt( B, getSizeRect(1,n,k,P) );
          fillInt( C, getSizeRect(m,n,1,P));
      }

      MPI_Barrier( MPI_COMM_WORLD );
      auto start = std::chrono::steady_clock::now();
      if (rank < P) {
          rectMM( A, B, C, m, n, k, P, r, &pattern[0], divPattern.data() );
      }
      MPI_Barrier( MPI_COMM_WORLD );
      auto end = std::chrono::steady_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
      if (rank == 0) {
          times.push_back(duration);
      }
  }

  if( getRank() == 0 ) {
    std::cout << "OLD_CARMA TIMES [ms] = ";
    for (int i = 0; i < n_iter; ++i) {
        std::cout << times[i] << " ";
    }
    std::cout << std::endl;
  }

  if (rank < P) {
      free(A);
      free(B);
      free(C);
  }

  MPI_Finalize();
  return 0;
}
