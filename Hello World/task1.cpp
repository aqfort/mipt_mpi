#include "mpi.h"
#include <iostream>

using namespace std;

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int size = 0;
    int rank = 0;
    
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    printf("hello from (rank = %02d) process\nsize = %02d\n", rank, size);
    
    MPI_Finalize();
    
    return 0;
}
