#include "mpi.h"
#include <iostream>

using namespace std;

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  cout << "hello_friend\n";
  int size = 0;
  int rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  cout << "size = " << size << endl << "rank = " << rank << endl << endl;
  MPI_Finalize();
  return 0;
}
