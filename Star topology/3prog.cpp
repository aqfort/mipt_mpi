#include "mpi.h"
#include <iostream>
#include <cstdlib>
#include <iomanip>

using namespace std;

int main(int argc, char **argv) {

  MPI_Init(&argc, &argv);

  MPI_Status status;

  int size = 0;
  int rank = 0;
  
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(argc != (size + 1)) { //check if argv has proper number of arguments
    cout << "argc: invalid number of arguments" << endl;
    MPI_Finalize();
    exit(10);
  }

  int *index = (int *) malloc(sizeof(int) * size); //number of arcs proceeding from each vertex {2, 1, 1}
  for(int i = 1; i < size; i++) {
    index[i] = 1;
  }
  index[0] = size - 1;

  int *edges = (int *) calloc((size - 1) * 2, sizeof(int)); //sequential list of graph arcs {1, 2, 0, 0}
  for(int i = 0; i < (size - 1); i++) {
    edges[i] = i + 1;
  }

  MPI_Comm COMM_STAR; //communicator with graph type (star) topology
  MPI_Graph_create(MPI_COMM_WORLD, size, index, edges, 1, &COMM_STAR);

  free(index);
  free(edges);

  if(rank == 0) {
    int *DATA = (int *) malloc(sizeof(int) * size); //input data (array of numbers)
    for(int i = 0; i < size; i++) { //send DATA to other processes
      DATA[i] = atoi(argv[i + 1]); // + 1; //DATA (+1)
      MPI_Send(&DATA[i], 1, MPI_INT, i, 10, COMM_STAR);
    }
    for(int i = 1; i < size; i++) { //get replies from other processes
      MPI_Recv(&DATA[i], 1, MPI_INT, i, 10, COMM_STAR, &status);
      if(DATA[i] != 100) {
        cout << "ERROR" << endl;
      } else {
        cout << "SUCCESS(" << setw(2) << setfill('0') << i << ")" << endl;
      }
    }
    //DATA[0]++; //DATA (+2) only root
    //for(int i = 0; i < size; i++) {//print result (should be DATA + 2)
    //  cout << "DATA[" << i << "] = " << DATA[i] << endl;
    //}
    //for(int i = 0; i < size + 1; i++) {//print result (should be DATA + 2)
    //  cout << "argv[" << i << "] = " << argv[i] << endl;
    //}
    free(DATA);
  } else {
    int ARG = 0;
    MPI_Recv(&ARG, 1, MPI_INT, 0, 10, COMM_STAR, &status);
    //ARG++; //DATA (+2) except root
    cout << "DATA[" << setw(2) << setfill('0') << rank << "] = " << ARG << endl;
    ARG = 100;
    MPI_Send(&ARG, 1, MPI_INT, 0, 10, COMM_STAR);
  }

  //int A = 0;
  //if(rank == 0) {
  //  A = 10;
  //}
  //MPI_Bcast(&A, 1, MPI_INT, 0, COMM_STAR);
  //cout << "A (rank = " << rank << ") = " << A << endl;

  MPI_Comm_free(&COMM_STAR);

  MPI_Finalize();

  return 0;
}