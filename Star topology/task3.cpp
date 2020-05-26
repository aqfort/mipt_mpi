#include "mpi.h"
#include <iostream>
#include <cstdlib>
#include <iomanip>

using namespace std;

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    MPI_Status STATUS;

    int size = 0;
    int rank = 0;
    
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(argc != (size + 1)) { //check does argv have proper number of arguments
        cout << "argc: invalid number of arguments" << endl;
        MPI_Finalize();
        exit(10);
    }

    int *index = new int[size]; //number of arcs proceeding from each vertex {2, 3, 4}
    index[0] = size - 1;
    for(int i = 1; i < size; i++) {
        index[i] = index[i - 1] + 1;
    }

    int *edges = new int[(size - 1) * 2]; //sequential list of graph arcs {1, 2, 0, 0}
    for(int i = 0; i < (size - 1); i++) {
        edges[i] = i + 1;
    }
    for(int i = (size - 1); i < (size - 1) * 2; i++) {
        edges[i] = 0;
    }

    MPI_Comm COMM_STAR; //communicator with graph type (star) topology
    MPI_Graph_create(MPI_COMM_WORLD, size, index, edges, 1, &COMM_STAR);

    delete[] index;
    delete[] edges;

    if(rank == 0) {
        int *DATA = new int[size]; //input data (array of numbers)

        for(int i = 1; i < size; i++) { //send DATA to other processes
            DATA[i] = atoi(argv[i + 1]);
            MPI_Send(&DATA[i], 1, MPI_INT, i, 10, COMM_STAR);
        }

        for(int i = 1; i < size; i++) { //get replies from other processes
            int ARG = 0;
            MPI_Recv(&ARG, 1, MPI_INT, i, 10, COMM_STAR, &STATUS);
            if(ARG != DATA[i]) {
                printf("ERROR\n");
            } else {
                printf("SUCCESS(%02d)\n", i);
            }
        }

        delete[] DATA;
    } else {
        int ARG = 0;
        MPI_Recv(&ARG, 1, MPI_INT, 0, 10, COMM_STAR, &STATUS);

        printf("DATA[%02d] = %d\n", rank, ARG);

        MPI_Send(&ARG, 1, MPI_INT, 0, 10, COMM_STAR);
    }

    MPI_Comm_free(&COMM_STAR);

    MPI_Finalize();

    return 0;
}