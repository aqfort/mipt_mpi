#include "mpi.h"
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <sstream>

using namespace std;

void my_error(const string &str) { //catch errors
    cout << str << endl;
    exit(10);
}

int main(int argc, char **argv) { //calculating partial sum of numbers from 1 to N
    MPI_Init(&argc, &argv);

    MPI_Status STATUS;

    if(argc < 2) {
        my_error("argc: not enough arguments");
    }

    long N = atol(argv[1]);

    double time_0 = MPI_Wtime(); //begining of whole program (timer)

    int size = 0;
    int rank = 0;
    
    if(MPI_Comm_size(MPI_COMM_WORLD, &size) != MPI_SUCCESS) {
        MPI_Finalize();
        my_error("MPI_Comm_size error");
    }
    if(MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS) {
        MPI_Finalize();
        my_error("MPI_Comm_rank error");
    }

    long A = N / ((long) (size - 1)); //each process [0; (size - 1)) sums A elements, process (size - 1) sums the rest then finds result
    long partial_sum = 0; //sum of A elements

    if(rank == (size - 1)) { //root
        long result = 0;

        for(long i = ((long) rank) * A + ((long) 1); i <= N; i++) { //calculating the rest
            result += i;
        }

        for(int i = 0; i < (size - 1); i++) { //calculating the result
            MPI_Recv(&partial_sum, 1, MPI_LONG, i, 10, MPI_COMM_WORLD, &STATUS);
            result += partial_sum;
        }

        cout << "----- ----- ----- ----- -----" << endl;
        cout << A << " elements (in each process)" << endl;
        cout << "size = " << size << endl;
        cout << "result: sum of " << N << " numbers = " << result << endl;
        cout << "----- ----- ----- ----- -----" << endl;
    }
    else{
        long partial_end = (((long) rank) + (long) 1) * A;

        for(long i = ((long) rank) * A + ((long) 1); i <= partial_end; i++) { 
            partial_sum += i;
        }

        MPI_Send(&partial_sum, 1, MPI_LONG, (size - 1), 10, MPI_COMM_WORLD);
    }

    double time_spent = MPI_Wtime() - time_0;





    MPI_Barrier(MPI_COMM_WORLD);





    if(rank == (size - 1)) { //timer stops at root
        MPI_Barrier(MPI_COMM_WORLD);
        printf("time spent = %.5lf\n", time_spent);
    } else {
        printf("time(%02d) = %.5lf\n", rank, time_spent);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Finalize();

    return 0;
}
