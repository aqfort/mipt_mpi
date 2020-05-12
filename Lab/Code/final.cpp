#include "mpi.h"
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <fstream>

#define T 1 // t in [0; T]
#define X 1 // x in [0; X]

#define A 0.5

#define K 300
#define M 100

using namespace std;

// ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

double U_EXACT(const int &k,
               const int &m);

double F_EXACT(const int &k,
               const int &m);

double U_EXPLICIT(const int &m,
                  double* DATA_LINE,
                  double* fooo_LINE);

void PRINT(double* DATA);

// ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

// ROOT == 0

// ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);

    if(argc != 1) {
        cout << "No argument is needed!\n";
        MPI_Finalize();
        exit(1);
    }

    MPI_Status STATUS;

    int SIZE = 0;
    int RANK = 0;
    
    MPI_Comm_size(MPI_COMM_WORLD, &SIZE);
    MPI_Comm_rank(MPI_COMM_WORLD, &RANK);

// ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

    double* DATA_EXACT;
    double* DATA;

    double* fooo;

    double* SUBDATA_EXACT;
    double* SUBDATA;

    if(RANK == 0) {
        DATA_EXACT = (double*) calloc((K + 1) * (M + 1), sizeof(double));
        DATA = (double*) calloc((K + 1) * (M + 1), sizeof(double));
        fooo = (double*) calloc((K + 1) * (M + 1), sizeof(double));
    }

    int P = (M + 1) / SIZE; // x-axis splitting: P nodes for each process

    SUBDATA_EXACT = (double*) calloc(P, sizeof(double));
    SUBDATA = (double*) calloc(P, sizeof(double));

// ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

    // Initial conditions: U(0, x) = PHI(x)

    MPI_Barrier(MPI_COMM_WORLD);

    for(int p = 0; p < P; p++) {
        SUBDATA[p] = U_EXACT(0, P * RANK + p);
    }

    if(RANK == 0) {
        for(int p = (M + 1) - (M + 1) % SIZE; p <= M; p++) {
            DATA[p] = U_EXACT(0, p);
        }
    }

    MPI_Gather(SUBDATA, P, MPI_DOUBLE, DATA, P, MPI_DOUBLE, 0, MPI_COMM_WORLD);

// ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

    // Initial conditions: U(t, 0) = PSI(t), U(t, X) = PSI'(t), F(t, x)

    for(int k = 0; k <= K; k++) {
        MPI_Barrier(MPI_COMM_WORLD);

        for(int p = 0; p < P; p++) {
            SUBDATA_EXACT[p] = U_EXACT(k, P * RANK + p);
            SUBDATA[p] = F_EXACT(k, P * RANK + p); // F(t, x) -- the majority
        }

        if(RANK == 0) {
            for(int p = (M + 1) - (M + 1) % SIZE; p <= M; p++) {
                DATA_EXACT[(M + 1) * k + p] = U_EXACT(k, p);
                fooo[(M + 1) * k + p] = F_EXACT(k, p); // F(t, x) -- the rest
            }

            DATA[(M + 1) * k] = U_EXACT(k, 0); // Initial conditions: U(t,0) = PSI(t)
            DATA[(M + 1) * k + M] = U_EXACT(k, M); // Initial conditions: U(t, X) = PSI'(t)
        }

        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Gather(SUBDATA_EXACT, P, MPI_DOUBLE, &DATA_EXACT[(M + 1) * k], P, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(SUBDATA, P, MPI_DOUBLE, &fooo[(M + 1) * k], P, MPI_DOUBLE, 0, MPI_COMM_WORLD); // gathering F(t, x) in root process
    }

// ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

    double TIME_BEGIN = 0;
    double TIME_END = 0;

    if(RANK == 0) {
        TIME_BEGIN = MPI_Wtime();
    }

    // Explicit method: calculating U(t, x)

    free(SUBDATA);

    SUBDATA = (double*) calloc((M - 1), sizeof(double));

    double* DATA_LINE = (double*) calloc((M + 1), sizeof(double)); // Preceding line of DATA
    double* fooo_LINE = (double*) calloc((M + 1), sizeof(double)); // Preceding line of fooo

    P = (M - 1) / SIZE;

    for(int k = 1; k <= K; k++) {
        MPI_Barrier(MPI_COMM_WORLD);

        if(RANK == 0) {
            for(int m = 0; m <= M; m++ ) {
                DATA_LINE[m] = DATA[(M + 1) * (k - 1) + m];
                fooo_LINE[m] = fooo[(M + 1) * (k - 1) + m];
            }
        }

        MPI_Bcast(DATA_LINE, (M + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(fooo_LINE, (M + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);

        for(int p = 0; p < P; p++) {
            SUBDATA[p] = U_EXPLICIT(1 + P * RANK + p, DATA_LINE, fooo_LINE);
        }

        if(RANK == 0) {
            for(int p = (M - 1) - (M - 1) % SIZE + 1; p < M; p++) {
                DATA[(M + 1) * k + p] = U_EXPLICIT(p, DATA_LINE, fooo_LINE);
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Gather(SUBDATA, P, MPI_DOUBLE, &DATA[(M + 1) * k + 1], P, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

// ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

    if(RANK == 0) {
        // cout << "░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░\n";
        // cout << "░░░░░░░░░░░░░░░░░░░░░░░░░DATA_EXACT░░░░░░░░░░░░░░░░░░░░░░░░░\n";
        // cout << "░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░\n";
        // PRINT(DATA_EXACT);
        // cout << "░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░\n";
        // cout << "░░░░░░░░░░░░░░░░░░░░░░░░░░░░DATA░░░░░░░░░░░░░░░░░░░░░░░░░░░░\n";
        // cout << "░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░\n";
        // PRINT(DATA);
    }

// ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

    MPI_Barrier(MPI_COMM_WORLD);

    if(RANK == 0) {
        free(DATA_EXACT);
        free(DATA);
        free(fooo);

        TIME_END = MPI_Wtime();

        cout << "░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░\n";
        cout << "░░░░░░░░░░░░░░░░░░░░TIME SPENT: " << fixed << setprecision(7) << TIME_END - TIME_BEGIN << "...\n";
        cout << "░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░\n";
    }

    free(SUBDATA);
    free(SUBDATA_EXACT);

    free(DATA_LINE);
    free(fooo_LINE);

    MPI_Finalize();

    return 0;

}

// ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

double U_EXACT(const int &k,
               const int &m) {
    return sin(M_PI * (0.5 + k * T / ((double) K))) *
           sin(M_PI * m * X / ((double) M));
}

double F_EXACT(const int &k,
               const int &m) {
    return (U_EXACT(k + 1, m) - 0.5 * (U_EXACT(k, m + 1) + U_EXACT(k, m - 1))) * K / ((double) T) + 
           (U_EXACT(k, m + 1) - U_EXACT(k, m - 1)) * A * M * 0.5 / ((double) X);
}

double U_EXPLICIT(const int &m,
                  double* DATA_LINE,
                  double* fooo_LINE) {
    return 0.5 * (DATA_LINE[m + 1] + DATA_LINE[m - 1]) +
           (fooo_LINE[m] - (DATA_LINE[m + 1] - DATA_LINE[m - 1]) * A * M * 0.5 / ((double) X)) * T / ((double) K);
}

void PRINT(double* DATA) {
    long flag = cout.precision();
    cout << fixed << setprecision(7);

    for(int k = 0; k <= K; k++) {
        for(int m = 0; m <= M; m++) {
            cout << k << ' ' << m << ' ' << DATA[(M + 1) * k + m] << endl;
        }
        cout << endl;
    }

    cout.precision(flag);
    cout << resetiosflags(ios::fixed);
}