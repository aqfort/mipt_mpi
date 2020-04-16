#include "mpi.h"
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector>

using namespace std;

//each function calculates elements of proper sequence and saves them in array

void SEQUENCE_E(const long &SIZE, double *ELEMENTS) { //e (mathematical constant)
  if(SIZE != 0) {
    ELEMENTS[0] = (double) 1;
  }
  for(long i = 1; i < SIZE; i++) {
    ELEMENTS[i] = ELEMENTS[i - 1] / ((double) i);
  }
}

void SEQUENCE_PI(const long &SIZE, double *ELEMENTS) { //pi (mathematical constant)
  for(long i = 0; i < SIZE; i++) {
    ELEMENTS[i] = ((double) 4) / ((double) (2 * i + 1));
    if(i % 2 != 0) {
      ELEMENTS[i] = -ELEMENTS[i];
    }
  }
}

void SEQUENCE_PI_FASTER(const long &SIZE, double *ELEMENTS) { //pi (mathematical constant)
  if(SIZE != 0) {
    ELEMENTS[0] = 3;
  }
  for(long i = 1; i < SIZE; i++) {
    ELEMENTS[i] = ((double) 4) / ((double) ((2 * i) * (2 * i + 1) * (2 * (i + 1))));
    if(i % 2 == 0) {
      ELEMENTS[i] = -ELEMENTS[i];
    }
  }
}

void SEQUENCE_1(const long &SIZE, double *ELEMENTS) { //sum of ln(1 + 3/(n^2))
  for(long i = 0; i < SIZE; i++) {
    ELEMENTS[i] = log(1 + ((double) 3) / ((double) pow((i + 1), 2)));
  }
}

int main(int argc, char **argv) {

  MPI_Init(&argc, &argv);

  double time_begin = MPI_Wtime();

  if(argc != 2) {
    cout << "argc: invalid number of arguments" << endl;
    MPI_Finalize();
    exit(1);
  }

  MPI_Status status;

  int _size = 0;
  int _rank = 0;
  
  MPI_Comm_size(MPI_COMM_WORLD, &_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &_rank);

  long SIZE = atol(argv[1]); //gets number of elements
  long SUBSIZE = SIZE / _size; //number of elements for each process

  if(_rank == 0) {
    cout << "----- ----- ----- ----- -----" << endl;
    cout << SIZE << " members in each sequence" << endl;
    cout << "----- ----- ----- ----- -----" << endl;
  }

  double *ELEMENTS = (double *) malloc(sizeof(double) * SIZE); //array of elements
  double *ELEMENTS_PI = (double *) malloc(sizeof(double) * SIZE);
  double *ELEMENTS_PI_FASTER = (double *) malloc(sizeof(double) * SIZE);
  double *ELEMENTS_1 = (double *) malloc(sizeof(double) * SIZE);
  double *SUBELEMENTS = (double *) malloc(sizeof(double) * SUBSIZE); //elements in each process
  double *SUBELEMENTS_PI = (double *) malloc(sizeof(double) * SUBSIZE);
  double *SUBELEMENTS_PI_FASTER = (double *) malloc(sizeof(double) * SUBSIZE);
  double *SUBELEMENTS_1 = (double *) malloc(sizeof(double) * SUBSIZE);

  double SUBSUM = 0; //sum of elements in each process
  double SUBSUM_PI = 0, SUBSUM_PI_FASTER = 0, SUBSUM_1 = 0;
  double SUM = 0; //resultunt sum
  double SUM_PI = 0, SUM_PI_FASTER = 0, SUM_1 = 0;

  if(_rank == 0) {
    SEQUENCE_E(SIZE, ELEMENTS); //needed function is to be written here
    SEQUENCE_PI(SIZE, ELEMENTS_PI);
    SEQUENCE_PI_FASTER(SIZE, ELEMENTS_PI_FASTER);
    SEQUENCE_1(SIZE, ELEMENTS_1);
    // for(long i = 0; i < SIZE; i++) {
    //   cout << ELEMENTS[i] << endl;
    // }
    // cout << endl;
  }

  MPI_Scatter(ELEMENTS, SUBSIZE, MPI_DOUBLE, SUBELEMENTS, SUBSIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD); //distributes elements to processes
  MPI_Scatter(ELEMENTS_PI, SUBSIZE, MPI_DOUBLE, SUBELEMENTS_PI, SUBSIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatter(ELEMENTS_PI_FASTER, SUBSIZE, MPI_DOUBLE, SUBELEMENTS_PI_FASTER, SUBSIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatter(ELEMENTS_1, SUBSIZE, MPI_DOUBLE, SUBELEMENTS_1, SUBSIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // cout << "here(" << _rank << "): " << endl;
  // for(long i = 0; i < SUBSIZE; i++) {
  //   cout << SUBELEMENTS[i] << endl;
  // }
  // cout << endl;

  for(long i = 0; i < SUBSIZE; i++) { //subsum in each process
    SUBSUM += SUBELEMENTS[i];
    SUBSUM_PI += SUBELEMENTS_PI[i];
    SUBSUM_PI_FASTER += SUBELEMENTS_PI_FASTER[i];
    SUBSUM_1 += SUBELEMENTS_1[i];
  }

  if(_rank == 0) { //sum of rest elements
    for(long i = SUBSIZE * _size; i < SIZE; i++) {
      SUBSUM += ELEMENTS[i];
      SUBSUM_PI += ELEMENTS_PI[i];
      SUBSUM_PI_FASTER += ELEMENTS_PI_FASTER[i];
      SUBSUM_1 += ELEMENTS_1[i];
    }
  }

  MPI_Reduce(&SUBSUM, &SUM, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); //finds sum of subsums calculated in each process in root (0)
  MPI_Reduce(&SUBSUM_PI, &SUM_PI, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&SUBSUM_PI_FASTER, &SUM_PI_FASTER, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&SUBSUM_1, &SUM_1, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  // cout << "subsum(" << _rank << ") = " << SUBSUM << endl;

  if(_rank == 0) {
    //cout << "result = " << fixed << setprecision(15) << SUM << endl;
    cout << fixed << setprecision(15);
    cout << "e:\n\t" << SUM << endl;
    cout << "pi:\n\t" << SUM_PI << endl;
    cout << "pi:\n\t" << SUM_PI_FASTER << endl;
    cout << "sum of ln(1 + 3/(n^2)):\n\t" << SUM_1 << endl;
  }

  cout << "time(" << setw(2) << setfill('0') << _rank << "): " << fixed << setprecision(5) << MPI_Wtime() - time_begin << endl;

  free(ELEMENTS);
  free(SUBELEMENTS);
  free(ELEMENTS_PI);
  free(SUBELEMENTS_PI);
  free(ELEMENTS_PI_FASTER);
  free(SUBELEMENTS_PI_FASTER);
  free(ELEMENTS_1);
  free(SUBELEMENTS_1);

  MPI_Finalize();

  return 0;
}