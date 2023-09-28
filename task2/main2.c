#include<mpi.h>
#include<stdio.h>
// #include<windows.h> // remove
#include<unistd.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
//000000000000000000000000000000000000000000000000000000000000000000000000000000 - 80 columns
void generateMatrix(int size, float* matrix){
    for (int x = 0; x < size; x++) {
        matrix[x] = 1;
        matrix[x + (size-1)*size] = 1;
        matrix[size * x] =1;
        matrix[size-1 + size * x] = 1;
    }
}



int main(int argc, char **argv){
  double start, end;
  MPI_Init(&argc, &argv);
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  float precision = 0.01;
  int rows_per_proc=1;
  if(argc > 1){  if(argv[1] && atoi(argv[1]) > 0) {
    rows_per_proc = atoi(argv[1]);
  }}

  int row_size = rows_per_proc*world_size + 2;
  int elements_per_proc = rows_per_proc*row_size;
  
  float* matrix;

  start = MPI_Wtime();

  if(world_rank == 0) {
    matrix=malloc(row_size*row_size*sizeof(float));
    generateMatrix(row_size, matrix);


  }
  

  int* elements_per_proc_list = NULL;
  int* recv_elements_per_proc_list = NULL;

  int* proc_displacement = NULL;
  int* recv_proc_displacement = NULL;

  if (world_rank == 0) {
    elements_per_proc_list = malloc(sizeof(int)*world_size);
    recv_elements_per_proc_list = malloc(sizeof(int)*world_size);
    proc_displacement = malloc(sizeof(int)*world_size);
    recv_proc_displacement = malloc(sizeof(int)*world_size);
    for (int p = 0; p < world_size; p++) {

        recv_elements_per_proc_list[p] = elements_per_proc;
        recv_proc_displacement[p] = elements_per_proc * p + row_size;
        elements_per_proc_list[p] = elements_per_proc + 2 * row_size;
        proc_displacement[p] = elements_per_proc * p;
    }

    // printf("Starting Matrix: \n");
    // for (int i = 0; i < row_size*row_size; i+=row_size ) {
    //   for (int j = 0; j < row_size; j++ ) {
    //     printf("%f, ", matrix[i + j]);
    //   }
    //   printf("\n");
    // }
    // printf("\n");
    // fflush(stdout);

  }

  int elements_in_my_proc = elements_per_proc + 2 * row_size;
  float* subsection=malloc(sizeof(float)*elements_in_my_proc);
  float* newsubsection=malloc(sizeof(float)*elements_in_my_proc);
  float* topRow=malloc(sizeof(float)*row_size);
  float* bottomRow=malloc(sizeof(float)*row_size);
  float* swap;

  
  MPI_Scatterv(matrix, elements_per_proc_list, proc_displacement, MPI_FLOAT, subsection, elements_in_my_proc, MPI_FLOAT, 0, MPI_COMM_WORLD);

  for (int i = 0; i < elements_in_my_proc; i+=row_size ) {
    for (int j = 0; j < row_size; j++ ) {
      newsubsection[i + j] = subsection[i + j];
    }
  }

  int globalPrecision = 1;
  int localPrecision = 0;
  int iterations = 0;
  int pos; // A variable used to store the calculation to work out precision (Readability)


  // MPI_Scatter(matrix, elements_per_proc, MPI_FLOAT, subsection, elements_in_my_proc, MPI_FLOAT, 0, MPI_COMM_WORLD);
  while (globalPrecision > 0) {
    iterations++;

    for (int i = row_size; i < elements_in_my_proc - row_size; i+=row_size ) {
        for (int j = 1; j < row_size - 1; j++ ) {
        pos = i+j;
        newsubsection[pos] = (subsection[pos+1] + subsection[pos-1] + subsection[pos+row_size] + subsection[pos-row_size])/4;
        }
    }




    localPrecision = 0;
    for (int i = row_size; (localPrecision == 0) && (i < elements_in_my_proc - row_size); i+=row_size ) {
      for (int j = 1; j < row_size - 1; j++ ) {
        pos = i+j;
        if (fabs(newsubsection[pos] - subsection[pos]) > precision) {
          localPrecision = 1;
          break;
        }
      }
    }
    MPI_Allreduce(&localPrecision, &globalPrecision, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    // printf("%d: PROC: %d: local: %d, global %d\n",iterations, world_rank, localPrecision, globalPrecision);
    // fflush(stdout);


    // printf("%d: PROC: %d:\n",iterations, world_rank);
    // for (int i = world_rank == 0 ? 0 : row_size; i < ( world_rank != world_size - 1 ? elements_in_my_proc - row_size : elements_in_my_proc); i+=row_size ) {
    //   for (int j = 0; j < row_size; j++ ) {
    //     printf("%f, ", newsubsection[i + j]);
    //   }
    //   printf("\n");
    // }
    // printf("\n");
    // fflush(stdout);

    if (globalPrecision == 0) {
      break;
    }







    if (world_rank % 2 == 0) {
      if (world_rank != 0) {
        MPI_Send(&(newsubsection[row_size]), row_size, MPI_FLOAT, world_rank-1, 0, MPI_COMM_WORLD);
        MPI_Recv(topRow, row_size, MPI_FLOAT, world_rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
      if (world_rank != world_size - 1) {
        MPI_Send(&(newsubsection[elements_in_my_proc - 2 * row_size]), row_size, MPI_FLOAT, world_rank+1, 0, MPI_COMM_WORLD);
        MPI_Recv(bottomRow, row_size, MPI_FLOAT, world_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
    }
    else {
      if (world_rank != world_size - 1) {
        MPI_Recv(bottomRow, row_size, MPI_FLOAT, world_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&(newsubsection[elements_in_my_proc - 2 * row_size]), row_size, MPI_FLOAT, world_rank+1, 0, MPI_COMM_WORLD);
      }
      if (world_rank != 0) {
        MPI_Recv(topRow, row_size, MPI_FLOAT, world_rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&(newsubsection[row_size]), row_size, MPI_FLOAT, world_rank-1, 0, MPI_COMM_WORLD);
      }
      
    }

    if (world_rank != 0) {
      for (int col = 0; col < row_size; col++) {
        newsubsection[col] = topRow[col];
        subsection[col] = topRow[col];
      }
    }
    if (world_rank != world_size - 1) {
      for (int col = 0; col < row_size; col++) {
        newsubsection[col + elements_in_my_proc - row_size] = bottomRow[col];
        subsection[col + elements_in_my_proc - row_size] = bottomRow[col];
      }
    }
    swap = newsubsection;
    newsubsection = subsection;
    subsection = swap;

    


  }
  
  MPI_Gatherv(&(newsubsection[row_size]), elements_per_proc, MPI_FLOAT, matrix, recv_elements_per_proc_list, recv_proc_displacement, MPI_FLOAT, 0, MPI_COMM_WORLD);

  end = MPI_Wtime();
  if (world_rank == 0) {
    printf("Time taken: %lf seconds\n\n", end - start);
    fflush(stdout);
  }
  MPI_Finalize();
  free(matrix);
  
}

