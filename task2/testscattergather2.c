#include<mpi.h>
#include<stdio.h>
#include<windows.h> // remove
#include<unistd.h>
#include<stdlib.h>
#include<time.h>

float* generateMatrix(int size, float* array){
  for(int i=0; i<size; i++){
      for(int j=0; j<size; j++){
        array[i * size + j] = i + (float)j/100;
      }
  }
  return array;
}



int main(int argc, char **argv){
  double start, end;
  MPI_Init(&argc, &argv);
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  printf("The number of processes are: %d\n",world_size);
  fflush(stdout);

  printf("This is processor: : %d\n",world_rank);
  fflush(stdout);


  int rows_per_proc=atoi(argv[1]) ? atoi(argv[1]) : 1;
  int row_size = rows_per_proc*world_size + 2;
  int elements_per_proc = rows_per_proc*row_size;
  
  float* matrix=malloc(row_size*row_size*sizeof(float));

  MPI_Barrier(MPI_COMM_WORLD); // remove
  start = MPI_Wtime();

  if(world_rank == 0) {
    generateMatrix(row_size, matrix);

  }
  

  int* elements_per_proc_list = malloc(sizeof(int)*world_size);
  int* recv_elements_per_proc_list = malloc(sizeof(int)*world_size);

  int* proc_displacement = malloc(sizeof(int)*world_size);
  int* recv_proc_displacement = malloc(sizeof(int)*world_size);

  if (world_rank == 0) {
    for (int p = 0; p < world_size; p++) {
        recv_elements_per_proc_list[p] = elements_per_proc;
        recv_proc_displacement[p] = elements_per_proc * p + row_size;
        elements_per_proc_list[p] = elements_per_proc + 2 * row_size;
        proc_displacement[p] = elements_per_proc * p;
    }

    printf("Starting Matrix: \n");
    for (int i = 0; i < row_size*row_size; i+=row_size ) {
      for (int j = 0; j < row_size; j++ ) {
        printf("%f, ", matrix[i + j]);
      }
      printf("\n");
    }
    printf("\n");
    fflush(stdout);

  }

  int elements_in_my_proc = elements_per_proc + 2 * row_size;
  float* subsection=malloc(sizeof(float)*elements_in_my_proc);
  float* newsubsection=malloc(sizeof(float)*elements_in_my_proc);
  float* topRow=malloc(sizeof(float)*row_size);
  float* bottomRow=malloc(sizeof(float)*row_size);
  float* swap;
  MPI_Barrier(MPI_COMM_WORLD); // remove
  MPI_Scatterv(matrix, elements_per_proc_list, proc_displacement, MPI_FLOAT, subsection, elements_in_my_proc, MPI_FLOAT, 0, MPI_COMM_WORLD);

  for (int i = 0; i < elements_in_my_proc; i+=row_size ) {
    for (int j = 0; j < row_size; j++ ) {
      newsubsection[i + j] = subsection[i + j];
    }
  }


  // MPI_Scatter(matrix, elements_per_proc, MPI_FLOAT, subsection, elements_in_my_proc, MPI_FLOAT, 0, MPI_COMM_WORLD);
  for (int u = 0; u < 5; u++) {


    printf("Before PROC %d:\n",world_rank);
    for (int i = 0; i < elements_in_my_proc; i+=row_size ) {
      for (int j = 0; j < row_size; j++ ) {

        printf("%f, ", subsection[i + j]);
      }
      printf("\n");
    }
    printf("\n");


    for (int i = row_size; i < elements_in_my_proc - row_size; i++ ) {
        newsubsection[i] = subsection[i] + subsection[i - row_size];   
    }

    printf("After PROC %d:\n",world_rank);

    for (int i = 0; i < elements_in_my_proc; i+=row_size ) {
      for (int j = 0; j < row_size; j++ ) {
        printf("%f, ", newsubsection[i + j]);
      }
      printf("\n");
    }
    printf("\n");
    fflush(stdout);
    if (world_rank % 2 == 0) {
      if (world_rank != 0) {
        MPI_Send(&(newsubsection[row_size]), row_size, MPI_FLOAT, world_rank-1, 0, MPI_COMM_WORLD);
        // MPI_Recv(newsubsection, row_size, MPI_FLOAT, world_rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(topRow, row_size, MPI_FLOAT, world_rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
      if (world_rank != world_size - 1) {
        MPI_Send(&(newsubsection[elements_in_my_proc - 2 * row_size]), row_size, MPI_FLOAT, world_rank+1, 0, MPI_COMM_WORLD);
        // MPI_Recv(&(newsubsection[elements_in_my_proc - row_size]), row_size, MPI_FLOAT, world_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(bottomRow, row_size, MPI_FLOAT, world_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
    }
    else {
      if (world_rank != world_size - 1) {
        // MPI_Recv(&(newsubsection[elements_in_my_proc - row_size]), row_size, MPI_FLOAT, world_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(bottomRow, row_size, MPI_FLOAT, world_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&(newsubsection[elements_in_my_proc - 2 * row_size]), row_size, MPI_FLOAT, world_rank+1, 0, MPI_COMM_WORLD);
      }
      if (world_rank != 0) {
        // MPI_Recv(newsubsection, row_size, MPI_FLOAT, world_rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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

  
  MPI_Gatherv(&(subsection[row_size]), elements_per_proc, MPI_FLOAT, matrix, recv_elements_per_proc_list, recv_proc_displacement, MPI_FLOAT, 0, MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD); // remove


  end = MPI_Wtime();
  if (world_rank == 0) {
    printf("Time taken: %lf seconds\n\n", end - start);
    fflush(stdout);
    printf("Final Matrix: \n");
    for (int i = 0; i < row_size*row_size; i+=row_size ) {
      for (int j = 0; j < row_size; j++ ) {
        printf("%f, ", matrix[i + j]);
      }
      printf("\n");
    }
    fflush(stdout);

  }

  MPI_Finalize();
  free(matrix);
  
}

