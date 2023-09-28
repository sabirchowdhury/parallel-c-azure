#include<mpi.h>
#include<stdio.h>
// #include<windows.h> // remove
#include<unistd.h>
#include<stdlib.h>
#include<time.h>

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
  printf("The number of processes are: %d\n",world_size);
  fflush(stdout);

  printf("This is processor: : %d\n",world_rank);
  fflush(stdout);


  int rows_per_proc=1;
  if(argc > 1){  if(argv[1] && atoi(argv[1]) > 0) {
    rows_per_proc = atoi(argv[1]);
  }}
  printf("%d",argc);
  fflush(stdout);
  int row_size = rows_per_proc*world_size + 2;
  int elements_per_proc = rows_per_proc*row_size;
  
  float* matrix=NULL;

  MPI_Barrier(MPI_COMM_WORLD); // remove
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
  MPI_Barrier(MPI_COMM_WORLD); // remove

  for (int u = 0; u < 38; u++) {
    MPI_Scatterv(matrix, elements_per_proc_list, proc_displacement, MPI_FLOAT, subsection, elements_in_my_proc, MPI_FLOAT, 0, MPI_COMM_WORLD);

    int pos;
    for (int i = row_size; i < elements_in_my_proc - row_size; i+=row_size ) {
        newsubsection[i] = subsection[i];
        newsubsection[i+row_size-1] = subsection[i+row_size-1];
        for (int j = 1; j < row_size - 1; j++ ) {
        pos = i+j;
        newsubsection[pos] = (subsection[pos+1] + subsection[pos-1] + subsection[pos+row_size] + subsection[pos-row_size])/4;
        }
    }
    MPI_Gatherv(&(newsubsection[row_size]), elements_per_proc, MPI_FLOAT, matrix, recv_elements_per_proc_list, recv_proc_displacement, MPI_FLOAT, 0, MPI_COMM_WORLD);

  }

  

  MPI_Barrier(MPI_COMM_WORLD); // remove


  end = MPI_Wtime();
  if (world_rank == 0) {

    // printf("Final Matrix: \n");
    // // int rowsPrint = row_size;
    // // int colsPrint = row_size;

    // int rowsPrint = row_size;
    // int colsPrint = row_size;

    // for (int i = 0; i < row_size*rowsPrint; i+=row_size ) {
    //   for (int j = 0; j < colsPrint; j++ ) {
    //     printf("%f, ", matrix[i + j]);
    //   }
    //   printf("\n");
    // }
    // fflush(stdout);
    printf("Time taken: %lf seconds\n\n", end - start);
    fflush(stdout);
  }

  MPI_Finalize();
  free(matrix);
  
}

