#include<mpi.h>
#include<stdio.h>
#include<unistd.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include <string.h>

//000000000000000000000000000000000000000000000000000000000000000000000000000000 - 80 columns

// Generate matrix takes in an array size and a pointer to the matrix
// uses the pointer to modify the original matrix and sets up the initial state
void generateMatrix(int size, float* matrix);

// Helper to display the matrix given column width and number of elements
void displayMatrix(float* matrix, int columns, int nElements);

// This is a test function: to display the allocation of the rows
// Arguments take in both send and recieve displacement (beginning position)
// and the number of elements. Also take the number of columns (matrix width)
// and the processor number (for printf/displaying purposes)
void testAllocateRows(int sendDisplacement,
                      int sendElementsPerProc,
                      int recvDisplacement,
                      int recvElementsPerProc,
                      int nColumns,
                      int p);




int main(int argc, char **argv) {



  double startTime, endTime;
  int rank;
  int nProc;
  float precision = 0.01;
  // We are using a float to store values, this limits us to 7 decimal
  // places of accuracy, to increase this we could switch to using
  // doubles or even long doubles, but I've decided not to, so that my
  // matrix uses less space at large matrix sizes and can allocate memory.
  int nColumns = 20;
  int displayLevel = 0; // comment about what the levels mean

  // Initialise MPI, use library to extract the rank and total processors
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProc);
  // Get current time as a start time stamp used later to measure time
  startTime = MPI_Wtime(); ////////////////////////////////////// Should we do this differently ///////////////////////////////////////////////////



  // This for loop sweeps through all the arguments and collects the values
  // from the flags. To set a value from terminal you would use a flag:
  // e.g. -size 2
  // The algorithm below sniffs all the arguments for any recognised flags
  // and sets the value after it to the internal variable

  // e.g. mpiexec filename -size 10 -precision 0.01 -display 2
  for (int argIndex = 1; argIndex < argc - 1; argIndex++ ) {
      if (!strcmp(argv[argIndex], "-size")) {
          nColumns = atoi(argv[argIndex + 1]);
      }
      if (!strcmp(argv[argIndex], "-precision")) {
          float precisionArg = (float) atof(argv[argIndex + 1]);
          // if the precision is invalid i.e. negative, set it to default 0.01
          precision = precisionArg < 0 ? 0.01 : precisionArg;
      }
      if (!strcmp(argv[argIndex], "-display")) {
          displayLevel = atoi(argv[argIndex + 1]);
      }
  }

  // Unit testing for arguments and checking all processors are communicating
  if (displayLevel > 0) {
    if (rank == 0) {
      printf("Arguments: size: %d, precision %f, display level: %d\n",
            nColumns, precision, displayLevel);
      printf("The number of processes are: %d\n",nProc);
      fflush(stdout); // Allows prints to happen in parallel instead of serially
    }
    // Barriers ensures synchronisation in prints, groups them together
    // and makes it easier to read
    MPI_Barrier(MPI_COMM_WORLD);

    printf("This is processor: : %d\n",rank);
    fflush(stdout);

    MPI_Barrier(MPI_COMM_WORLD);


  }

  // Validation statement to ensure there are enough rows for 1 per processor
  // note -2 is to ignore the border pixels as part of calculations
  nColumns = nProc < nColumns - 2 ? nColumns : nProc + 2;

  int rowsPerProc = floor((nColumns - 2)/nProc);
  int nProcsWithExtraRow = (nColumns - 2) % nProc;
  
//000000000000000000000000000000000000000000000000000000000000000000000000000000 - 80 columns

  // these are variables initialised for scatterv and gatherv they will be
  // arrays that store the no of elements and displacement for each core. They
  // are only needed in the control processor and initialised empty to begin
  // with the memeory allocate is within the if statement below to save memory

  int* sendElementsPerProc;
  int* recvElementsPerProc;
  int* sendDisplacement;
  int* recvDisplacement;

  // Q: Why do we have different variables (displacement + no elements)
  // for both the send and recieve?
  // A: The send and the recieve variables will have different displacements
  // and number of elements the reason is that when sending out each division
  // of the matrix to each processor to calculate the neighbour avaerage the
  // processor needs the rows above and below when we recieve all the data at
  // then end, we only care about the non border values.
  // Hence for the differnt sizes we create extra variables for the recv

  // Processor rank 0 is our control processor
  if (rank == 0) {
    // memory allocation only happens to rank 0 to save memory
    // allocating an integer for each processor (nProc)
    sendElementsPerProc = malloc(sizeof(int)*nProc);
    recvElementsPerProc = malloc(sizeof(int)*nProc);

    sendDisplacement = malloc(sizeof(int)*nProc);
    recvDisplacement = malloc(sizeof(int)*nProc);

    if (!sendElementsPerProc || !recvElementsPerProc ||
        !sendDisplacement || !recvDisplacement){
        printf("Error: Memory allocation failed");
        exit(1);
    }

    // for loop, where p represents the processor rank
    for (int p = 0; p < nProc; p++) {
      // to divide up ONLY WHEN rows do not perfectly divide the array size
      // figures if code should take/add an extra row to balance out workload
      int rowsAssignedToProcessorP = p < nProcsWithExtraRow ?
                                          (rowsPerProc + 1) : rowsPerProc;

      // Calculate row displacements:
      // As some processors will take an extra row, we need to take this into
      // account when thinking about displacement. If the processor rank
      // is less than the "nProcsWithExtraRow" number of processors that need to
      // take an extra row, than the processor will take an extra row. That
      // being said, it also means all n "p" processors before it also took an
      // extra row, so the extra displacement is "p". In the case where
      // "p" > "nProcsWithExtraRow", the extra displacment is "nProcsWithExtraRow"

      int extraRowsFromPrevProcs = p < nProcsWithExtraRow ? p : nProcsWithExtraRow;
      int rowDisplacement = p * rowsPerProc + extraRowsFromPrevProcs;

      // + 2 for upper neighbour row and lower neighbour row
      sendElementsPerProc[p] = nColumns * (rowsAssignedToProcessorP + 2);
      recvElementsPerProc[p] = nColumns * rowsAssignedToProcessorP;

      sendDisplacement[p] = nColumns * (rowDisplacement);
      recvDisplacement[p] = nColumns * (rowDisplacement + 1);
      // We add 1 to the row displacement to ignore the top border and move all
      // calculated rows down by 1.

      // the send uses two extra rows, the upper and lower, these extra rows
      // are required for the above and below neighbour averaging in relaxation

      // Unit testing
      // checking the level of information to display for testing
      if (displayLevel > 1) {
        testAllocateRows(sendDisplacement[p],
                        sendElementsPerProc[p],
                        recvDisplacement[p],
                        recvElementsPerProc[p],
                        nColumns, p);
        // Flush the print buffer so outputs happen in relative order
        fflush(stdout);
      }

    }

  }




//000000000000000000000000000000000000000000000000000000000000000000000000000000 - 80 columns

  // variable to store matrix, initially set empty to save memory
  // only to be set on processor 0
  float* matrix;
  if(rank == 0) {
    // allocate memory square matrix (columns = rows) (columns x columns)
    matrix=malloc(nColumns*nColumns*sizeof(float));

    // checking if memory allocation passed
    if (!matrix) {
      printf("Error: Memory allocation failed - matrix");
      exit(1);
    }

    // sends pointer of matrix to generate matrix to get values filled in
    generateMatrix(nColumns, matrix);

    // if display level is 2 or greater show the matrix at the start
    if (displayLevel > 1) {
      printf("\nStarting Matrix:\n");
      displayMatrix(matrix, nColumns, nColumns*nColumns);
      printf("---------------------------------------------\n\n");
      fflush(stdout);
    }
  }



  // similarly to above, we are calculating how many rows this processor is
  // designated, conditional statement checks if it needs to take an extra row
  // These variables are for scatter and gather communication count
  // number of rows are multiplied by number of columns (width)
  int rowsInProc = (rank < nProcsWithExtraRow) ?
                                 (rowsPerProc + 1) : rowsPerProc; 

  // Number of elements being edited by this processor (i.e. number it will send
  // back to control processor (gather) once reached precision)
  int nElementsToEdit = rowsInProc * nColumns;

  // Number of elements being sent to the processor from matrix, these are the
  // rows it was assigned + 2 extra rows for the top and bottom neighbour rows
  int elementsInProc = (rowsInProc + 2) * nColumns;

  // subsection is the array that each processor will store thier share of the
  // matrix
  float* subsection=malloc(sizeof(float)*elementsInProc);
  // newSubsection is the resultant array where calculations are stored
  float* newSubsection=malloc(sizeof(float)*elementsInProc);
  // This is readability purposes but is an array at the size of a row,
  // used as a communication buffer to store the top row on recieve
  float* topRow=malloc(sizeof(float)*nColumns);
  // This is readability purposes but is an array at the size of a row,
  // used as a communication buffer to store the bottom row on recieve
  float* bottomRow=malloc(sizeof(float)*nColumns);
  // float pointer created as a temporary medium to swap two pointers
  float* pointerSwap;

  // checking all memory allocations succeeded
  if (!subsection || !newSubsection ||
      !topRow || !bottomRow) {
    printf("Error: Memory allocation failed - subsection\n");
    printf("Proc: %d, elements in proc %d\n", rank, elementsInProc);
    printf("subsection %p newsubsection %p top row %p bottom row %p",subsection,newSubsection,topRow,bottomRow);
    exit(1);
  }


  if (displayLevel > 0) {
    // Barriers ensures synchronisation in prints, groups them together
    // and makes it easier to read - barrier will not be run in normal execution
    MPI_Barrier(MPI_COMM_WORLD);
    printf("Proc %d: Elements in proc: %d\n",rank, elementsInProc);
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

  }

  // Uses all logic described above to distribute the shares of the matrix to
  // processors. scatterV is used as the shares aren't equal, blocking
  // communication allows synchronisation and for all processors to start
  // together and avoid race conditions.
  // This scatter only needs to happen once (hence outside the while loop)
  // as we only communicate neighbour borders between processors (further down),
  // not the entire sections of the matrix.
  MPI_Scatterv(matrix, sendElementsPerProc, sendDisplacement, MPI_FLOAT, 
                subsection, elementsInProc, MPI_FLOAT, 0, MPI_COMM_WORLD);

  // Above the scatter recieved its section of the matrix in "subsection", we
  // need to copy the same values into newSubsection to ensure correctness when
  // swapping the matrix sections. This only needs to happen once at the start
  // as the algorithm updates the values from here on.
  // This is quicker than creating an identical scatter but just with
  // "newSubsection" as the recieve buffer.
  for (int i = 0; i < elementsInProc; i+=nColumns ) {
    for (int j = 0; j < nColumns; j++ ) {
      newSubsection[i + j] = subsection[i + j];
    }
  }

  // This is used as a flag/buffer to store the precision of all other
  // processors, when this is 0, it tells us that all processors are precise
  int globalPrecision = 1;

  // This is a local precision flag, (at 0 processor is precise), this will
  // be sent to all other processors to indiciate the precision for this
  // processor
  int localPrecision = 0;

  // Stores the iteration count (how many times we have had to do relaxation)
  int iterations = 0;

  // A variable used to store the linear coordinates (for readability)
  int pos;

  // linear element position for start of bottom border
  int bottomBorderPosition = elementsInProc - nColumns;



//000000000000000000000000000000000000000000000000000000000000000000000000000000 - 80 columns

  // Keep doing the iteration as long as atleast one processor is not precise
  while (globalPrecision > 0) {
    iterations++;

    // For Testing: Displays section of matrix for this processor
    if (displayLevel > 2) {
      printf("%d) Processor %d:\n", iterations, rank);
      displayMatrix(subsection, nColumns, elementsInProc);
      printf("\n");
      fflush(stdout);
    }


    // RELAXATION TECHNIQUE:

    // This is the same as 1 for loop with an if statement inside to see
    // if the cell is not a border value, but I chose this way as I find it
    // more readable. Furthermore it runs slightly quicker than running
    // a "single for loop linearly with an if statement" as this method doesnt
    // require an if statement to be run for EACH individual cell 
    for (int i = nColumns; i < bottomBorderPosition; i+=nColumns ) {
        for (int j = 1; j < nColumns - 1; j++ ) {
        pos = i+j;
        newSubsection[pos] = (subsection[pos+1] + 
                              subsection[pos-1] + 
                              subsection[pos+nColumns] + 
                              subsection[pos-nColumns]) / 4;
        }
    }

    // reset local precision flag to precise
    localPrecision = 0;
    // Similar for loop style to for loop above, goes through each
    // cell and compares the result matrix to the pre result matrix
    // Breaks instantly when the loop finds a cell that is not precise
    // to save iterations
    // this flag checking could have been joined with the loop above, however
    // we would miss out on the time saving of breaking the loop early
    // hence i split this out in to two loops
    for (int i = nColumns; i < bottomBorderPosition; i+=nColumns ) {
      for (int j = 1; j < nColumns - 1; j++ ) {
        pos = i+j;
        if (fabs(newSubsection[pos] - subsection[pos]) > precision) {
          localPrecision = 1;
          break;
        }
      }
      if (localPrecision != 0) {
        break;
      }
    }
//000000000000000000000000000000000000000000000000000000000000000000000000000000 - 80 columns

    // This is blocking communication and acts as a synchronisation barrier
    // The purpose this serves is that it communicates with all other threads,
    // sums all the local precision values and adds it to the global precision
    // flag. All reduce does this in one step
    // Summing all the flags works here as its 0 for precise and >0 for not
    // precise, hence if all cores are precise, then the sum would be 0 and we
    // can exit
    MPI_Allreduce(&localPrecision, &globalPrecision, 1,
                   MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    // For Testing: Displays section of matrix for this processor
    if (displayLevel > 3) {
      printf("%d) (result) Processor %d:\n", iterations, rank);
      displayMatrix(newSubsection, nColumns, elementsInProc);
      printf("\n");
      fflush(stdout);
    }

    // Exits the loop if all cores are precise, (explained above)
    if (globalPrecision == 0) {
      break;
    }


//000000000000000000000000000000000000000000000000000000000000000000000000000000 - 80 columns


    // Send and recieve the updated neighbour rows:

    // Instead of giving each processor an updated version of the whole matrix
    // we can optimise and reduce communication to only recieve the neighbouring
    // rows. This processor needs the values of the row above to calculate the
    // average for its top rows, and the values of the row below to calulate the
    // neighbour avergages at the bottom row.

    //              |___...___|     _________ 
    //  subsection  |_________| -> |_________|  subsection
    //   proc: n    |_________| <- |___...___|   proc: n+1 
    //                  //         |         |

    // To ensure we do not get a deadlock we use even processors to send first
    // and recieve next, and odd processors to recieve first and send next
    // NB: all if statements skip if there is only one processor

    // split out even and odd processors
    if (rank % 2 == 0) {
      // Send top row to prev proc to be recieved as its new bottom row
      // Recieve new top row from the prev proc (bottom row)
      // (skip rank 0, row above is matrix border i.e. constant)
      if (rank != 0) {
        MPI_Send(&(newSubsection[nColumns]), nColumns,
                  MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD);

        MPI_Recv(topRow, nColumns,
                  MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }

      // Send bottom row to next proc to be recieved as its new top row
      // Recieve new bottom row from the next proc (top row)
      // (skip end processor, row below is matrix border i.e. constant)
      if (rank != nProc - 1) {
        MPI_Send(&(newSubsection[elementsInProc - 2 * nColumns]), nColumns,
                  MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD);

        MPI_Recv(bottomRow, nColumns,
                  MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
    }
    else {
      // Recieve new bottom row from the next proc (top row)
      // Send bottom row to next proc to be recieved as its new top row
      // (skip end processor, row below is matrix border i.e. constant)
      if (rank != nProc - 1) {
        MPI_Recv(bottomRow, nColumns,
                  MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Send(&(newSubsection[elementsInProc - 2 * nColumns]), nColumns,
                  MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD);
      }

      // Recieve new top row from the prev proc (bottom row)
      // Send top row to prev proc to be recieved as its new bottom row
      // (skip rank 0, row above is matrix border)
      if (rank != 0) {
        MPI_Recv(topRow, nColumns,
                  MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Send(&(newSubsection[nColumns]), nColumns,
                  MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD);
      }
    }


    // Copy the new bottom and top rows to both our matrix sections
    // this is required 
    // NB: all if statements skip if there is only one processor
    if (rank != 0) {
      for (int col = 0; col < nColumns; col++) {
        newSubsection[col] = topRow[col];
        subsection[col] = topRow[col];
      }
    }
    if (rank != nProc - 1) {
      for (int col = 0; col < nColumns; col++) {
        newSubsection[col + bottomBorderPosition] = bottomRow[col];
        subsection[col + bottomBorderPosition] = bottomRow[col];
      }
    }
    pointerSwap = newSubsection;
    newSubsection = subsection;
    subsection = pointerSwap;

    

//000000000000000000000000000000000000000000000000000000000000000000000000000000 - 80 columns

  }
  


  MPI_Gatherv(&(newSubsection[nColumns]), nElementsToEdit, MPI_FLOAT,
              matrix, recvElementsPerProc, recvDisplacement, MPI_FLOAT,
                                                      0, MPI_COMM_WORLD);
  endTime = MPI_Wtime();
  if (rank == 0) {
    printf("Iterations: %d\n", iterations);
    printf("Time taken: %lf seconds\n\n", endTime - startTime);
    if (displayLevel > 1) {
      printf("---------------------------------------------\n\n");
      printf("\nFinal Matrix:\n");
      displayMatrix(matrix, nColumns, nColumns*nColumns);
      printf("\n");
    }
    fflush(stdout);
  }

  MPI_Finalize();
  free(matrix);
  free(sendElementsPerProc);
  free(recvElementsPerProc);
  free(sendDisplacement);
  free(recvDisplacement);
  free(subsection);
  free(newSubsection);
  free(topRow);
  free(bottomRow);
  
}

// ------------------------------------------------------
// ------------------ Help functions --------------------
// ------------------------------------------------------


void generateMatrix(int size, float* matrix){
    for (int x = 0; x < size; x++) {
        matrix[x] = 1;                  // top border
        matrix[x + (size-1)*size] = 1;  // bottom border
        matrix[size * x] =1;            // left border
        matrix[size-1 + size * x] = 1;  // right border
    }
}


// ------------------------------------------------------
// --------------------- Testing ------------------------
// ------------------------------------------------------


void displayMatrix(float* matrix, int columns, int nElements) {
  for (int i = 0; i < nElements; i+=columns ) {
    for (int j = 0; j < columns; j++ ) {
      printf("%f, ", matrix[i + j]);
    }
    printf("\n");
  }
  printf("\n");
  fflush(stdout);
}

void testAllocateRows(int sendDisplacement,
                      int sendElementsPerProc,
                      int recvDisplacement,
                      int recvElementsPerProc,
                      int nColumns,
                      int p) {
  // Unit testing

  // We are using the simple rule that
  // (x = linear position % width)
  // (y = linear position // width)

  // the linear start position is the displacement
  // the linear end position is the displacement + number of elements

  int sendX = sendDisplacement % nColumns;
  int sendY = floor(sendDisplacement / nColumns);
  int sendXEnd = (sendDisplacement + sendElementsPerProc) % nColumns;
  int sendYEnd = floor((sendDisplacement + sendElementsPerProc) / nColumns);
  int recvX = recvDisplacement% nColumns;
  int recvY = floor(recvDisplacement / nColumns);
  int recvXEnd = (recvDisplacement + recvElementsPerProc) % nColumns;
  int recvYEnd = floor((recvDisplacement + recvElementsPerProc) / nColumns);


  printf("Processor %d send starts at (x,y) (%d,%d) "
         "and ends at (x,y) (%d,%d), rows: %d \n",
         p, sendX, sendY, sendXEnd, sendYEnd, sendYEnd-sendY );

  printf("Processor %d recv starts at (x,y) (%d,%d) "
         "and ends at (x,y) (%d,%d), rows: %d \n\n",
         p, recvX, recvY, recvXEnd, recvYEnd, recvYEnd-recvY );
}
//000000000000000000000000000000000000000000000000000000000000000000000000000000 - 80 columns
