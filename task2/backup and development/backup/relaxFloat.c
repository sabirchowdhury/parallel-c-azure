#include<mpi.h>
#include<stdio.h>
#include<unistd.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include <string.h>
#include <time.h>

//000000000000000000000000000000000000000000000000000000000000000000000000000000 - 80 columns

// Generate matrix takes in an array size and a pointer to the matrix
// uses the pointer to modify the original matrix and sets up the initial state
void generateMatrix(int size, float* matrix);

// Helper to display the matrix given column width and number of elements
void displayMatrix(float* matrix, int columns, int nElements);

// This is a test function: to display the row allocations of each processor.
// Arguments take in both send and receive displacement (beginning position)
// and the number of elements. Also take the number of columns (matrix width)
// and the processor number/rank
void testAllocateRows(int sendDisplacement,
                      int sendElementsPerProc,
                      int recvDisplacement,
                      int recvElementsPerProc,
                      int nColumns,
                      int p);



// We are using a float to store values, this limits us to 7 decimal
// places of accuracy, to increase this we could switch to using
// doubles or even long doubles, but I've decided not to, so that my
// matrix uses less space at large matrix sizes and hence avoid failures
// in memory allocations. This allows us to test much larger problems.


int main(int argc, char **argv) {

  // Variables to store the times - used to calculate processing time
  double startTime, endTime;
  int rank; // Rank of this processor
  int nProc; // World size/total number of processors
  int nColumns = 20; // Array width/height default set to 20

  // Variable for storing the precision, default is set to 0.01.
  float precision = 0.01f;

  // Variable used for testing, so we can unit test and print out the inner
  // workings of the program and extra stats.
  int displayLevel = 0;
  // Display level meanings:
  // Level 1: shows input arguments, no of processors and each processor
  // communicates that they are alive
  // Level 2: Level 1 + the allocation split for each processor - both from
  // the control processor POV and each individual processors POV.
  // Also prints full starting matrix and the full final matrix
  // Level 3: Level 2 + displays the "pre-result" section (allocated rows)
  // delegated to each processor at each iteration.
  // Level 4: Level 3 + displays the result section matrix for each
  // processor at each iteration.

  // Initialise MPI, use library to extract the rank and total processors
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProc);
  // Get current time as a start time stamp used later to measure time
  startTime = MPI_Wtime();


  // This for loop sweeps through all the arguments and collects the values
  // from the flags. To set a value from terminal you would use a flag:
  // e.g. -size 2
  // The algorithm below sniffs all the arguments for any recognised flags
  // and sets the value after it to the corresponding variable

  // e.g. mpiexec filename -size 10 -precision 0.01 -display 2
  for (int argIndex = 1; argIndex < argc - 1; argIndex++ ) {
      if (!strcmp(argv[argIndex], "-size")) {
          nColumns = atoi(argv[argIndex + 1]);
      }
      if (!strcmp(argv[argIndex], "-precision")) {
          float precisionArg = (float) atof(argv[argIndex + 1]);
          // if the precision is invalid i.e. negative, set it to default 0.01
          precision = precisionArg < 0 ? 0.01f : precisionArg;
      }
      if (!strcmp(argv[argIndex], "-display")) {
          displayLevel = atoi(argv[argIndex + 1]);
      }
  }

  // Validation statement to ensure there are enough rows for 1 per
  // processor note "- 2" is to ignore the top & bottom border rows.
  nColumns = nProc < nColumns - 2 ? nColumns : nProc + 2;


  // Unit testing for arguments and checking all processors are communicating
  if (displayLevel > 0) {
    if (rank == 0) {
      printf("Arguments: size: %d, precision %f, display level: %d\n",
            nColumns, precision, displayLevel);
      printf("The number of processes are: %d\n",nProc);
      fflush(stdout); // Allows prints to happen in parallel instead of serially
    }
    // Barriers ensures synchronisation in prints, groups them together
    // and makes it easier to read - without the display flag these barriers
    // wont run and impede performance.
    MPI_Barrier(MPI_COMM_WORLD);

    printf("This is processor: : %d\n",rank);
    fflush(stdout);

    MPI_Barrier(MPI_COMM_WORLD);


  }


  // variable to store matrix, initially set empty to save memory
  // memory allocation only should be set on the control processor (rank 0)
  float* matrix;
  if(rank == 0) {
    // Allocate memory square matrix (columns = rows) (columns x columns)
    // Using calloc because we want it to be 0 initialised
    matrix=calloc((unsigned int) (nColumns*nColumns), sizeof(float));

    // checking if memory allocation passed
    if (!matrix) {
      printf("Error: Memory allocation failed - matrix");
      MPI_Finalize();
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


  // Calculate the base amount of rows per processor
  int rowsPerProc = (int) floor((nColumns - 2)/nProc);

  // To share uneven divisions some of the processors take an extra row. This
  // variable tells us how many remainder rows there are. If there are 5
  // remaining rows, 5 processors will take an extra row, these will be the
  // first 5 processors. Hence the name, number of processors with an extra row
  int nProcsWithExtraRow = (nColumns - 2) % nProc;

  // these are variables initialised for scatterv and gatherv they will be
  // arrays that store the no of elements and displacement for each core. They
  // are initialised empty to begin with. The if statement below ensure memory
  // allocation only happens for the processor with rank 0
  // NB: these variables are only specific/local to the control processor, as
  // it will be controlling the gather and scatter. These aren't to be available
  // to all processors (++ saves a small amount of memory + computation)

  int* sendElementsPerProc;
  int* recvElementsPerProc;
  int* sendDisplacement;
  int* recvDisplacement;

  // Q: Why do we have different variables (displacement + no elements)
  // for both the send and receive?
  // A: The send and the receive variables will have different displacements
  // and number of elements the reason is that when sending out each division
  // of the matrix to each processor to calculate the neighbour avaerage the
  // processor needs the rows above and below when we receive all the data at
  // then end, we only care about the non border values.
  // Hence for the differnt sizes we create extra variables for the recv

  // Processor rank 0 is our control processor
  if (rank == 0) {
    // memory allocation only happens to rank 0 to save memory
    // We need an integer value for each processor (nProc)
    sendElementsPerProc = malloc(sizeof(int) * (unsigned int) nProc);
    recvElementsPerProc = malloc(sizeof(int) * (unsigned int) nProc);

    sendDisplacement = malloc(sizeof(int) * (unsigned int) nProc);
    recvDisplacement = malloc(sizeof(int) * (unsigned int) nProc);

    // Checking if memory allocation failed
    if (!sendElementsPerProc || !recvElementsPerProc ||
        !sendDisplacement || !recvDisplacement){
        printf("Error: Memory allocation failed");
        MPI_Finalize();
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



  // similarly to above, we are calculating how many rows this processor is
  // designated, conditional statement checks if it needs to take an extra row
  // These variables are for scatter and gather communication count
  // number of rows are multiplied by number of columns (width)
  // NB: these variables below are LOCAL to the each processor, and tells iteslf
  // how many elements it should control. 
  int rowsInProc = (rank < nProcsWithExtraRow) ?
                                 (rowsPerProc + 1) : rowsPerProc; 

  // Number of elements being edited by this processor (i.e. number it will send
  // back to control processor (in the gather) once reached precision)
  int nElementsToEdit = rowsInProc * nColumns;

  // No. of elements being sent to this processor from original matrix, these are
  // the rows it was assigned + 2 extra rows for the top&bottom neighbour rows
  int elementsInProc = (rowsInProc + 2) * nColumns;

  // subsection is the array that each processor will store thier share of the
  // matrix
  float* subsection=malloc(sizeof(float) * (unsigned int) elementsInProc);
  // newSubsection is the resultant array where calculations are stored
  float* newSubsection=malloc(sizeof(float) * (unsigned int) elementsInProc);
  // float pointer created as a temporary medium to swap two pointers
  float* pointerSwap;

  // checking all memory allocations succeeded
  if (!subsection || !newSubsection) {
    printf("Error: Memory allocation failed - subsection\n");
    printf("Proc: %d, elements in proc %d\n", rank, elementsInProc);
    printf("subsection %p newsubsection %p",subsection,newSubsection);
    MPI_Finalize();
    exit(1);
  }


  if (displayLevel > 1) {
    // Barriers ensures synchronisation in prints, groups them together
    // and makes it easier to read - barrier will not be run in normal execution
    MPI_Barrier(MPI_COMM_WORLD);
    printf("Proc %d: Elements in proc: %d\n",rank, elementsInProc);
    // Flush ensures prints happen in parallel
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

  }

  // Scatter to distribute sections of the matrix to each processor:

  // Uses all logic described above to distribute the shares of the matrix to
  // processors. scatterV is used as the shares aren't equal, blocking
  // communication allows synchronisation and for all processors to start
  // together and avoid race conditions i.e. stops it from starting computation
  // before its even recieved any matrix data.
  // This scatter only needs to happen once (hence outside the while loop)
  // as we only communicate neighbour borders between processors (further down),
  // not the entire sections of the matrix.
  MPI_Scatterv(matrix, sendElementsPerProc, sendDisplacement, MPI_FLOAT, 
                subsection, elementsInProc, MPI_FLOAT, 0, MPI_COMM_WORLD);

  // Above the scatter received its section of the matrix in "subsection", we
  // need to copy the same values into newSubsection to ensure correctness when
  // swapping the matrix sections (explained later). This only needs to happen
  // once at the start - the algorithm keeps the values in sync from here on.
  // NB: This is quicker than creating an identical scatter but just with
  // "newSubsection" as the receive buffer.
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

  // A temporary variable used to store the linear coordinates (for readability)
  int pos;

  // The linear element position for start of bottom border (for readability)
  // also means we don't have to keep recalculating.
  int bottomBorderPosition = elementsInProc - nColumns;


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
    // more readable. Furthermore it runs slightly quicker than running a
    // "single for loop linearly with an if statement to avoid borders" as
    // this method doesnt require running an if statement for EACH individual
    // cell millions of times.
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
    // this flag checking could have been joined with the relaxation loop above,
    // however we would miss out on the time saving of breaking the loop early
    // hence i split this out in to two seperate nested loops
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

    // This is blocking communication and acts as a synchronisation barrier
    // The purpose this serves is that it communicates with all other threads,
    // sums all the local precision values and adds it to the global precision
    // flag. "Allreduce" does this in one step/line of code.
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


    // Communicate adjacent rows between adjacent processprs:
    // i.e. send and receive the updated border rows

    // Instead of giving each processor an updated version of the whole matrix
    // we can optimise and reduce communication to only receive the adjacent
    // rows. This processor needs the values of the row above to calculate the
    // average for its top row, and the values of the row below to calulate the
    // neighbour avergages at the bottom row.

    //              |___...___|     _________ 
    //  subsection  |_________| -> |_________|  subsection
    //   proc: n    |_________| <- |___...___|   proc: n+1 
    //                  //         |         |

    // To ensure we do not get a deadlock we use even processors to send first
    // and receive next, and odd processors to receive first and send next
    // The communication is also blocking so we can safely not worry about
    // the processor relooping before it has received its updated borders
    // NB: all if statements skip if there is only one processor

    // Variables here are created for readability
    float* topBorderRow = newSubsection;
    float* topRow = newSubsection + nColumns;
    float* bottomRow = newSubsection + elementsInProc - 2 * nColumns;
    float* bottomBorderRow = newSubsection + bottomBorderPosition;

    // split out even and odd processors
    if (rank % 2 == 0) {
      // Send top row to prev proc to be received as its new bottom border row
      // receive new top border row from the prev proc (bottom row)
      // (skip if rank 0, row above is matrix border i.e. constant)
      if (rank != 0) {
        MPI_Send(topRow, nColumns,
                  MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD);

        MPI_Recv(topBorderRow, nColumns,
                  MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }

      // Send bottom row to next proc to be received as its new top border row
      // receive new bottom border row from the next proc (top row)
      // if end processor: skip (bc row below it is matrix border i.e. constant)
      if (rank != nProc - 1) {
        MPI_Send(bottomRow, nColumns,
                  MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD);

        MPI_Recv(bottomBorderRow, nColumns,
                  MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
    }
    else {
      // receive new bottom border row from the next proc (top row)
      // Send bottom row to next proc to be received as its new top border row
      // if end processor: skip (bc row below it is matrix border i.e. constant)
      if (rank != nProc - 1) {
        MPI_Recv(bottomBorderRow, nColumns,
                  MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Send(bottomRow, nColumns,
                  MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD);
      }

      // receive new top border row from the prev proc (bottom row)
      // Send top row to prev proc to be received as its new bottom border row
      // (skip if rank 0, row above is matrix border)
      if (rank != 0) {
        MPI_Recv(topBorderRow, nColumns,
                  MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Send(topRow, nColumns,
                  MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD);
      }
    }


    // We are doing a pointer swap, it does not matter what is inside the result
    // array i.e. newSubsection, as long as subsection is completely up to date
    // Pointer swap is relatively quick compared to copying all the values
    // manually.
    pointerSwap = newSubsection;
    newSubsection = subsection;
    subsection = pointerSwap;

    


  }
  

  // Combine all the results from each processor into one final matrix
  // GatherV allows us to collect all the data from each processor and store
  // them all in one resultant array. Here the data from each processor is the
  // calculation/edit zone.
  // The calibration of where the results will be stored is calculated in the
  // top of the program

  // The start of the processors edit zone (i.e. calculated values)
  float* processorEditZoneStart = newSubsection + nColumns;
  MPI_Gatherv(processorEditZoneStart, nElementsToEdit, MPI_FLOAT,
              matrix, recvElementsPerProc, recvDisplacement, MPI_FLOAT,
                                                      0, MPI_COMM_WORLD);

  // Get end time
  endTime = MPI_Wtime();
  if (rank == 0) {
    // Output format = Processors, Precision, Size, Time, Iterations,
    // This allows us to easily convert into CSV for easy data handling
    double timeTaken = endTime - startTime;
    printf("\n%d, %f, %d, %f, %d,",
        nProc, precision, nColumns, timeTaken, iterations);
    if (displayLevel > 1) {
      printf("Iterations: %d\n", iterations);
      printf("Time taken: %lf seconds\n\n", endTime - startTime);
      printf("---------------------------------------------\n\n");
      printf("\nFinal Matrix:\n");
      displayMatrix(matrix, nColumns, nColumns*nColumns);
      printf("\n");
    }
    fflush(stdout);
  }

  // Exiting MPI and freeing all the allocated memory
  // to return resources back to the cpu
  MPI_Finalize();
  if (rank == 0) {
    // These were only allocated on processor 0 so must only
    // be freed on processor 0
    free(matrix);
    free(sendElementsPerProc);
    free(recvElementsPerProc);
    free(sendDisplacement);
    free(recvDisplacement);
  }
  free(subsection);
  free(newSubsection);
  
}

// ------------------------------------------------------
// ------------------ Help functions --------------------
// ------------------------------------------------------


void generateMatrix(int size, float* matrix){
    for (int x = 0; x < size; x++) {
        matrix[x] = 1.0f;                  // top border
        matrix[x + (size-1)*size] = 1.0f;  // bottom border
        matrix[size * x] = 1.0f;            // left border
        matrix[size-1 + size * x] = 1.0f;  // right border
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
  int sendY = (int) floor(sendDisplacement / nColumns);
  int sendXEnd = (sendDisplacement + sendElementsPerProc) % nColumns;
  int sendYEnd = (int) floor(
                 (sendDisplacement + sendElementsPerProc) / nColumns);
  int recvX = recvDisplacement% nColumns;
  int recvY = (int) floor(recvDisplacement / nColumns);
  int recvXEnd = (recvDisplacement + recvElementsPerProc) % nColumns;
  int recvYEnd = (int) floor(
                 (recvDisplacement + recvElementsPerProc) / nColumns);


  printf("Processor %d send starts at (x,y) (%d,%d) "
         "and ends at (x,y) (%d,%d), rows: %d \n",
         p, sendX, sendY, sendXEnd, sendYEnd, sendYEnd-sendY );

  printf("Processor %d recv starts at (x,y) (%d,%d) "
         "and ends at (x,y) (%d,%d), rows: %d \n\n",
         p, recvX, recvY, recvXEnd, recvYEnd, recvYEnd-recvY );
}
//000000000000000000000000000000000000000000000000000000000000000000000000000000 - 80 columns
