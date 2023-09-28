// Include libraries

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
// OPENMPI library responsible of distributed memory architecture
#include <mpi.h>
#include <omp.h>
#include <math.h>
#include <sys/time.h>

// ______________________________________Initialisation of variables_____________________________________
// initialisation of functions used
struct matrices create_matrix(int arraySize); 
void print_matrix(double* matrix, int arraySize);
void relaxTechnique(double* buffer_recv, double* buffer_send, int arraySize, int currentReceive);

// Structure to hold the two Matrices
typedef struct matrices {
    double *matrix1;  // Matrix for reading
    double *matrix2; // Matrix for writing
} matrices;

// Struct for timing variables intialisation
double start,end;

// variables for testing all numbers
int total;

// ______________________________________Main_Program_____________________________________
// Main function
int main(int argc, char *argv[]) {
    
    // Number of cores intialisation will be determined by the system with the OPENMPI Function
    int nCores;

    // MPI variables intialisation
    int errorMPI; // This variable hold if an error occurs with the OPENMPI function
    int ranks; // Current processor being used

    // MPI basics initialisation
    printf("1111");
    errorMPI = MPI_Init(&argc, &argv); 
    printf("2222");

	errorMPI = MPI_Comm_rank(MPI_COMM_WORLD, &ranks); // Communication of ranks
        printf("3333");

	errorMPI = MPI_Comm_size(MPI_COMM_WORLD, &nCores);// Communication of number of cores
    printf("4444");


    // Main variables intialisation
    int overallPrecision = 0; // Hold the precision of the system default 0 precision is not reached
    int arraySize; // matrix of ArraySize * ArraySize
    int testing; // Unit Testing
    double precision; // Hold the required precision

    // This condition check for given arguments in command line
    if (argc < 2) {
        // In case nothing given provide default values
      arraySize = 5000;
      precision = 0.1;
      testing = 0;
    }
    else { // Required arguments run for those values
      arraySize = atoi(argv[1]);
      precision = atof(argv[2]);
      testing = atoi(argv[3]);
    }

    //  Creating the two matrix we want to analyse
    matrices Matrices;
    Matrices = create_matrix(arraySize); //Create the 1D matrices 
   
    // Printing the initial matrices
    if (ranks == 0) { // Only using the main core to print initial matrix
        start = MPI_Wtime(); // Start the timing of the program done by the main core
        
        if (testing == 1) {
            printf("\nWe are using %i \n",nCores); // Check the number of Cores used for this analysis useful for testing
            printf(" \n \t Initial Matrices\n");
            print_matrix(Matrices.matrix1,arraySize);
            print_matrix(Matrices.matrix2,arraySize);
            printf(" \n __________________________________________________________________________________________________________________\n");
        }
    }

  // Each core has an equal number of rows to analyse
    int rowsPerCore = ( arraySize-2)/nCores; // we removed two because we do not want the edges values to be analysed

    // // These hold the communication arrays between the cores to be given to the broadcastin function of OPENMPI
    // Values to be send
    int *startingPos_Send = malloc(nCores*sizeof(int)); // Starting position to be sent or displacement in documentation
	int *nValues_Send     = malloc(nCores*sizeof(int)); // Number of values to be sent or send_count in documentation

    // Values to be received
	int *startingPos_Received = malloc(nCores*sizeof(int)); // Starting position to be received or gather_displacement in documentation
	int *nValues_Received   = malloc(nCores*sizeof(int)); // Number of values  to be received or recv_count in documentation

    //current size to receive and send 
    int currentSend; //Current buffer sent opposite of the send above
    int currentReceive; //Current buffer received opposite of the received above

    // We want to divide the numbers to each core 
    // If the number of cells divided by the number of cores we are using is even give them each same number
    // If the number of cells divided by the number of cores we are using is odd give the last one the rest of the cells

    // Determining sizing 
    int sizeReceive = (rowsPerCore * arraySize); // Smaller size to receive
    int sizeSend = ((rowsPerCore + 2) * arraySize);// Larger size to receive for neighboring values
    int sizeReceive_LastCore = ((rowsPerCore + ((arraySize-2) % nCores)) * arraySize); // Giving reminders to the last core Received
    int sizeSend_LastCore = ((rowsPerCore + ((arraySize-2) % nCores) + 2) * arraySize);// Giving reminders to the last core Send
    
    // Number of elements current buffer( check for even or odd as explained above) to be done by all cores

        if (ranks != (nCores-1)) {
            currentSend = sizeReceive ; // Not last core sending is the amount they received(maincore)
            currentReceive = sizeSend ;// Not last core receiving is the amount they received(maincore)
            
        }
        else {
            currentSend = sizeReceive_LastCore; //  Last core sending is the amount they received(maincore)
            currentReceive = sizeSend_LastCore;//  Last core receiving is the amount they sending(maincore
        }

        // Buffer arrays to be given as documentation hilighted in scatterv
        double *buffer_recv = malloc(currentReceive * sizeof(double)); // Use documentation buffer array receive
	    double *buffer_send = malloc(currentSend * sizeof(double));



    // Looping through all the cores to determine the data intialised above
    if (ranks == 0) { // Only main core perform this operation

        for(int i = 0; i < nCores;i++){


            // Number of elements to give ( check for even or odd as explained above)
            if ( i != (nCores-1)){ // If is the last core to check
                 nValues_Send[i] = sizeSend; // Gives send numbers for all other cores 
                nValues_Received[i] = sizeReceive;// Gives received numbers for all other cores 
            
            }
            else {
                nValues_Send[i] = sizeSend_LastCore; // Give reminders send + rows calculated above
                nValues_Received[i] = sizeReceive_LastCore; // Give reminders received + rows calculated above
               
            }

            // Starting position same process for all Core no need to check for last core as it doesnt change
            startingPos_Send[i] = i * rowsPerCore * arraySize;
            startingPos_Received[i] = arraySize + (i * rowsPerCore * arraySize);

            if (testing == 1){
                // Unit testing total number of cells
                total = total + nValues_Received[i];
            }

            }
        }
        else {
            //  If its not main core broadcast NULL
            nValues_Send = NULL;
            startingPos_Send = NULL;
        }

    
     if (ranks == 0){ // Testing done by main core only
        if (testing == 1){
            printf("\n Buffer current send: %i",currentSend);
            printf("\n Buffer current received: %i",currentReceive);
            printf("\n nValuesSend: %i",nValues_Send[3]);
            printf("\n nValuesReceived: %i",nValues_Received[3]);
            printf("\n Starting pos send: %i",startingPos_Send[3]);
            printf("\n Starting pos received: %i",startingPos_Received[3]);

            printf("Total Number of cells: %i",total);
            
        }   
     }

    // Infine loop break when precision has been reached
    while (1){

        // Using scatterv because not all core have the same amount of data
        MPI_Scatterv(Matrices.matrix1, nValues_Send, startingPos_Send, MPI_DOUBLE, buffer_recv, currentReceive, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Relax technique done by all the cores
        relaxTechnique(buffer_recv, buffer_send, arraySize, currentReceive);

        // Gather function 
        MPI_Gatherv(buffer_send, currentSend, MPI_DOUBLE, Matrices.matrix2, nValues_Received,startingPos_Received, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Calculation of overall precision
        if (ranks == 0){

            // Temporary overallPrecision checker its the opposite of overallPrecision
			int test = 0;
            // Looping through all the 1D array to determine if precision has been reached by comparing the matrices
            for (int i = 0; i < arraySize*arraySize; i++){
		        if (fabs(Matrices.matrix1[i] - Matrices.matrix2[i]) > precision){
			        test = 1;
		        }
	        }

            // Updating the overallPrecision  to be broadcast for other core to know when to stop
            if (test == 1){
                overallPrecision = 0;
            }
            else {
                overallPrecision = 1;
            }
		}

        // Letting all core know the result of the threshold
        MPI_Bcast(&overallPrecision, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // As in cW1 swapping the matrices using a temp matrix
        double *holdingMatrix = Matrices.matrix1;
        Matrices.matrix1 = Matrices.matrix2;
        Matrices.matrix2 = holdingMatrix;

        if (overallPrecision){
			break; // Here all the core will be break from the loop as we broadcasted the overallPrecision pointer to all available cores
		}
   
    }

    //  Done program close the MPI
    errorMPI = MPI_Finalize();

    if (ranks == 0) { //Done by the main core
        // Determining timing of processing time 
        end = MPI_Wtime();
        printf("Timing: %f seconds\n",end - start ); // Printing the timing it took outside of testing area

        if (testing == 1) {
            printf(" \n \t Final Matrices\n");
            print_matrix(Matrices.matrix1,arraySize);
            print_matrix(Matrices.matrix2,arraySize);
            printf(" \n __________________________________________________________________________________________________________________\n");
        }

    }
}

// Relaxation Technique as coursework
void relaxTechnique(double* buffer_recv, double* buffer_send, int arraySize, int currentReceive){

    int initPosition = arraySize; 

    int endPosition = currentReceive - arraySize - 1 ;

    int j = 0;

    // Analysing the values of chunk of data sent to each cores
    for (int i = initPosition; i <= endPosition; i++){

        //  Removing edges values
        if (  ((i+1) % arraySize != 0) && (i % arraySize != 0) ){

            // Neighboring average
            double averageCalc =  (buffer_recv[i+1]+buffer_recv[i+arraySize]+buffer_recv[i-1]+buffer_recv[i-arraySize])/4;
            // Assign it
            buffer_send[j] = averageCalc;

        }
        else {
            // Edges values keep it the same = 1
            buffer_send[j] = buffer_recv[i];
        }
        // Increment counter
        j++;
    }
}   
// ______________________________________Function_For_Testing_____________________________________



//Function to create a 2D matrix
struct matrices create_matrix(int arraySize) {
    // We have now a pointer to array of pointers 
    matrices Matrices;
    double *matrix1 = malloc(arraySize * arraySize * sizeof(double*));
    double *matrix2 = malloc(arraySize * arraySize * sizeof(double*));

    // Failed malloc function
    if (matrix1 == NULL || matrix2 == NULL){
        printf("Error: Failed with malloc");
        exit(1);
    }

    // We are creating a 1D array as an illusion for a 2D array where each row has arraySize means arraySize*2 = 2ndRow
    for(int i = 0; i < arraySize*arraySize; i++) {

        // Filling up the data
            if (i < arraySize) {// First row
                matrix1[i] = 1;
                matrix2[i] = 1;
            }
            else if (((i+1) % arraySize) == 0) {// First column
                matrix1[i] = 1;
                matrix2[i] = 1;
            }
            else if ((i % arraySize) == 0) { // Last column
                matrix1[i] = 1;
                matrix2[i]= 1;
            }
            
            else if ((i < arraySize * arraySize) && (i >= arraySize * (arraySize-1))) {// First column
                matrix1[i] = 1;
                matrix2[i] = 1;
            }
            else {
                matrix1[i] = 0.0005;
                matrix2[i] = 0;
            }
    }

    Matrices.matrix1 = matrix1;
    Matrices.matrix2 = matrix2;
    
    return Matrices;

}


// Function to print the 2D matrix based on a 1D storage
void print_matrix(double* matrix, int arraySize) {

    // printing the matrix
    printf("\n______________2D_Matrix______________\n\n");
    for(int i = 0; i < arraySize*arraySize; i++) {
        if (i % arraySize == 0){
            printf("\n");
        }
        printf("%f ", matrix[i]);
    }
    printf("\n");
}
