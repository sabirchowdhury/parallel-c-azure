#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <pthread.h>
#include <sys/time.h>
#include <time.h>
#include <string.h>


pthread_barrier_t barrier;
pthread_mutex_t nPreciseThreadsLock;

typedef struct PositionXY {
    int x;
    int y;
} TPosition;
struct timeval startTime, endTime;
typedef struct Matrices {
    double ** matrix;
    double ** matrixNew;
} TMatrices;
typedef struct ArgOptions {
    double fillValue;
    bool hasFillValue;
    double borderValue;
    bool hasBorderValue;
    bool display;
    int seed;
    bool hasSeed;

} ArgOptions;
typedef struct RelaxMatrixArgs {
    int threadIndex;
    double precision;
    int nThreads;
    int size;
    int start;
    int end;;
    ArgOptions *argOptions;
    TMatrices *matrices;
    int * notPrecise;
} TRelaxMatrixArgs;


void multiThread(int nThreads, int size, double precision, ArgOptions argOptions, TMatrices *matrices);
void relaxMatrix();
void displayMatrix(int size, double ** matrix);
double calcAverage(int x, int y, double ** matrix);
TPosition getPositionXY(int position, int size);
TMatrices generateMatrices(int size, ArgOptions argOptions);


 


int main (int argc, char * argv[])
{
    // Set default arguments
    int size = 20;
    int nThreads = 1;
    double precision = 0.01;
    ArgOptions argOptions = {
        .borderValue = 0,
        .display = false,
        .fillValue = 0,
        .hasBorderValue = false,
        .hasFillValue = false,
        .hasSeed = false,
        .seed = 0,
    };


    // Process command line arguments
    for (int argIndex = 1; argIndex < argc - 1; argIndex++ ) {
        if (!strcmp(argv[argIndex], "-size")) {
            int sizeArg = atoi(argv[argIndex + 1]);
            size = sizeArg < 3 ? 3 : sizeArg;
        }
        if (!strcmp(argv[argIndex], "-threads")) {
            nThreads = atoi(argv[argIndex + 1]);
        }
        if (!strcmp(argv[argIndex], "-precision")) {
            double precisionArg = (double) atof(argv[argIndex + 1]);
            precision = precisionArg < 0 ? 0.01 : precisionArg;
        }
        if (!strcmp(argv[argIndex], "-display") || !strcmp(argv[argIndex + 1], "-display")) {
            argOptions.display = true;
        }
        if (!strcmp(argv[argIndex], "-border")) {
            argOptions.borderValue = atof(argv[argIndex + 1]);
            argOptions.hasBorderValue = true;
        }
        if (!strcmp(argv[argIndex], "-fill")) {
            argOptions.fillValue = atof(argv[argIndex + 1]);
            argOptions.hasFillValue = true;
        }
        if (!strcmp(argv[argIndex], "-seed")) {
            argOptions.seed = atoi(argv[argIndex + 1]);
            argOptions.hasSeed= true;

        }
    }

    // validate nThhrads
    int totalCells = (size - 2)*(size - 2);
    nThreads = nThreads > totalCells ? totalCells : (nThreads < 1 ? 1 : nThreads);

    // Start timer and initialise random seed and barrier
    gettimeofday(&startTime, NULL);
    srand((unsigned int) (argOptions.hasSeed ? argOptions.seed : time(0)));
    pthread_barrier_init(&barrier, NULL, (unsigned int) nThreads);

    // create matrices and start multithreading process
    TMatrices matrices = generateMatrices(size, argOptions);
    multiThread(nThreads, size, precision, argOptions, &matrices);

    // display option to print initial and final matrices
    if (argOptions.display) {
        displayMatrix(size, matrices.matrix);
    }

;

    // end timer and print results
    gettimeofday(&endTime, NULL);
    double timeTakenSeconds = (double) endTime.tv_sec - startTime.tv_sec;
    double timeTakenUSeconds = (double) ( endTime.tv_usec - startTime.tv_usec ) / 1000000  ;
    printf("%d, %d, %f, %f,\n", size, nThreads, precision, timeTakenUSeconds + timeTakenSeconds);

    // clean up memory allocations
    pthread_barrier_destroy(&barrier);
    for(int i = 0; i < size; i++){
        free(matrices.matrix[i]);
        free(matrices.matrixNew[i]);
    }
    free(matrices.matrix);
    free(matrices.matrixNew);

    return 0;
}

void multiThread (int nThreads, int size, double precision, ArgOptions argOptions, TMatrices *matrices) {

    // Individual memory and thread variables for each thread
    pthread_t threads[nThreads];
    TRelaxMatrixArgs relaxMatrixArgs[nThreads];

    // calculations to assist with sharing of cells
    int totalCells = (size - 2) * (size - 2);
    int cellsPerThread = (int) floor(totalCells / nThreads);
    int remainderCells = totalCells % nThreads;
    int extraCell = 0;
    int notPrecise = 1;

    int prevEnd = 0;

    for (int threadIndex = 0; threadIndex < nThreads; threadIndex++) {

        // initialise common varibales to be passed
        relaxMatrixArgs[threadIndex] = (TRelaxMatrixArgs) {
            .threadIndex = threadIndex,
            .matrices = matrices,
            .notPrecise = &notPrecise,
            .size = size,
            .precision = precision,
            .nThreads = nThreads,
            .argOptions = &argOptions,
         };

        // start cells at previous thread end
        relaxMatrixArgs[threadIndex].start = prevEnd;

        // checks if the thread should pick up an extra cell to distribute the remainder
        extraCell = (threadIndex<remainderCells) ? 1 : 0;
        relaxMatrixArgs[threadIndex].end = prevEnd + cellsPerThread + extraCell;
        prevEnd = relaxMatrixArgs[threadIndex].end;
        
        // create thread
        pthread_create(&(threads[threadIndex]), NULL, (void*(*)(void *)) relaxMatrix, (void *) &(relaxMatrixArgs[threadIndex]));


    }

    // wait for thread to end
    for (int threadPos = 0; threadPos < nThreads; threadPos++) {
        pthread_join(threads[threadPos], NULL);

    }
};



TMatrices generateMatrices(int size, ArgOptions argOptions) {

    // allocate matrix memory
    double * values = calloc((unsigned int) (size*size), sizeof(double));
    double * valuesNew = calloc((unsigned int) (size*size), sizeof(double));
    double ** matrix = malloc((unsigned int) size * sizeof(double *));
    double ** matrixNew = malloc((unsigned int) size * sizeof(double *));
    for (int x = 0; x < size; x++) {
      matrix[x] = values + x*size;
      matrixNew[x] = valuesNew + x*size;
    }

    // if there is a border value, fill all the borders with that value, else fill with random values
    for (int x = 0; x < size; x++) {
        matrix[x][0] = matrixNew[x][0] = argOptions.hasBorderValue ? argOptions.borderValue : (double) ( rand() % 1000) / 1000 ;
        matrix[x][size-1] = matrixNew[x][size-1] = argOptions.hasBorderValue ? argOptions.borderValue : (double) ( rand() % 1000) / 1000 ;
        matrix[0][x] = matrixNew[0][x] = argOptions.hasBorderValue ? argOptions.borderValue : (double) ( rand() % 1000) / 1000 ;
        matrix[size-1][x] = matrixNew[size-1][x] = argOptions.hasBorderValue ? argOptions.borderValue : (double) ( rand() % 1000) / 1000 ;

    }

    // return both matrices
    TMatrices matrices = {
        .matrix = matrix,
        .matrixNew = matrixNew,
    }; 
    return matrices;
};


void displayMatrix (int size, double ** matrix) {
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            printf ("%12.5f, ",matrix[i][j]);
        }
        printf ("\n");
    }
    
}


void relaxMatrix (TRelaxMatrixArgs *relaxMatrixArgs) {

    double ** matrix = relaxMatrixArgs->matrices->matrix;
    double ** matrixNew = relaxMatrixArgs->matrices->matrixNew;

    // check if there is a seed value for randomiser
    // seed the matrix (done in parallel)
    bool hasSeedArg = !!relaxMatrixArgs->argOptions->seed;
    int seedValue = hasSeedArg ? relaxMatrixArgs->argOptions->seed : time(0);
    srand((unsigned int) (seedValue & (relaxMatrixArgs->threadIndex+1)));
    double cellValue = 0;
    for (int position = relaxMatrixArgs->start; position < relaxMatrixArgs->end; position++) {
        TPosition xy = getPositionXY(position,relaxMatrixArgs->size);
        cellValue = (double) ( rand() % 1000) / 1000 ;
        matrix[xy.y][xy.x] = relaxMatrixArgs->argOptions->hasFillValue ? relaxMatrixArgs->argOptions->fillValue : cellValue;
        matrixNew[xy.y][xy.x] = 0;
    }


    int iterationCount = 0;


    // exit code: break loop if not precise

    while (*relaxMatrixArgs->notPrecise == 1) {



        // barrier to stop reset happening before a thread lagging behind reaches loop
        // can cause thread to exit if flag is reset before thread has looped
        pthread_barrier_wait(&barrier);

        // update pointer for global pointer of arrays to be consistent
        // reset precise flag
        // if statement prevents all threads doing this
        // only one therad needs to execute this command
        if (relaxMatrixArgs->threadIndex == 0) {
            *relaxMatrixArgs->notPrecise = 0;
            relaxMatrixArgs->matrices->matrix = matrix;
            relaxMatrixArgs->matrices->matrixNew = matrixNew;

            // display options
            if (relaxMatrixArgs->argOptions->display) {
                if (iterationCount == 0) {
                    printf("Initial matrix:\n\n");
                    displayMatrix(relaxMatrixArgs->size, matrix);
                    printf("\n\nFinal matrix:\n\n");


                }
            }
        }

        // barrier to ensure matrix pointer is updated with new values before all threads start executing
        // barrier also stops reset happening midway through execution
        pthread_barrier_wait(&barrier);


        // for loop sweeps through linear positions
        for (int position = relaxMatrixArgs->start; position < relaxMatrixArgs->end; position++) {
            // converts linear position to x and y coordinate
            TPosition xy = getPositionXY(position,relaxMatrixArgs->size);
            // calculate average
            matrixNew[xy.y][xy.x]=( matrix[xy.y][xy.x + 1] + matrix[xy.y + 1][xy.x] + matrix[xy.y - 1][xy.x] + matrix[xy.y][xy.x - 1] ) / 4;
            // set not precise flag if difference between curreent matrix and previous result is more than precision
            if (fabs(matrix[xy.y][xy.x] - matrixNew[xy.y][xy.x]) > relaxMatrixArgs->precision) {
                *relaxMatrixArgs->notPrecise = 1; // setting to 1 (constant) so do not need mutex
            }
        }

        // swap local matrix pointers
        double ** tempMatrixNewSwap = matrixNew;
        matrixNew = matrix;
        matrix = tempMatrixNewSwap;
        // iterate count increment
        iterationCount++;
        // barrier ensures that all threads have finished checking precision before thread decides to exit 
        pthread_barrier_wait(&barrier);
    }

    // Final update of global pointer if first thread
    // if statement prevents all threads doing this
    // only one therad needs to execute this command
    if (relaxMatrixArgs->threadIndex == 0) {
        printf("%d, ",iterationCount - 1);
        relaxMatrixArgs->matrices->matrix = matrix;
        relaxMatrixArgs->matrices->matrixNew = matrixNew;
    }

    pthread_exit(NULL);
};


TPosition getPositionXY(int position, int size) {
    // find x y position from linear position
    TPosition xy = {
        .y = (int) floor(position/(size - 2)) + 1,
        .x = position % (size - 2) + 1,
    };
    return xy;
}
