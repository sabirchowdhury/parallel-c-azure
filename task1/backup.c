#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <pthread.h>
#include <time.h>
#include <string.h>

// TODO - catch the failed memory allocs
// TODO - catch the errors from the thread create

pthread_barrier_t barrier;
pthread_mutex_t nPreciseThreadsLock;

typedef struct PositionXY {
    int x;
    int y;
} TPosition;
typedef struct Matricies {
    double ** matrix;
    double ** matrixNew;
} TMatricies;
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
    TMatricies *matricies;
    int * nPreciseThreads;
} TRelaxMatrixArgs;



void multiThread(int nThreads, int size, double precision, ArgOptions argOptions, TMatricies *matricies);
void relaxMatrix();
void displayMatrix(int size, double ** matrix);
double calcAverage(int x, int y, double ** matrix);
TPosition getPositionXY(int position, int size);
TMatricies generateMatricies(int size, ArgOptions argOptions);


 


int main (int argc, char * argv[])
{

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
    
    for (int argIndex = 1; argIndex < argc - 1; argIndex++ ) {
        if (!strcmp(argv[argIndex], "-size")) {
            int sizeArg = atoi(argv[argIndex + 1]);
            size = sizeArg < 3 ? 3 : sizeArg;
        }
        if (!strcmp(argv[argIndex], "-threads")) {
            nThreads = atoi(argv[argIndex + 1]);
        }
        if (!strcmp(argv[argIndex], "-precision")) {
            double precisionArg = (double) atoi(argv[argIndex + 1]);
            precision = precisionArg < 0 ? 0.01 : precisionArg;
        }
        if (!strcmp(argv[argIndex], "-display")) {
            int displayArg = atoi(argv[argIndex + 1]);
            argOptions.display = displayArg == 1 ? true : false;
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

    int totalCells = (size - 2)*(size - 2);
    srand((unsigned int) (argOptions.hasSeed ? argOptions.seed : time(0)));
    nThreads = nThreads > totalCells ? totalCells : (nThreads < 1 ? 1 : nThreads);












    // if (argc>1) {
    //     int sizeArg = atoi(argv[1]);
    //     size = sizeArg < 3 ? size : sizeArg;
        
    // }
    // if (argc>2) {
    //     int nThreadsArg = atoi(argv[2]);
    //     int totalCells = (size - 2)*(size - 2);
    //     nThreads = nThreadsArg > totalCells ? totalCells : (nThreads < 1 ? 1 : nThreads);
    // }
    // if (argc>3) {
    //     double precisionArg = (double) atoi(argv[3]);
    //     precision = precisionArg < 0 ? 0.01 : precisionArg;
    // }
    // if (argc>4) {
    //     int displayArg = atoi(argv[4]);
    //     display = displayArg == 1 ? true : false;
    // }


    // switch (argc) {
    //     case 1:
    //         break;
    //     case 2:
    //         printf("%s",argv[1]);
    //         break;
    //     case 3:
    //         printf("%s",argv[1]);
    //         printf("%s",argv[2]);

    //         break;
    //     case 4:
    //         printf("%s",argv[1]);
    //         printf("%s",argv[2]);
    //         printf("%s",argv[3]);
    //         break;

    //     default:
    //         if (argc > 4) {
    //         printf("%s",argv[1]);
    //         printf("%s",argv[2]);
    //         printf("%s",argv[3]);
    //         printf("%s",argv[4]);
    //         }
    //         break;
    // }

    // exit(0);
    


    pthread_mutex_init(&nPreciseThreadsLock, NULL);
    pthread_barrier_init(&barrier, NULL, (unsigned int) nThreads); // check if you need to initilise it with one higher

    TMatricies matricies = generateMatricies(size, argOptions);

    multiThread(nThreads, size, precision, argOptions, &matricies);

    if (argOptions.display) {
        displayMatrix(size, matricies.matrix);
    }

    pthread_mutex_destroy(&nPreciseThreadsLock);
    pthread_barrier_destroy(&barrier);

    return 0;
}


void multiThread (int nThreads, int size, double precision, ArgOptions argOptions, TMatricies *matricies) {
    pthread_t threads[nThreads];
    TRelaxMatrixArgs relaxMatrixArgs[nThreads];

    int totalCells = (size - 2) * (size - 2);
    int cellsPerThread = (int) floor(totalCells / nThreads);
    int remainderCells = totalCells % nThreads;
    int extraCell = 0;
    int nPreciseThreads = 0;
    // **printf("cellsoerthred %d %d %d\n",cellsPerThread, remainderCells, totalCells);

    int prevEnd = 0;


    // todo: if more threads than cores do not create extras 
    for (int threadIndex = 0; threadIndex < nThreads; threadIndex++) {
        relaxMatrixArgs[threadIndex] = (TRelaxMatrixArgs) {
            .threadIndex = threadIndex,
            .matricies = matricies,
            .nPreciseThreads = &nPreciseThreads,
            .size = size,
            .precision = precision,
            .nThreads = nThreads,
            .argOptions = &argOptions,
         };

        relaxMatrixArgs[threadIndex].start = prevEnd;
        extraCell = (threadIndex<remainderCells) ? 1 : 0;
        relaxMatrixArgs[threadIndex].end = prevEnd + cellsPerThread + extraCell;
        prevEnd = relaxMatrixArgs[threadIndex].end;
        
        pthread_create(&(threads[threadIndex]), NULL, (void*(*)(void *)) relaxMatrix, (void *) &(relaxMatrixArgs[threadIndex]));


    }


    for (int threadPos = 0; threadPos < nThreads; threadPos++) {
        pthread_join(threads[threadPos], NULL);

    }
};



TMatricies generateMatricies(int size, ArgOptions argOptions) {
    double * values = calloc((unsigned int) (size*size), sizeof(double));
    double * valuesNew = calloc((unsigned int) (size*size), sizeof(double));
    double ** matrix = malloc((unsigned int) size * sizeof(double *));
    double ** matrixNew = malloc((unsigned int) size * sizeof(double *));
    for (int x = 0; x < size; x++) {
      matrix[x] = values + x*size;
      matrixNew[x] = valuesNew + x*size;
    }

    // Fill in the borders with 1s
    for (int x = 0; x < size; x++) {

        // Switch all these from rand to 1 and say that the matrix should be symmetric


        matrix[x][0] = matrixNew[x][0] = argOptions.hasBorderValue ? argOptions.borderValue : (double) ( rand() % 1000) / 100 ;
        matrix[x][size-1] = matrixNew[x][size-1] = argOptions.hasBorderValue ? argOptions.borderValue : (double) ( rand() % 1000) / 100 ;
        matrix[0][x] = matrixNew[0][x] = argOptions.hasBorderValue ? argOptions.borderValue : (double) ( rand() % 1000) / 100 ;
        matrix[size-1][x] = matrixNew[size-1][x] = argOptions.hasBorderValue ? argOptions.borderValue : (double) ( rand() % 1000) / 100 ;

    }

    TMatricies matricies = {
        .matrix = matrix,
        .matrixNew = matrixNew,
    }; 
    return matricies;
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
    // printf ("start of relax thread: %d\n",relaxMatrixArgs->threadIndex);
    bool hasSeedArg = !!relaxMatrixArgs->argOptions->seed;
    int seedValue = hasSeedArg ? relaxMatrixArgs->argOptions->seed : time(0);

    srand((unsigned int) (seedValue & (relaxMatrixArgs->threadIndex+1)));

    int iterationCount = 0;


    double ** matrix = relaxMatrixArgs->matricies->matrix;
    double ** matrixNew = relaxMatrixArgs->matricies->matrixNew;

    int totalPreciseCells = 0;
    bool matrixSeeded = false;


    while (*relaxMatrixArgs->nPreciseThreads < relaxMatrixArgs->nThreads) {
        // **printf ("start of loop thread: %d\n",relaxMatrixArgs->threadIndex);
        // **printf("%d %d\n",totalPreciseCells,relaxMatrixArgs->nThreads);
        // **printf("init: matrix ptr: %p matrixNew ptr %p\n",*matrix,*matrixNew);

        pthread_barrier_wait(&barrier);
        if (relaxMatrixArgs->threadIndex == 0) {
            *relaxMatrixArgs->nPreciseThreads = 0;
        }
        totalPreciseCells = 0;
    
        for (int position = relaxMatrixArgs->start; position < relaxMatrixArgs->end; position++) {
            TPosition xy = getPositionXY(position,relaxMatrixArgs->size);

            if (matrixSeeded == false) {

                double cellValue = (double) ( rand() % 1000) / 100 ;
                matrix[xy.y][xy.x] = relaxMatrixArgs->argOptions->hasFillValue ? relaxMatrixArgs->argOptions->fillValue : cellValue;//cellValue; // todo check if its okay to reference a pointer array like this
                matrixNew[xy.y][xy.x] = 0; // TODO - check if this assignment works


            } else {

                // printf("Before: Matrix New: %f Matirx Old: %f\n",matrixNew[xy.y][xy.x], matrix[xy.y][xy.x]);

                matrixNew[xy.y][xy.x]=calcAverage(xy.x,xy.y,matrix);
                // printf("Matrix New: %f Matirx Old: %f\n",matrixNew[xy.y][xy.x], matrix[xy.y][xy.x]);

                // printf("diff: %f precision: %f \n", fabs(matrix[xy.y][xy.x] - matrixNew[xy.y][xy.x]), relaxMatrixArgs->precision);
                if (fabs(matrix[xy.y][xy.x] - matrixNew[xy.y][xy.x]) < relaxMatrixArgs->precision) {
                    totalPreciseCells++;
                    // printf("mmmmm  cell no: %d  preciseCells %d  mmmm\n\n", position-relaxMatrixArgs->start, totalPreciseCells);
                }
            }
            // printf("to next %d\n", position);
        }
        // **printf("I am here 2, thread: %d\n",relaxMatrixArgs->threadIndex);

        pthread_barrier_wait(&barrier);

        pthread_mutex_lock(&nPreciseThreadsLock);
        // **printf("I am here 3, thread: %d\n",relaxMatrixArgs->threadIndex);

        if (totalPreciseCells >= relaxMatrixArgs->end - relaxMatrixArgs->start) {
            (*relaxMatrixArgs->nPreciseThreads)++;
        }
        pthread_mutex_unlock(&nPreciseThreadsLock);
        // **printf("I am here 4, thread: %d\n",relaxMatrixArgs->threadIndex);

        if(matrixSeeded == true) {
            double ** tempMatrixNewSwap = matrixNew;
            matrixNew = matrix;
            matrix = tempMatrixNewSwap; 
        }

        pthread_barrier_wait(&barrier);
        // **printf("I am here 5, thread: %d\n",relaxMatrixArgs->threadIndex);


        // if there is a bug its probably with this part below

        if (relaxMatrixArgs->threadIndex == 0) {
             
                // printf("wwwwwwwwwwwwwwwwwwwwwwwwwwwwww\n");
                // printf("%p %p dis \n",relaxMatrixArgs->matricies->matrixNew, relaxMatrixArgs->matricies->matrix );
                // printf("%d %d dos \n",relaxMatrixArgs->matricies->matrixNew[1][1], relaxMatrixArgs->matricies->matrix[1][1] );
                // printf("%p %p dat \n", matrixNew, matrix );

                relaxMatrixArgs->matricies->matrix = matrix;
                relaxMatrixArgs->matricies->matrixNew = matrixNew;
                iterationCount++;

                // printf("%p %p dis afta \n",relaxMatrixArgs->matricies->matrixNew, relaxMatrixArgs->matricies->matrix );
                // printf("%d %d dos afta  \n",relaxMatrixArgs->matricies->matrixNew[1][1], relaxMatrixArgs->matricies->matrix[1][1] );
                // printf("%p %p dat afta \n", matrixNew, matrix );
            
            if (matrixSeeded == false && relaxMatrixArgs->argOptions->display) { // if testing

                printf("begin: %f \n", matrix[1][1]);
                displayMatrix(relaxMatrixArgs->size, matrix);
                printf("\n\nStart^\n\nEND:\n\n");

            }
        }
        matrixSeeded = true;

        // pthread_mutex_lock(&nPreciseThreadsLock);
        // printf("pointer matrixNew: %p, matrixNew: %p pointer matrix: %p, matrix %p  , thread %d:: \n",relaxMatrixArgs->matricies->matrixNew,matrixNew ,relaxMatrixArgs->matricies->matrix, matrix, relaxMatrixArgs->threadIndex );
        // int input, i;

        // i = scanf("%d", &input);
        // pthread_mutex_unlock(&nPreciseThreadsLock);



        // **printf("end of loop, thread: %d\n",relaxMatrixArgs->threadIndex);
        }
    // **printf("pthread exit, thread: %d\n",relaxMatrixArgs->threadIndex);

    if (relaxMatrixArgs->threadIndex==0 && relaxMatrixArgs->argOptions->display) {
        printf("\n--------------------------------------------------\n\nNumber of iterations: %d\n\n",iterationCount);
    }
    pthread_exit(NULL);
};

double calcAverage (int x, int y, double ** matrix) {
    // if(x==1&&y==1){
    // printf ("calc av: %f + %f + %f + %f = %f\n",matrix[y][x + 1], matrix[y + 1][x], matrix[y - 1][x], matrix[y][x - 1] , (matrix[y][x + 1] + matrix[y + 1][x] + matrix[y - 1][x] + matrix[y][x - 1] ) / 4);
    // }
    return ( matrix[y][x + 1] + matrix[y + 1][x] + matrix[y - 1][x] + matrix[y][x - 1] ) / 4;
};

TPosition getPositionXY(int position, int size) {
    // fix the border workouts
    TPosition xy = {
        .y = (int) floor(position/(size - 2)) + 1,
        .x = position % (size - 2) + 1,
    };
    return xy;
}

// seed matrix run
// while could check the precision met
// create two matricies initially copy the borders in the main thread.
// do the relaxation
// lock the total precision
// update total precision + 1 for this thread if all are precise
// unlock total precision
// barrier
// if first thread: swap the pointers
// barrier
// check if the precision has been met - if yes thead exit, if no, total precise threads = 0
// loop
// thread exit/kill