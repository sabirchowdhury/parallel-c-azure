#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <pthread.h>
#include <sys/time.h>
#include <time.h>
#include <string.h>
//
// TODO - catch the failed memory allocs
// TODO - catch the errors from the thread create

pthread_barrier_t barrier;
pthread_mutex_t nPreciseThreadsLock;

typedef struct PositionXY {
    int x;
    int y;
} TPosition;
struct timeval startTime, endTime, TIS, TIE, TS, TE;
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
    int * notPrecise;
} TRelaxMatrixArgs;



double totalInner = 0;

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


    int totalCells = (size - 2)*(size - 2);
    srand((unsigned int) (argOptions.hasSeed ? argOptions.seed : time(0)));
    nThreads = nThreads > totalCells ? totalCells : (nThreads < 1 ? 1 : nThreads);

    gettimeofday(&startTime, NULL);

    pthread_barrier_init(&barrier, NULL, (unsigned int) nThreads);

    TMatricies matricies = generateMatricies(size, argOptions);

    multiThread(nThreads, size, precision, argOptions, &matricies);

    if (argOptions.display) {
        displayMatrix(size, matricies.matrix);
    }

    pthread_barrier_destroy(&barrier);

    gettimeofday(&endTime, NULL);

    double timeTakenSeconds = (double) endTime.tv_sec - startTime.tv_sec;
    double timeTakenUSeconds = (double) ( endTime.tv_usec - startTime.tv_usec ) / 1000000  ;
    double Tis = (double) TIE.tv_sec - TIS.tv_sec;
    double Tius = (double) ( TIE.tv_usec - TIS.tv_usec ) / 1000000  ;
    //printf("tIME INSIDE %f,\n", Tis + Tius );
    double Ts = (double) TE.tv_sec - TS.tv_sec;
    double Tus = (double) ( TE.tv_usec - TS.tv_usec ) / 1000000  ;
    //printf("TIME OUTSIDE %f,\n", Ts + Tus );
    //printf("total inner %f\n",totalInner);
    printf("%d, %d, %f, %f,\n", size, nThreads, precision, timeTakenUSeconds + timeTakenSeconds);
    free(matricies.matrix);
    free(matricies.matrixNew);


    return 0;
}

void multiThread (int nThreads, int size, double precision, ArgOptions argOptions, TMatricies *matricies) {
    pthread_t threads[nThreads];
    TRelaxMatrixArgs relaxMatrixArgs[nThreads];

    int totalCells = (size - 2) * (size - 2);
    int cellsPerThread = (int) floor(totalCells / nThreads);
    int remainderCells = totalCells % nThreads;
    int extraCell = 0;
    int notPrecise = 1;

    int prevEnd = 0;

    for (int threadIndex = 0; threadIndex < nThreads; threadIndex++) {
        relaxMatrixArgs[threadIndex] = (TRelaxMatrixArgs) {
            .threadIndex = threadIndex,
            .matricies = matricies,
            .notPrecise = &notPrecise,
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

    for (int x = 0; x < size; x++) {

        matrix[x][0] = matrixNew[x][0] = argOptions.hasBorderValue ? argOptions.borderValue : (double) ( rand() % 1000) / 1000 ;
        matrix[x][size-1] = matrixNew[x][size-1] = argOptions.hasBorderValue ? argOptions.borderValue : (double) ( rand() % 1000) / 1000 ;
        matrix[0][x] = matrixNew[0][x] = argOptions.hasBorderValue ? argOptions.borderValue : (double) ( rand() % 1000) / 1000 ;
        matrix[size-1][x] = matrixNew[size-1][x] = argOptions.hasBorderValue ? argOptions.borderValue : (double) ( rand() % 1000) / 1000 ;

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
    if (relaxMatrixArgs->threadIndex==1) {
            gettimeofday(&TS, NULL);
    }
    bool hasSeedArg = !!relaxMatrixArgs->argOptions->seed;
    int seedValue = hasSeedArg ? relaxMatrixArgs->argOptions->seed : time(0);
    srand((unsigned int) (seedValue & (relaxMatrixArgs->threadIndex+1)));
    int iterationCount = 0;


    double ** matrix = relaxMatrixArgs->matricies->matrix;
    double ** matrixNew = relaxMatrixArgs->matricies->matrixNew;

    bool matrixSeeded = false;
    struct timeval LE, LO;


    while (true) {

        //////////////////////////////////////////////////////////////////////
        // if (relaxMatrixArgs->threadIndex==1) {
        //     gettimeofday(&LE, NULL);
        //     gettimeofday(&TIS, NULL); // measurer /////////////

        // }
        ////////////////////////////////////////////////////////////////////////////

    
        pthread_barrier_wait(&barrier);
        if (*relaxMatrixArgs->notPrecise == 0) {
            break;
        }
        pthread_barrier_wait(&barrier);


        if (relaxMatrixArgs->threadIndex == 0) {
            if (matrixSeeded == true) {
                *relaxMatrixArgs->notPrecise = 0;
                relaxMatrixArgs->matricies->matrix = matrix;
                relaxMatrixArgs->matricies->matrixNew = matrixNew;
            }
            if (relaxMatrixArgs->argOptions->display && matrixSeeded == false) {
                printf("Initial matrix:\n\n");
                displayMatrix(relaxMatrixArgs->size, matrix);

            }
        }
        pthread_barrier_wait(&barrier);

        //////////////////////////////////////////////////////////////////////////////////
        // if (relaxMatrixArgs->threadIndex==1) {// End measure //////////////
        //     gettimeofday(&LO, NULL);
        //     gettimeofday(&TIE, NULL);
        //     double Tis = (double) LO.tv_sec - LE.tv_sec;
        //     double Tius = (double) ( LO.tv_usec - LE.tv_usec ) / 1000000  ;
        //     totalInner += Tis + Tius;
        // }
        /////////////////////////////////////////////////////////////////////////////////

    
        for (int position = relaxMatrixArgs->start; position < relaxMatrixArgs->end; position++) {
            TPosition xy = getPositionXY(position,relaxMatrixArgs->size);

            if (matrixSeeded == false) {

                double cellValue = (double) ( rand() % 1000) / 1000 ;
                matrix[xy.y][xy.x] = relaxMatrixArgs->argOptions->hasFillValue ? relaxMatrixArgs->argOptions->fillValue : cellValue;
                matrixNew[xy.y][xy.x] = 0;
                *relaxMatrixArgs->notPrecise = 1;


            } else {

                matrixNew[xy.y][xy.x]=calcAverage(xy.x,xy.y,matrix);
                if (fabs(matrix[xy.y][xy.x] - matrixNew[xy.y][xy.x]) > relaxMatrixArgs->precision) {
                    *relaxMatrixArgs->notPrecise = 1;
                }

            }
        }

        if(matrixSeeded == true) {
            double ** tempMatrixNewSwap = matrixNew;
            matrixNew = matrix;
            matrix = tempMatrixNewSwap;

        }
        matrixSeeded = true;
        iterationCount++;
    }

    if (relaxMatrixArgs->threadIndex == 0) {
        relaxMatrixArgs->matricies->matrix = matrix;
        relaxMatrixArgs->matricies->matrixNew = matrixNew;
    }

    if (relaxMatrixArgs->threadIndex==0 && (relaxMatrixArgs->argOptions->display || true)) {
        printf("%d, ",iterationCount - 1);

    }

    if (relaxMatrixArgs->threadIndex==1) {
        gettimeofday(&TE, NULL);
    }

    pthread_exit(NULL);
};

double calcAverage (int x, int y, double ** matrix) {
    return ( matrix[y][x + 1] + matrix[y + 1][x] + matrix[y - 1][x] + matrix[y][x - 1] ) / 4;
};

TPosition getPositionXY(int position, int size) {
    TPosition xy = {
        .y = (int) floor(position/(size - 2)) + 1,
        .x = position % (size - 2) + 1,
    };
    return xy;
}
