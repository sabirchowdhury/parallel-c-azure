#include<stdio.h>
#include<unistd.h>
#include<stdlib.h>
#include<time.h>
#include <math.h>




int main(){

  int nProc = 176;

  int nColumns = 1000;
  int rowsPerProc = floor((nColumns - 2)/nProc);
  int nProcsWithExtraRow = (nColumns - 2) % nProc;


  

  int* sendElementsPerProc = malloc(sizeof(int)*nProc);
  int* recvElementsPerProc = malloc(sizeof(int)*nProc);

  int* sendDisplacement = malloc(sizeof(int)*nProc);
  int* recvDisplacement = malloc(sizeof(int)*nProc);

  for (int p = 0; p < nProc; p++) {
    // to divide up when rows do not perfectly divide the array size
    // figures if code should take/add an extra row to balance out workload
    int rowsAssignedToProcessorP = p < nProcsWithExtraRow ? (rowsPerProc + 1) : rowsPerProc;

    // Calculate row displacements
    // calculates the extra displacement from previous rows taking up extra rows (above /\)
    // if processor number is less than n processors that need an extra row
    // that means that ALL preceding processors have taken an extra row - no. of prev rows = "p"
    // ELSE if p doesnt need to take an extra row, that means all n processors that needed to take
    // an extra row have taken thier extra rows. So that means the total extra displacements is
    // the number of proccesors that needed an extra row i.e. "nProcsWithExtraRow" 
    int extraRowsFromPrevProcs = p < nProcsWithExtraRow ? p : nProcsWithExtraRow;
    int rowDisplacement = p * rowsPerProc + extraRowsFromPrevProcs + 1; // + 1 to displace border row


    sendElementsPerProc[p] = nColumns * (rowsAssignedToProcessorP + 2); // + 2 for upper and lower rows
    recvElementsPerProc[p] = nColumns * rowsAssignedToProcessorP;

    sendDisplacement[p] = nColumns * (rowDisplacement - 1); // - 1 to include the row above (for calculation)
    recvDisplacement[p] = nColumns * rowDisplacement;

    // the send uses two extra rows as oppose to the recieve
    // these extra rows are required for the above and below averaging in relaxation


    // Unit testing
    int sendX = sendDisplacement[p] % nColumns;
    int sendY = floor(sendDisplacement[p] / nColumns);
    int sendXEnd = (sendDisplacement[p] + sendElementsPerProc[p]) % nColumns;
    int sendYEnd = floor((sendDisplacement[p] + sendElementsPerProc[p]) / nColumns);
    int recvX = recvDisplacement[p] % nColumns;
    int recvY = floor(recvDisplacement[p] / nColumns);
    int recvXEnd = (recvDisplacement[p] + recvElementsPerProc[p]) % nColumns;
    int recvYEnd = floor((recvDisplacement[p] + recvElementsPerProc[p]) / nColumns);
    printf("Processor %d send starts at (x,y) (%d,%d) and ends at (x,y) (%d,%d), rows: %d \n", p, sendX, sendY, sendXEnd, sendYEnd, sendYEnd-sendY );
    printf("Processor %d recv starts at (x,y) (%d,%d) and ends at (x,y) (%d,%d), rows: %d \n\n", p, recvX, recvY, recvXEnd, recvYEnd, recvYEnd-recvY );

  }



}

