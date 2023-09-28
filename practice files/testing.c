#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>


void neighbourAverage (double image[4][4], double newImage[4][4], int i, int j);

double ** generate_matrix (int size);

int main (void)
{
  double image[4][4] = {
    {1.0, 2.0, 3.0, 4.0},
    {5.0, 0, 0, 0},
    {9.0, 0, 0, 0},
    {13.0, 0, 0, 0},
  };
  
  double newImage[4][4];
  
  // copy border
  
  for(int imageLen = 0; imageLen < 4; imageLen++) {
      newImage[imageLen][0] = image[imageLen][0];
      newImage[imageLen][3] = image[imageLen][3];
      newImage[0][imageLen] = image[0][imageLen];
      newImage[3][imageLen] = image[3][imageLen];
  }
  
  int precision = 0;
  

  while(precision < 3) {
      for(int i= 1; i < 3; i++) {
          for(int j = 1; j < 3; j++) {
              neighbourAverage(image, newImage, i, j);
          }
      }
      
      //copy new image to old image
      precision = 0;
      
      for(int i = 1; i < 3; i++) {
          for(int j = 1; j < 3; j++) {
              // printf ("%f, %f, %f\n",fabs(image[i][j] - newImage[i][j]), image[i][j], newImage[i][j]);

              if (fabs(image[i][j] - newImage[i][j]) < 0.0001) {
                  precision++;
              }
              image[i][j] = newImage[i][j];
          }
      }
      
  }
  
  
  printf ("\n\n");

  // replace these is and js

  int size = 9;

  double** matrix = generate_matrix(size);
  
  for(int i = 0; i < size; i++) {
      for(int j = 0; j < size; j++) {
          printf ("%f, ", matrix[i][j]);
      }
      printf ("\n");
  }
  //neighbourAverage(image);
  

  
//   pthread_create (&thr1, NULL, (void *(*)(void *)) neighbourAverage,
//           (void *) *image);
//   pthread_join (thr1, NULL);
  return 0;
}


double ** generate_matrix (int size) {
    // some sort of functionality yo check if its a posotive integer greater than 0

    double* values = calloc(size*size, sizeof(double));
    double ** matrix = malloc( size * sizeof(double *));

    for (int x = 0; x < size; x++) {
      matrix[x] = values + x*size;
    }

    for (int x = 0; x < size; x++) {
      matrix[x][0] = 1;
      matrix[x][size-1] = 1;
      matrix[0][x] = 1; 
      matrix[size-1][x] = 1;
    }

    return matrix;
}


void neighbourAverage (double image[4][4], double newImage[4][4], int i, int j)
{
    newImage[i][j] = ( image[i][j + 1] + image[i + 1][j] + image[i - 1][j] + image[i][j - 1] ) / 4;
}