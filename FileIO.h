#ifndef FILEIO
#define FILEIO
#include <iostream>
#include <fstream>
#include <stdio.h> 

void Write1DFunc(double* X, double* Y, int N, const char* fname){
    FILE* fp = fopen(fname, "w");
    if(fp == NULL){
        printf("Cannot open a new file for writing functions\n");
        exit(0);
    }
    for(int i=0; i<N; i++){
        fprintf(fp, "%f\t%f\n", X[i], Y[i]);
    }
}

#endif