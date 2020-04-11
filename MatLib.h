#ifndef MATLIB
#define MATLIB

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "engine.h"

Engine* StartMatlabEngine(){
    Engine *ep;
    if (!(ep = engOpen(""))) {
		fprintf(stderr, "\nCan't start MATLAB engine\n");
		exit(0);
	}
    return ep;
}

void CloseMatlabEngine(Engine* ep){
    engClose(ep);
}

// shrinked 1D column major matrix: d_array
void MatrixMemcpyCToMatlab(Engine* ep, double* d_array,int M, int N, const char* vname){
    mxArray *var_mlab = NULL;
    var_mlab = mxCreateDoubleMatrix(M, N, mxREAL);
    memcpy((void *)mxGetPr(var_mlab), (void *)d_array, M*N*sizeof(double));
    engPutVariable(ep, vname, var_mlab);
    if(engGetVariable(ep, vname) == NULL){
        printf("Matrix transfer fails!\n");
        exit(0);
    }
}


#endif

