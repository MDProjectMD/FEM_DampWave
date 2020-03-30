#ifndef MATRIX
#define MATRIX
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#define MATLAB_ENGINE

#ifdef MATLAB_ENGINE
#include "engine.h"
#endif


// if Template is USED, then the implementation and the declaration should be palced in the SAME FILE !!!!
template <typename Scalar> 
class Matrix{
private:
    unsigned int M; // row
    unsigned int N; // column
    Scalar** ptr;
    Matrix(); // default empty object
    void AllocateMemory(){
        ptr = (Scalar**) calloc(M, sizeof(Scalar*));
        for(int i=0; i<M; i++){
            ptr[i] = (Scalar*) calloc(N, sizeof(Scalar));
        }
    }
public:
    Matrix(unsigned int m, unsigned int n); // initialize a matrix with all elements to be 0
    Matrix(unsigned int m, unsigned int n, Scalar* darray, char mode = 'r'); // assemble a given 2D array, 'r' -- row-major & 'c' -- column-major
    // for *.mat file, matlab stores matrix in column-major sequence
    void Free(); // HERE IF using deconstructor, then 'operator*' CANNOT return a object!
    Scalar* const operator[] (const unsigned int idx);
    Matrix<Scalar> operator* (Matrix<Scalar> M); // pairwise element product
    #ifdef MATLAB_ENGINE
    Matrix<Scalar> operator/ (Matrix<Scalar> b); // Only defined for double-type, based on Matlab library
    #endif
    unsigned int getRowNum();
    unsigned int getColumnNum();
    Scalar sum();
    void Show();
};

template <typename Scalar>
Scalar Matrix<Scalar>::sum(){
    Scalar ac = 0.;
    for(int i=0; i<M; i++){
        for(int j=0; j<N; j++){
            ac += ptr[i][j];
        }
    }
    return ac;
}

template <typename Scalar>
void Matrix<Scalar>::Show(){
    for(int i=0; i<M; i++){
        for(int j=0; j<N; j++){
            printf("%.3f\t", ptr[i][j]);
        }
        std::cout << std::endl;
    }
}

#ifdef MATLAB_ENGINE
// A*x = b  --->    x = A/b
template <typename Scalar>
Matrix<Scalar> Matrix<Scalar>::operator/ (Matrix<Scalar> b_vec){
    // first check the dimension consistency
    int Mb = b_vec.getRowNum();
    int Nb = b_vec.getColumnNum();
    if(M!=N || Nb!=1 || M!=Mb){
        std::cout << "Error in dimension match! A is 2D matrix and b is a column vector" << std::endl;
        exit(0);
    }
    // start a Matlab engine 
    Engine *ep;
	mxArray *A_mlab = NULL, *b_mlab = NULL, *x_mlab = NULL;
    if (!(ep = engOpen(""))) {
		fprintf(stderr, "\nCan't start MATLAB engine\n");
		exit(0);
	}
    A_mlab = mxCreateDoubleMatrix(M, N, mxREAL);
    b_mlab = mxCreateDoubleMatrix(Mb, Nb, mxREAL);
    double* ptr1D_col_A = (double*) calloc(M*N, sizeof(double));
    double* ptr1D_col_b = (double*) calloc(Mb, sizeof(double));
    for(int j=0; j<N; j++){
        for(int i=0; i<M; i++){
            int idx = j*M + i;
            ptr1D_col_A[idx] = ptr[i][j];
        }
    }
    for(int i=0; i<Mb; i++){
        ptr1D_col_b[i] = b_vec[i][0];
    }
    memcpy((void *)mxGetPr(A_mlab), (void *)ptr1D_col_A, M*N*sizeof(double));
    memcpy((void *)mxGetPr(b_mlab), (void *)ptr1D_col_b, Mb*sizeof(double));
    engPutVariable(ep, "A", A_mlab);
    engPutVariable(ep, "b", b_mlab);
    free(ptr1D_col_A);
    free(ptr1D_col_b);
    engEvalString(ep, "x = mldivide(A,b);");

    // retrieve result from Matlab
    double* ptr1D_col_x = (double*) calloc(Mb, sizeof(double));
    x_mlab = engGetVariable(ep,"x");
    memcpy((void *)ptr1D_col_x, (void *)mxGetPr(x_mlab), Mb*sizeof(double));
    Matrix<double> result_x(Mb, Nb, ptr1D_col_x, 'c');
    
    // free temp memory
    free(ptr1D_col_x);
    mxDestroyArray(A_mlab);
    mxDestroyArray(b_mlab);
    mxDestroyArray(x_mlab);
    return result_x;
}
#endif

template <typename Scalar>
Matrix<Scalar> Matrix<Scalar>::operator* (Matrix<Scalar> Mat){
    unsigned int Nr = Mat.getRowNum();
    unsigned int Nc = Mat.getColumnNum();
    if (M!=Nr || N!=Nc){
        std::cout<< "element pairwise product must have same dimensions " << std::endl;
        exit(0);
    }
    Matrix<Scalar> R(M, N);
    for(int i=0; i<Nr; i++){
        for(int j=0; j<Nc; j++){
            R[i][j] = ptr[i][j] * Mat[i][j];
        }
    }
    return R;
}

template <typename Scalar>
Matrix<Scalar>::Matrix(unsigned int m, unsigned int n, Scalar* darray, char mode){
    M = m;
    N = n;
    AllocateMemory();
    if(mode == 'r'){
        for(int i=0; i<M; i++){
            for(int j=0; j<N; j++){
                ptr[i][j] = darray[i*N+j];
            }
        }
    }else{
        for(int i=0; i<M; i++){
            for(int j=0; j<N; j++){
                ptr[i][j] = darray[j*M+i];
            }
        }
    }
}

template <typename Scalar>
Matrix<Scalar>::Matrix(unsigned int m, unsigned int n){
    M = m;
    N = n;
    AllocateMemory();
}

template<typename Scalar>
void Matrix<Scalar>::Free(){
    //std::cout<<"matrix ptr free"<<std::endl;
    for(int i=0; i<M; i++){
        free(ptr[i]);
    }
    free(ptr);
}

template<typename Scalar>
Scalar* const Matrix<Scalar>::operator[](const unsigned int idx){
    return ptr[idx];
}

template<typename Scalar>
unsigned int Matrix<Scalar>::getRowNum(){
    return M;
}

template<typename Scalar>
unsigned int Matrix<Scalar>::getColumnNum(){
    return N;
}


//class SparseMatrix



#endif