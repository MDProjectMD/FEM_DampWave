#include "Matrix.h"

int main(){
    Matrix<double> A(4, 4);
    Matrix<double> b(4, 1);
    for(int i=0;i<4;i++){
        for(int j=0;j<4; j++){
            if(i==j){
                A[i][j] = 1.;
            }
        }
    }
    for(int i=0;i<4;i++){
        b[i][0] = 4.;
    }
    Matrix<double> x = A/b;
    x.Show();
}