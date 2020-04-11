#include "Matrix.h"
#include "MatLib.h"
#include "Stochastic.h"
#include "Define.h"
#include "FileIO.h"
#include "SystemDefine.h"

int main(){
    LinearWestervelt lw;
    lw.SetStoreSolution(false);
    //Eigen::Matrix<double, Dynamic, Dynamic, Eigen::RowMajor> vel_matrix = Eigen::Matrix<double, Dynamic, Dynamic, Eigen::RowMajor> (3, 5);
    Eigen::MatrixXd MtxA = Eigen::MatrixXd::Identity(4,4);
    MtxA(0,3) = 5;  MtxA(0, 2) = 2.;
    std::cout << MtxA.row(0) << std::endl;
    std::cout << MtxA.row(1) << std::endl;
    std::cout << MtxA.row(2) << std::endl;
    std::cout << MtxA.row(2) << std::endl;
    double a[4] = {0, 1, 2, 3};
    Eigen::VectorXd vec = Eigen::Map<Eigen::VectorXd> (a, 4);
    //MtxA.row(0) = vec;
    //std::cout << vec.size() << std::endl;
    //std::cout << MtxA.size() << std::endl;
    //std::cout << MtxA.row(0) << std::endl;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Acopy = MtxA;
    /*for(int i=0; i<16; i++){
        std::cout << *(Acopy.data()+i) << std::endl;
    }
    std::cout<< std::endl;
    for(int i=0; i<16; i++){
        std::cout << *(MtxA.data()+i) << std::endl;
    }*/
    // exchange row 0 and 1
    Eigen::VectorXd tmp = Acopy.row(0);
    Acopy.row(0) = Acopy.row(1);
    Acopy.row(1) = tmp;
    std::cout << Acopy << std::endl;

    //tmp = tmp * 3;
    std::cout << tmp*3 << std::endl;

    a[3] = 3.3333;
    std::cout << vec(3) << "--" << *(vec.data()+3) << std::endl;    // DOES NOT CHANGE
    std::cout << a << "--" << (vec.data()) << std::endl;


}