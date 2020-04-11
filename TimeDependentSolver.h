#ifndef SOLVER
#define SOLVER

#include "SystemDefine.h"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <vector>

// Based on Eigen C++ library
class TimeEvolutionSolver{
private:
    Eigen::SparseMatrix<double, Eigen::RowMajor> stiff_matrix;
    Eigen::SparseMatrix<double, Eigen::RowMajor> mass_matrix;
    Eigen::VectorXd g_vector;
    void InitialRHSVector();
    void UpdateRHSVector();
    void AssembleStiffnessMatrix();
    void AssembleMassMatrix();
    LinearWestervelt* linWest;
    unsigned int stepIdx;   // indicates the latest time index to be solved
public:
    TimeEvolutionSolver(LinearWestervelt* lwptr);   
    void BuildUp();     // prepare mass/stifness matrix
    void InitializeScheme();  // initializing first few steps using BCs
    void UpdateOneStep();   // update ahead one time step
};

#endif