#ifndef SOLVER
#define SOLVER

#include "SystemDefine.h"
#include "MatEng.h"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <vector>


class SolverConfig{
public:
    unsigned int stepIdx;
    unsigned int niter;
    double err;
public:
    void initial_solver_log(){
        printf("time_step\titerations\terrors\n");
    }
    void print_solver_log(){
        printf("%d\t%d\t%e\n", stepIdx, niter, err);
    }
    void animation_record(LinearWestervelt* linWest);   // ifStoreSolution must be set True !!
};

// Based on Eigen C++ library
class TimeEvolutionSolver{
private:
    Eigen::VectorXd u1_vector;  // solution at step n - 1
    Eigen::VectorXd u2_vector;  // solution at step n - 2
    Eigen::VectorXd u3_vector;  // solution at step n - 3

    Eigen::SparseMatrix<double, Eigen::RowMajor> stiff_matrix;
    Eigen::SparseMatrix<double, Eigen::RowMajor> mass_matrix;
    Eigen::VectorXd g_vector;

    void InitialRHSVector();
    void UpdateRHSVector();
    void AssembleStiffnessMatrix();
    void AssembleMassMatrix();
    LinearWestervelt* linWest;
    SolverConfig sol_conf;
    unsigned int stepIdx;   // indicates the latest time index to be solved
public:
    TimeEvolutionSolver(LinearWestervelt* lwptr);   
    void BuildUp();     // prepare mass/stifness matrix
    void InitializeScheme();  // initializing first few steps using BCs
    void UpdateOneStep();   // update ahead one time step
    void RecordAnimation();
};

#endif