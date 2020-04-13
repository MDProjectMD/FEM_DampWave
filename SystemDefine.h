#ifndef SYSTEMDEFINE
#define SYSTEMDEFINE
#include "Define.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Stochastic.h"

class LinearWestervelt{
private:
    bool ifStoreSolution;
    double TMAX, ZMAX;
    unsigned int N_t, N_z;
    void AllocateMemory(double** ptr, int N);
public:
    double A;
    double B;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> vel_matrix;    // discrete solutions to be stored, (time, spatial coordinate), allocated according to system size
    
    double* t_array;        // discrete time points
    Eigen::VectorXd t_vector;

    double* z_array;        // discrete spatial evaluated points, of length 2*Nz - 1
    Eigen::VectorXd z_vector;

    double* vel_array;
    Eigen::VectorXd vel_vector;  // discrete variables to be solved

    double* V0_array;
    Eigen::VectorXd V0_vector;   // discrete time series for function V0(t)

    double* V1_array;
    Eigen::VectorXd V1_vector;   // discrete time series for function V1(t)

    double* g_array;
    Eigen::VectorXd g_vector;    // discrete spatial series for BC at t=0    g(z)    0th order, of length 2*Nz - 1

    double* h_array;
    Eigen::VectorXd h_vector;    // discrete spatial series for BC at t=0    h(z)    1st order, of length 2*Nz - 1

    bool ifStore();
    unsigned int getSpatialNodesNumber();
    unsigned int getTemporalNodesNumber();
    double getSpatialLength();
    double getTemporalLength();
public:
    // Ordered as Calling Sequence
    LinearWestervelt();
    void SetStoreSolution(bool tag);    // if "True" to be set, then vel_matrix will be allocated corresponding memory; Otherwise only first FEW initial steps will be stored for computation
    // No Calling Sequence prefered
    void SetDirichlet0(double(*func)(double t));    // allocate and prepare the V0(t), if function handle not needed, just pass NULL
    void SetDirichlet1(double(*func)(double t));    // allocate and prepare the V1(t)
    void SetInitial0(double(*func)(double z));      // allocate and prepare the g(z)
    void SetInitial1(double(*func)(double z));      // allocate and prepare the h(z)

    void BuildUp();                     // final step, prepare vel_matrix, vel_array/vel_vector, t_array/t_vector

    void Finish();
};

double V0(double t);

double V1(double t);

double g(double z);

double h(double z);


#endif