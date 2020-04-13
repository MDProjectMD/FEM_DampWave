#ifndef MATLIB
#define MATLIB

#define EIGENLIBRARY

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "engine.h"
#include <iostream>

#ifdef EIGENLIBRARY
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <vector>
#endif

Engine* StartMatlabEngine();

void CloseMatlabEngine(Engine* ep);

void MatrixCPlainToMatlab(Engine* ep, double** mat_c,int M, int N, const char* vname);

void MatTrajectoryEigenPlot(Engine* ep, Eigen::MatrixXd& mat_traj, Eigen::VectorXd& X);

#endif