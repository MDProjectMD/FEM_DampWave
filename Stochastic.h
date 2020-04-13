#ifndef STOCHASTIC
#define STOCHASTIC

#include "Random.h"
#include "Define.h"

// vel_array: C plain arrays    t_array:    C plain arrays for discrete time series , if t_array = NULL, t_array will be ignored 
void GenerateStochasticVelocitySeries(double** vel_array, double** t_array, int N, double dt, bool show_tag = true);

#endif