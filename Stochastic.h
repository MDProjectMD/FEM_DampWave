#ifndef STOCHASTIC
#define STOCHASTIC

#include "Random.h"
#include "Define.h"

void GenerateStochasticVelocitySeries(double** vel_array, double** t_array, int N, double dt, bool show_tag = true);

#endif