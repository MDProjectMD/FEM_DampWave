#include "Stochastic.h"
#include "Define.h"
#include "FileIO.h"

int main(){
    double* vel_array, * time_array;
    GenerateStochasticVelocitySeries(&vel_array, &time_array, Nt, dt);
    Write1DFunc(time_array, vel_array, Nt, "velz_short_large_dt");
}