#include "Stochastic.h"


void GenerateStochasticVelocitySeries(double** vel_array, double** t_array, int N, double dt, bool show_tag){
    (*vel_array) = (double*) calloc(N, sizeof(double));
    if((*t_array) == NULL){
        (*t_array) = (double*) calloc(N, sizeof(double));
    }
    double scalar = sqrt(T);
    (*vel_array)[0] = 0;
    double s1 = exp(0.-lbda*dt);
    double s2 = sqrt(1. - s1*s1);
    for(int i=1; i<N; i++){
        (*vel_array)[i] = s1*(*vel_array)[i-1] + s2*GaussianDistribution(0, 1);
    }

    for(int i=0; i<N; i++){
        (*vel_array)[i] *= scalar;
        (*t_array)[i] = i*dt;
    }

    if(show_tag){
        double mu = 0.;
        double sigma = 0.;
        for(int i=0; i<N; i++){
            mu += (*vel_array)[i];
            sigma += (*vel_array)[i]*(*vel_array)[i];
        }
        mu /= N;
        sigma /= N;
        printf("Mean mu:\t%f\tVariance sigma:\t%f\n", mu, sigma);
    }
}
