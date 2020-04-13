#include "SystemDefine.h"

LinearWestervelt::LinearWestervelt(){
    TMAX = TIME;
    ZMAX = L;
    N_t = Nt;
    N_z = Nz;
    A = c0*c0;
    B = 2*eta/rho;
}

unsigned int LinearWestervelt::getSpatialNodesNumber(){
    return N_z;
}

unsigned int LinearWestervelt::getTemporalNodesNumber(){
    return N_t;
}

double LinearWestervelt::getSpatialLength(){
    return ZMAX;
}

double LinearWestervelt::getTemporalLength(){
    return TMAX;
}

bool LinearWestervelt::ifStore(){
    return ifStoreSolution;
}

void LinearWestervelt::SetStoreSolution(bool tag){
    ifStoreSolution = tag;
}

void LinearWestervelt::AllocateMemory(double** ptr, int N){
    (*ptr) = (double*) calloc(N, sizeof(double));
}

void LinearWestervelt::BuildUp(){
    if(ifStoreSolution){
        vel_matrix = Eigen::MatrixXd (N_t, 2*N_z - 1);    // transfer to rowMajor matrix
    }else{
        vel_matrix = Eigen::MatrixXd (3, 2*N_z - 1);
    }

    AllocateMemory(&t_array, N_t);
    for(int i=0; i<N_t; i++){
        t_array[i] = i*dt;
    }
    t_vector = Eigen::Map<Eigen::VectorXd> (t_array, N_t); // must specify the size

    AllocateMemory(&z_array, 2*N_z - 1);
    for(int i=0; i<2*N_z - 1; i++){
        z_array[i] = i*dz/2.;
    }
    z_vector = Eigen::Map<Eigen::VectorXd> (z_array, 2*N_z - 1); // must specify the size

    AllocateMemory(&vel_array, 2*N_z - 1);
    memcpy(vel_array, g_array, (2*N_z - 1)*sizeof(double)); // velocity array at time step 0 is equal to BC Initial0 g(z)
    vel_vector = Eigen::Map<Eigen::VectorXd> (vel_array, 2*N_z - 1);
    // update vel_matrix information
    vel_matrix.row(0) = vel_vector;
}

// *********** Boundary Condition *************

void LinearWestervelt::SetDirichlet0(double(*func)(double t)){
    if(func == NULL){
        GenerateStochasticVelocitySeries(&V0_array, NULL, N_t, dt);
    }else{
        AllocateMemory(&V0_array, N_t);
        for(int i=0; i<N_t; i++){
            V0_array[i] = func(i*dt);
        }
    }
    printf("%f\t%f\t%f\n", V0_array[0], V0_array[1], V0_array[2]);
    V0_vector = Eigen::Map<Eigen::VectorXd> (V0_array, N_t);
}

void LinearWestervelt::SetDirichlet1(double(*func)(double t)){
    AllocateMemory(&V1_array, N_t);
    if(func == NULL){
        GenerateStochasticVelocitySeries(&V1_array, &t_array, N_t, dt);
    }else{
        AllocateMemory(&V1_array, N_t);
        for(int i=0; i<N_t; i++){
            V1_array[i] = func(i*dt);
        }
    }
    V1_vector = Eigen::Map<Eigen::VectorXd> (V1_array, N_t);
}

void LinearWestervelt::SetInitial0(double(*func)(double z)){
    AllocateMemory(&g_array, 2*N_z - 1);
    for(int i=0; i<2*N_z - 1; i++){
        g_array[i] = func(i*dz);
    }
    g_vector = Eigen::Map<Eigen::VectorXd> (g_array, 2*N_z - 1);
}

void LinearWestervelt::SetInitial1(double(*func)(double z)){
    AllocateMemory(&h_array, 2*N_z - 1);
    for(int i=0; i<2*N_z - 1; i++){
        h_array[i] = func(i*dz);
    }
    h_vector = Eigen::Map<Eigen::VectorXd> (h_array, 2*N_z - 1);
}

void LinearWestervelt::Finish(){
    // .........
    ;
}


double V0(double t){
    return sin(t);
}

double V1(double t){
    return 0.;
    //return sin(t);
}

double g(double z){
    return 0.;
}

double h(double z){
    return 0.;
}
