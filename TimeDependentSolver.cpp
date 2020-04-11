#include "TimeDependentSolver.h"

void TimeEvolutionSolver::AssembleStiffnessMatrix(){
    // allocate memory  0,1,2,...,Nz-1  0,  1/2,1,3/2,2,...,Nz-2,Nz-3/2,    Nz-1
    unsigned int Nstiff = 2*linWest->getSpatialNodesNumber() - 3;     // stiffness matrix edge size
    printf("Stiffness matrix size %d * %d\n", Nstiff, Nstiff);
    stiff_matrix = Eigen::SparseMatrix<double, Eigen::RowMajor>(Nstiff, Nstiff);
    std::vector< Eigen::Triplet<double> > coefficients;
    for(int i=0; i<Nstiff; i++){
        if(i%2 == 0){ // even: half points position under Dirichlet BCs
            unsigned int j = i - 1;
            if(j>=0 && j<Nstiff){
                coefficients.push_back(Eigen::Triplet<double>(i, j, 0.-8./3./dz));
            }
            j = i + 1;
            if(j>=0 && j<Nstiff){
                coefficients.push_back(Eigen::Triplet<double>(i, j, 0.-8./3./dz));
            }
            j = i;
            coefficients.push_back(Eigen::Triplet<double>(i, j, 16./3./dz));
        }else{     // odd: integer grid points
            unsigned int j = i - 2;
            if(j>=0 && j<Nstiff){
                coefficients.push_back(Eigen::Triplet<double>(i, j, 1./3./dz));
            }
            j = i - 1;
            if(j>=0 && j<Nstiff){
                coefficients.push_back(Eigen::Triplet<double>(i, j, 0.-8./3./dz));
            }
            j = i;
            coefficients.push_back(Eigen::Triplet<double>(i, j, 14./3./dz));
            j = i + 1;
            if(j>=0 && j<Nstiff){
                coefficients.push_back(Eigen::Triplet<double>(i, j, 0.-8./3./dz));
            }
            j = i + 2;
            if(j>=0 && j<Nstiff){
                coefficients.push_back(Eigen::Triplet<double>(i, j, 1./3./dz));
            }
        }
    }
    stiff_matrix.setFromTriplets(coefficients.begin(), coefficients.end());
}

void TimeEvolutionSolver::AssembleMassMatrix(){
    unsigned int Nmass = 2*linWest->getSpatialNodesNumber() - 3;
    printf("Mass matrix size %d * %d\n", Nmass, Nmass);
    mass_matrix = Eigen::SparseMatrix<double, Eigen::RowMajor>(Nmass, Nmass);
    std::vector< Eigen::Triplet<double> > coefficients;
    for(int i=0; i<Nmass; i++){
        if(i%2 == 0){   // even: half points position under Dirichlet BCs
            unsigned int j = i - 1;
            if(j>=0 && j<Nmass){
                coefficients.push_back(Eigen::Triplet<double>(i, j, 1./15.*dz));
            }
            j = i + 1;
            if(j>=0 && j<Nmass){
                coefficients.push_back(Eigen::Triplet<double>(i, j, 1./15.*dz));
            }
            j = i;
            coefficients.push_back(Eigen::Triplet<double>(i, j, 8./15.*dz));
        }else{  // odd: integer grid points
            unsigned int j = i - 2;
            if(j>=0 && j<Nmass){
                coefficients.push_back(Eigen::Triplet<double>(i, j, 0.-1./30.*dz));
            }
            j = i - 1;
            if(j>=0 && j<Nmass){
                coefficients.push_back(Eigen::Triplet<double>(i, j, 1./15.*dz));
            }
            j = i;
            coefficients.push_back(Eigen::Triplet<double>(i, j, 4./15.*dz));
            j = i + 1;
            if(j>=0 && j<Nmass){
                coefficients.push_back(Eigen::Triplet<double>(i, j, 1./15.*dz));
            }
            j = i + 2;
            if(j>=0 && j<Nmass){
                coefficients.push_back(Eigen::Triplet<double>(i, j, 0.-1./30.*dz));
            }
        }
    }
    mass_matrix.setFromTriplets(coefficients.begin(), coefficients.end());
}

void TimeEvolutionSolver::InitialRHSVector(){ // rhs vector g at time step 2 (initial time step is 0)
    stepIdx = 3;
    unsigned int Ng = 2*linWest->getSpatialNodesNumber() - 3;
    g_vector = Eigen::VectorXd (Ng);
    double V0 = linWest->V0_array[2];
    double dV0dt = (3.*linWest->V0_array[2] - 4.*linWest->V0_array[1] + linWest->V0_array[0])/2./dt;
    double ddV0dt = (linWest->V0_array[2] - 2.*linWest->V0_array[1] + linWest->V0_array[0])/dt/dt;
    for(int i=0; i<Ng; i++){
        if(i == 0){
            g_vector(i) = 0. - mass_matrix.coeff(1, 2)*ddV0dt - linWest->A * V0 * stiff_matrix.coeff(1, 2) - linWest->B * dV0dt * stiff_matrix.coeff(1, 2);
        }else if(i == 1){
            g_vector(i) = 0. - mass_matrix.coeff(1, 3)*ddV0dt - linWest->A * V0 * stiff_matrix.coeff(1, 3) - linWest->B * dV0dt * stiff_matrix.coeff(1, 3);
        }else{
            g_vector(i) = 0.;
        }
    }
}

void TimeEvolutionSolver::UpdateRHSVector(){
    unsigned int Ng = 2*linWest->getSpatialNodesNumber() - 3;
    double V0 = linWest->V0_array[stepIdx];
    double dV0dt = (3.*linWest->V0_array[stepIdx] - 4.*linWest->V0_array[stepIdx - 1] + linWest->V0_array[stepIdx - 2])/2./dt;
    double ddV0dt = (linWest->V0_array[stepIdx] - 2.*linWest->V0_array[stepIdx - 1] + linWest->V0_array[stepIdx - 2])/dt/dt;
    g_vector(0) = 0. - mass_matrix.coeff(1, 2)*ddV0dt - linWest->A * V0 * stiff_matrix.coeff(1, 2) - linWest->B * dV0dt * stiff_matrix.coeff(1, 2);
    g_vector(1) = 0. - mass_matrix.coeff(1, 3)*ddV0dt - linWest->A * V0 * stiff_matrix.coeff(1, 3) - linWest->B * dV0dt * stiff_matrix.coeff(1, 3);
}

TimeEvolutionSolver::TimeEvolutionSolver(LinearWestervelt* lwptr){
    linWest = lwptr;
}

void TimeEvolutionSolver::BuildUp(){
    AssembleStiffnessMatrix();
    AssembleMassMatrix();
    // Called after above two functions
    InitialRHSVector();
}

void TimeEvolutionSolver::InitializeScheme(){
    // compute velocity at time step 1
    linWest->vel_matrix.row(1) = linWest->h_vector * dt + linWest->vel_vector;

}