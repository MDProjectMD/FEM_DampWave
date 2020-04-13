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
    unsigned int Ng = 2*linWest->getSpatialNodesNumber() - 3;
    g_vector = Eigen::VectorXd (Ng);
    double V0 = linWest->V0_array[2];
    double dV0dt = (3.*linWest->V0_array[2] - 4.*linWest->V0_array[1] + linWest->V0_array[0])/2./dt;
    double ddV0dt = (linWest->V0_array[2] - 2.*linWest->V0_array[1] + linWest->V0_array[0])/dt/dt;
    double V1 = linWest->V1_array[2];
    double dV1dt = (3.*linWest->V1_array[2] - 4.*linWest->V1_array[1] + linWest->V1_array[0])/2./dt;
    double ddV1dt = (linWest->V1_array[2] - 2.*linWest->V1_array[1] + linWest->V1_array[0])/dt/dt;
    for(int i=0; i<Ng; i++){
        if(i == 0){
            g_vector(i) = 0. - mass_matrix.coeff(1, 2)*ddV0dt - linWest->A * V0 * stiff_matrix.coeff(1, 2) - linWest->B * dV0dt * stiff_matrix.coeff(1, 2);
        }else if(i == 1){
            g_vector(i) = 0. - mass_matrix.coeff(1, 3)*ddV0dt - linWest->A * V0 * stiff_matrix.coeff(1, 3) - linWest->B * dV0dt * stiff_matrix.coeff(1, 3);
        }else if(i == Ng - 1){
            g_vector(i) = 0. - mass_matrix.coeff(Ng - 3, Ng - 2)*ddV1dt - linWest->A * V1 * stiff_matrix.coeff(Ng - 3, Ng - 2) - linWest->B * dV1dt * stiff_matrix.coeff(Ng - 3, Ng - 2);
        }else if(i == Ng - 2){
            g_vector(i) = 0. - mass_matrix.coeff(Ng - 4, Ng - 2)*ddV1dt - linWest->A * V1 * stiff_matrix.coeff(Ng - 4, Ng - 2) - linWest->B * dV1dt * stiff_matrix.coeff(Ng - 4, Ng - 2);
        }else{
            g_vector(i) = 0.;
        }
    }
    stepIdx = 3;
}

void TimeEvolutionSolver::UpdateRHSVector(){
    unsigned int Ng = 2*linWest->getSpatialNodesNumber() - 3;
    double V0 = linWest->V0_array[stepIdx];
    double dV0dt = (3.*linWest->V0_array[stepIdx] - 4.*linWest->V0_array[stepIdx - 1] + linWest->V0_array[stepIdx - 2])/2./dt;
    double ddV0dt = (linWest->V0_array[stepIdx] - 2.*linWest->V0_array[stepIdx - 1] + linWest->V0_array[stepIdx - 2])/dt/dt;
    double V1 = linWest->V1_array[stepIdx];
    double dV1dt = (3.*linWest->V1_array[stepIdx] - 4.*linWest->V1_array[stepIdx - 1] + linWest->V1_array[stepIdx - 2])/2./dt;
    double ddV1dt = (linWest->V1_array[stepIdx] - 2.*linWest->V1_array[stepIdx - 1] + linWest->V1_array[stepIdx - 2])/dt/dt;
    g_vector(0) = 0. - mass_matrix.coeff(1, 2)*ddV0dt - linWest->A * V0 * stiff_matrix.coeff(1, 2) - linWest->B * dV0dt * stiff_matrix.coeff(1, 2);
    g_vector(1) = 0. - mass_matrix.coeff(1, 3)*ddV0dt - linWest->A * V0 * stiff_matrix.coeff(1, 3) - linWest->B * dV0dt * stiff_matrix.coeff(1, 3);
    g_vector(Ng - 1) = 0. - mass_matrix.coeff(Ng - 3, Ng - 2)*ddV1dt - linWest->A * V1 * stiff_matrix.coeff(Ng - 3, Ng - 2) - linWest->B * dV1dt * stiff_matrix.coeff(Ng - 3, Ng - 2);
    g_vector(Ng - 2) = 0. - mass_matrix.coeff(Ng - 4, Ng - 2)*ddV1dt - linWest->A * V1 * stiff_matrix.coeff(Ng - 4, Ng - 2) - linWest->B * dV1dt * stiff_matrix.coeff(Ng - 4, Ng - 2);
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
    unsigned int Ng = 2*linWest->getSpatialNodesNumber() - 3;
    // compute velocity at time step 1
    linWest->vel_matrix.row(1) = linWest->h_vector * dt + linWest->vel_vector;
    linWest->vel_matrix(1, 0) = linWest->V0_array[1];
    linWest->vel_matrix(1, Ng + 1) = linWest->V1_array[1];

    // compute velocity at time step 2
    Eigen::SparseMatrix<double, Eigen::RowMajor> MatA = 1./dt/dt * mass_matrix + linWest->A * stiff_matrix + 1./dt * linWest->B * stiff_matrix;
    Eigen::VectorXd Vecb = (2./dt/dt * mass_matrix + 1./dt * linWest->B * stiff_matrix)*(Eigen::VectorXd)(linWest->vel_matrix.row(1)).segment(1, Ng) - 1./dt/dt * mass_matrix * (Eigen::VectorXd)(linWest->vel_matrix.row(0)).segment(1, Ng) + g_vector;
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,  Eigen::Lower|Eigen::Upper, Eigen::DiagonalPreconditioner<double> > CG_SOLVER;
    CG_SOLVER.setMaxIterations(1000);
    CG_SOLVER.compute(MatA);
    Eigen::VectorXd X = CG_SOLVER.solve(Vecb);
    // computation details
    sol_conf.niter = CG_SOLVER.iterations();
    sol_conf.err = (MatA * X - Vecb).norm();
    sol_conf.stepIdx = 2;
    sol_conf.initial_solver_log();
    sol_conf.print_solver_log();

    Eigen::VectorXd x0(1);
    Eigen::VectorXd x_(1);
    x0 << linWest->V0_array[2];
    x_ << linWest->V1_array[2];
    Eigen::VectorXd u2_vec(2*linWest->getSpatialNodesNumber() - 1);
    u2_vec << x0, X, x_;
    linWest->vel_matrix.row(2) = u2_vec;

    // reserve the last 3 steps
    u1_vector = linWest->vel_matrix.row(2);
    u2_vector = linWest->vel_matrix.row(1);
    u3_vector = linWest->vel_matrix.row(0);
    //sol_conf.animation_record(linWest);
}

void SolverConfig::animation_record(LinearWestervelt* linWest){
    Engine* ep = StartMatlabEngine();
    Eigen:: MatrixXd MatC = (Eigen::MatrixXd) linWest->vel_matrix; // 如果是强转的话会产生临时变量，也就是右值，右值不能赋值给& ，只能赋值给CONST & 
    MatTrajectoryEigenPlot(ep, MatC, linWest->z_vector);
    CloseMatlabEngine(ep);
}

void TimeEvolutionSolver::UpdateOneStep(){
    unsigned int Ng = 2*linWest->getSpatialNodesNumber() - 3;
    UpdateRHSVector();  // Update the vector g_vector to the time step to be solved

    Eigen::SparseMatrix<double, Eigen::RowMajor> MatA = 2./dt/dt*mass_matrix + linWest->A*stiff_matrix + 3./2./dt*linWest->B*stiff_matrix;
    Eigen::SparseMatrix<double, Eigen::RowMajor> Mat1 = 5./dt/dt*mass_matrix + 2./dt*linWest->B*stiff_matrix;
    Eigen::SparseMatrix<double, Eigen::RowMajor> Mat2 = 4./dt/dt*mass_matrix + linWest->B/2./dt*stiff_matrix;
    Eigen::SparseMatrix<double, Eigen::RowMajor> Mat3 = 1./dt/dt*mass_matrix;
    Eigen::VectorXd Vecb = Mat1*u1_vector.segment(1, Ng) - Mat2*u2_vector.segment(1, Ng) + Mat3*u3_vector.segment(1, Ng) + g_vector;
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,  Eigen::Lower|Eigen::Upper, Eigen::DiagonalPreconditioner<double> > CG_SOLVER;
    CG_SOLVER.setMaxIterations(1000);
    CG_SOLVER.compute(MatA);
    Eigen::VectorXd X = CG_SOLVER.solve(Vecb);
    // computation details
    sol_conf.niter = CG_SOLVER.iterations();
    sol_conf.err = (MatA * X - Vecb).norm();
    sol_conf.stepIdx = stepIdx;
    sol_conf.print_solver_log();

    Eigen::VectorXd x0(1);
    Eigen::VectorXd x_(1);
    x0 << linWest->V0_array[stepIdx];
    x_ << linWest->V1_array[stepIdx];
    Eigen::VectorXd u2_vec(2*linWest->getSpatialNodesNumber() - 1);
    u2_vec << x0, X, x_;
    if(linWest->ifStore()){
        linWest->vel_matrix.row(stepIdx) = u2_vec;
    }
    u3_vector = u2_vector;
    u2_vector = u1_vector;
    u1_vector = u2_vec;
    stepIdx++;
}

void TimeEvolutionSolver::RecordAnimation(){
    sol_conf.animation_record(linWest);
}