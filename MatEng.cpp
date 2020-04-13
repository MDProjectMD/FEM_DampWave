#include "MatEng.h"

Engine* StartMatlabEngine(){
    Engine *ep;
    if (!(ep = engOpen(""))) {
		fprintf(stderr, "\nCan't start MATLAB engine\n");
		exit(0);
	}
    return ep;
}

void CloseMatlabEngine(Engine* ep){
    engClose(ep);
}

// shrinked 1D column major matrix: d_array
void MatrixCPlainToMatlab(Engine* ep, double** mat_c,int n_row, int n_col, const char* vname){
    mxArray *var_mlab = NULL;
    if(var_mlab == NULL){
        var_mlab = mxCreateDoubleMatrix(n_row, n_col, mxREAL);
    }
    double* mat_1D = (double*) calloc(n_row*n_col, sizeof(double));
    for(int j = 0; j<n_col; j++){
        for(int i=0; i<n_row; i++){
            mat_1D[j*n_row + i] = mat_c[i][j];
        }
    }
    memcpy((void *)mxGetPr(var_mlab), (void *)mat_1D, n_row*n_col*sizeof(double));
    free(mat_1D);
    engPutVariable(ep, vname, var_mlab);
    if(engGetVariable(ep, vname) == NULL){
        printf("Matrix transfer fails!\n");
        exit(0);
    }
}

void MatTrajectoryEigenPlot(Engine* ep, Eigen::MatrixXd& mat_traj, Eigen::VectorXd& X){
    // n_col must be equal to X length
    int n_row = mat_traj.rows();
    int n_col = mat_traj.cols();

    mxArray *var_mlab = NULL;
    if(var_mlab == NULL){
        var_mlab = mxCreateDoubleMatrix(n_row, n_col, mxREAL);
    }
    memcpy((void *)mxGetPr(var_mlab), (void *)mat_traj.data(), n_row*n_col*sizeof(double));
    engPutVariable(ep, "mat_traj", var_mlab);
    if(engGetVariable(ep, "mat_traj") == NULL){
        printf("Matrix transfer fails!\n");
        exit(0);
    }

    mxArray *var_X = NULL;
    if(var_X == NULL){
        var_X = mxCreateDoubleMatrix(1, n_col, mxREAL);
    }
    memcpy((void *)mxGetPr(var_X), (void *)X.data(), n_col*sizeof(double));
    engPutVariable(ep, "Z", var_X);
    if(engGetVariable(ep, "Z") == NULL){
        printf("Matrix transfer fails!\n");
        exit(0);
    }

    char matlab_cmd[] = "plt = plot(Z, mat_traj(i, :), '.r', 'MarkerSize', 10);";
    // for loop should be placed in one string
    engEvalString(ep, "writerObj=VideoWriter('test.avi');writerObj.FrameRate = 5; writerObj.Quality = 100; open(writerObj);");
    engEvalString(ep, "for i = 1:1:size(mat_traj, 1) plot(Z, mat_traj(i, :), '.r', 'MarkerSize', 10); axis([-0.5 22 -1.2 1.2]); frame = getframe(gcf); writeVideo(writerObj,frame); end");
    engEvalString(ep, "close(writerObj);");
    //engEvalString(ep, "saveas(plt, 'frame2.png');");
}

/*
"tidx = 1;
                        Y = mat_traj(tidx, :);  = plot(Z, Y, '.r', 'MarkerSize', 10); plt.YDataSource = 'Y';
                        for n = tidx:1:3
                            Y = mat_traj(tidx, :);
                            refreshdata;
                        end ";
*/


