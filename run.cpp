#include "Matrix.h"
#include "MatEng.h"
#include "Stochastic.h"
#include "Define.h"
#include "FileIO.h"
#include "SystemDefine.h"
#include "TimeDependentSolver.h"

int main(){
    LinearWestervelt lw;
    lw.SetStoreSolution(true);
    lw.SetDirichlet0(V0);
    lw.SetDirichlet1(V1);
    lw.SetInitial0(g); 
    lw.SetInitial1(h);
    lw.BuildUp();

    TimeEvolutionSolver solver(&lw);
    solver.BuildUp();
    solver.InitializeScheme();
    for(int i=3; i<lw.getTemporalNodesNumber(); i++){
        solver.UpdateOneStep();
    }
    //solver.UpdateOneStep();
    solver.RecordAnimation();
}