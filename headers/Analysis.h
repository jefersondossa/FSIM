#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "Element.h"
#include "Boundary.h"
#include "CompMesh.h"
#include <petscksp.h> 

// class Assemble;
class CompMesh;

class Analysis
{
private:
    
    CompMesh *fCompMesh;
    SolverType fSolverType;

    //PetscVariable
    Mat               fGlobalStiffness;
    Vec               fGlobalRhs, fGlobalSolution;
    KSP               ksp;
    PC                pc;

public:
    
    
    Analysis() = default;

    Analysis(CompMesh *cmesh, SolverType stype){
        fCompMesh = cmesh;
        fSolverType = stype;
    };

    Mat &Stiffness(){return fGlobalStiffness;}
    Vec &Rhs(){return fGlobalRhs;}
    Vec &Solution(){return fGlobalSolution;}
    CompMesh *Mesh(){return fCompMesh;}

    virtual void Compute();
    virtual void UpdateSolution();

    void Solve();
    void AllocateProblem();
    
    void Run(){
        AllocateProblem();
        Compute();
        Solve();
        UpdateSolution();
    };
    
    ~Analysis();
};





#endif