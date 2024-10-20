#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "Element.h"
#include "CompMesh.h"
#include "Arlequin.h"
#include "LinearSolver.h"
#include "MatrixType.h"
#include "PETScMatrix.h"

// class Assemble;
class CompMesh;
class Arlequin;
class LinearSolver;

enum class SolverType
{
                EMumps, EIterative, EUmfpack, ECholmod, EKLU, ESPQR, //PETSc supported solvers
                ELU, ELLt, ELDLt, EQR, ECG, EBiCGStab, ELSCG, //Eigen supported solvers
                AMGCLEBiCGStab // AMGCL-Demidov solvers
};

class Analysis
{
protected:
    
    std::vector<CompMesh *> fMeshVector;
    SolverType fSolverType;
    Arlequin* fArlequin;
    bool fReducedArlequin;

    //PetscVariable
    MatrixType       *fGlobalMatrix;
    LinearSolver     *fSolver;

public:
    
    
    Analysis() = default;

    Analysis(CompMesh *cmesh, SolverType stype);

    Analysis(Arlequin *arl, SolverType stype, bool reduced = true);
    
    int64_t NEquations();

    MatrixType * GlobalMatrix() {return fGlobalMatrix;}

    std::vector<CompMesh *> &MeshVector(){return fMeshVector;}

    virtual void Compute() = 0;
    virtual void ComputeJacobian() = 0;
    virtual void ComputeRhs() = 0;
    virtual void UpdateSolution() = 0;
    Arlequin *ArlequinModel(){return fArlequin;}

    void Solve();
    void AllocateMonomodel();
    void AllocateArlequin();

    SolverType &SType(){return fSolverType;}
    void SetSolverType(SolverType stype){
        fSolverType = stype;
    }
    
    virtual void Run();

    void PostProcessError(VecDouble &errorsTotal);
    
    ~Analysis() = default;

    LinearSolver *Solver(){return fSolver;}
};





#endif