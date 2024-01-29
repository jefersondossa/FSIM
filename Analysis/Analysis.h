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

enum SolverType{
                EMumps, EIterative, EUmfpack, ECholmod, EKLU, ESPQR, //PETSc supported solvers
                ELU, ELLt, ELDLt, EQR, ECG, EBiCGStab, ELSCG}; //Eigen supported solvers

class Analysis
{
protected:
    
    std::vector<CompMesh *> fMeshVector;
    SolverType fSolverType;
    Arlequin* fArlequin;

    //PetscVariable
    MatrixType       *fGlobalMatrix;
    LinearSolver     *fSolver;

public:
    
    
    Analysis() = default;

    Analysis(CompMesh *cmesh, SolverType stype);

    Analysis(Arlequin *arl, SolverType stype);
    
    int64_t NEquations();

    MatrixType * GlobalMatrix() {return fGlobalMatrix;}

    std::vector<CompMesh *> &MeshVector(){return fMeshVector;}

    virtual void Compute() = 0;
    virtual void UpdateSolution() = 0;
    Arlequin *ArlequinModel(){return fArlequin;}

    void Solve();
    void AllocateMonomodel();
    void AllocateArlequin();

    SolverType SType(){return fSolverType;}
    
    virtual void Run(){
        std::cout << "Allocating problem..." << std::endl;
        std::clock_t t3 = std::clock();
        Compute();
        std::clock_t t4 = std::clock();
        std::cout << "Time assembling = " << 1000.*(t4-t3)/CLOCKS_PER_SEC/1000. << "s \n";
        Solve();
        std::clock_t t5 = std::clock();
        std::cout << "Time Solving = " << 1000.*(t5-t4)/CLOCKS_PER_SEC/1000. << "s \n";
        std::cout << "Updating solution..." << std::endl;
        UpdateSolution();
    };

    void PostProcessError(VecDouble &errorsTotal);
    
    ~Analysis();
};





#endif