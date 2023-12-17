#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "Element.h"
#include "CompMesh.h"
#include "Arlequin.h"
#include <petscksp.h> 

// class Assemble;
class CompMesh;
class Arlequin;

enum SolverType{EMumps, EIterative, EUmfpack, ECholmod, EKLU, ESPQR};

class Analysis
{
private:
    
    std::vector<CompMesh *> fMeshVector;
    SolverType fSolverType;
    Arlequin* fArlequin;

    //PetscVariable
    Mat               fGlobalStiffness;
    Vec               fGlobalRhs, fGlobalSolution;
    KSP               ksp;
    PC                pc;

public:
    
    
    Analysis() = default;

    Analysis(CompMesh *cmesh, SolverType stype){
        fMeshVector.resize(1);
        fMeshVector[0] = cmesh;
        fSolverType = stype;
        AllocateMonomodel();
    };

    Analysis(Arlequin *arl, SolverType stype);
    
    int64_t NEquations();

    Mat &Stiffness(){return fGlobalStiffness;}
    Vec &Rhs(){return fGlobalRhs;}
    Vec &Solution(){return fGlobalSolution;}
    std::vector<CompMesh *> &MeshVector(){return fMeshVector;}

    virtual void Compute() = 0;
    virtual void UpdateSolution() = 0;
    Arlequin *ArlequinModel(){return fArlequin;}

    void Solve();
    void AllocateMonomodel();
    void AllocateArlequin();
    
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