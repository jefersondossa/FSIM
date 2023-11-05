#ifndef WEAKFORM_H
#define WEAKFORM_H

#include "DataTypes.h"
#include "IntPointData.h"
#include "Node.h"

class WeakForm
{
protected:
    int fMatId;
    int fNState = 1;
    int fDimension = 0;
    static double fBigNumber;
    std::function<void (const VecDouble &coord, VecDouble &u, MatrixDouble &gradU)> fExactSol = 0; 
    std::function<void (const VecDouble &coord, VecDouble &force)> fForceFunction = 0; 
    bool fPositionalFEM = false;

public:
    WeakForm() = default;

    ~WeakForm() = default;

    virtual void ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness) {};
    
    virtual void ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs) {};
    
    virtual void ComputeError(IntPointData &data, VecDouble &errors) {};
    
    void SetExactSolution(std::function<void (const VecDouble &coord, VecDouble &u, MatrixDouble &gradU)> exSol){
        fExactSol = exSol;
    }

    std::function<void (const VecDouble &coord, VecDouble &u, MatrixDouble &gradU)> &GetExactSolution(){
        return fExactSol;
    }

    void SetForcingFunction(std::function<void (const VecDouble &coord, VecDouble &force)> ffunction){
        fForceFunction = ffunction;
    }

    std::function<void (const VecDouble &coord, VecDouble &force)> &GetForcingFunction(){
        return fForceFunction;
    }

    int &Id() {return fMatId;}

    int &NState() {return fNState; }

    bool &IsPositionalFEM() {return fPositionalFEM;}

    virtual void ComputeStiffness(int &index, IntPointData &data, std::vector<MatrixDouble> &Stiffness) {};
    
    virtual void ComputeResidual(int &index, IntPointData &data, std::vector<VecDouble> &Rhs) {};
};


#endif
