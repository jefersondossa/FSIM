#ifndef ArlequinRedSolverEigen_H
#define ArlequinRedSolverEigen_H

#include "LinearSolver.h"
#include "ArlequinMatRedEigen.h"

class ArlequinRedSolverEigen : public LinearSolver
{
private:
    /* data */
public:
    ArlequinRedSolverEigen(Analysis *an);
    virtual ~ArlequinRedSolverEigen();

    void Solve() override;
};




#endif