#ifndef EIGEN_LINEARSOLVER_H
#define EIGEN_LINEARSOLVER_H

#include "LinearSolver.h"

class EigenLinearSolver : public LinearSolver
{
private:
    /* data */
public:
    EigenLinearSolver(Analysis *an);
    ~EigenLinearSolver();

    void Solve() override;
};




#endif