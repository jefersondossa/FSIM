#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H

#include "Analysis.h"
class Analysis;

class LinearSolver
{
protected:
    Analysis *fAnalysis;

public:
    LinearSolver() = default;
    LinearSolver(Analysis *an);
    virtual ~LinearSolver() {};

    virtual void Solve() = 0;
    virtual void SolveLumped() {
        PanicButton();
    };
};


#endif
