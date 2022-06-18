#ifndef DIRICHLET_BC_H
#define DIRICHLET_BC_H

#include "DegreeOfFreedom.h"

class DirichletBoundaryCondition
{
    private:
        DegreeOfFreedom* dof_;
        double value_;

    public:

        DirichletBoundaryCondition(DegreeOfFreedom* const dof, const double value);

        ~DirichletBoundaryCondition();

        DegreeOfFreedom* getDegreeOfFreedom() const;

        double getValue() const;
};

#endif