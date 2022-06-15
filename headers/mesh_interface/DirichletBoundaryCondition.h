#pragma once
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


DirichletBoundaryCondition::DirichletBoundaryCondition(DegreeOfFreedom* const dof, const double value)
    : dof_(dof), value_(value) {}

DirichletBoundaryCondition::~DirichletBoundaryCondition() {}

DegreeOfFreedom* DirichletBoundaryCondition::getDegreeOfFreedom() const
{
    return dof_;
}

double DirichletBoundaryCondition::getValue() const
{
    return value_;
}
