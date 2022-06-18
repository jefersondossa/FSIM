#include "DirichletBoundaryCondition.h"


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
