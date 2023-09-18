#ifndef CompMesh_H
#define CompMesh_H

#include "DataTypes.h"
#include "Boundary.h"
#include "FluidParameters.h"
#include "DomainIntegration.h"

class Boundary;

class CompMesh{
public:
    CompMesh() = default;

    int fDimension = 2;
    int fOrder = 2;
    int nBdNodes = 0;
    int nElNodes = 0;
    int nLocDOF = 0;

    /// Defines the vector of fluid nodes
    std::vector<Node *>       nodes_;

    /// Defines the vector of fluid boundaries mesh nodes
    std::vector<Boundary *> boundary_;

    FluidParameters fluidParameters;
    DomainIntegration* numIntegration; //Numerical integration
    ProblemType fProbType = ProblemType::ENavierStokes;

    /// Gets the fluid model nodes and export for solving the overlapping
    /// mesh problem with the Arlequin method
    /// @return fluid model nodes information
    std::vector<Node *> &getNodes(){return nodes_;}

    FluidParameters &getFluidParameters(){
        return fluidParameters;
    }

    DomainIntegration* getNumericalIntegration(){return numIntegration;}

    int &Dimension() {return fDimension;}
    int &GetDefaultOrder() {return fOrder;}
};

#endif