#ifndef CompMesh_H
#define CompMesh_H

#include "DataTypes.h"
#include "Element.h"
#include "Boundary.h"
#include "ProblemParameters.h"
#include "DomainIntegration.h"
#include "Analysis.h"

class Element;
class Boundary;
class Analysis;

class CompMesh{
public:
    CompMesh() = default;

    int fDimension = 2;
    int fOrder = 1;
    int nBdNodes = 0;
    int nElNodes = 0;
    int nLocDOF = 0;
    int* part_elem;      //Fluid Domain Decomposition - Elements
    int* part_nodes;     //Fluid Domain Decomposition - Nodes
    int numDOF;

    /// Defines the vector of fluid nodes
    std::vector<Node *>       nodes_;

    /// Defines the vector of fluid boundaries mesh nodes
    std::vector<Boundary *> boundary_;

    /// Defines the vector of fluid elements
    std::vector<Element *>   elements_;

    ProblemParameters fProbParameters;
    DomainIntegration* numIntegration; //Numerical integration
    ProblemType fProbType = ProblemType::ENavierStokes;

    /// Gets the fluid model nodes and export for solving the overlapping
    /// mesh problem with the Arlequin method
    /// @return fluid model nodes information
    std::vector<Node *> &NodeVec(){return nodes_;}

    /// Gets the fluid model elements and export for solving the overlapping
    /// mesh problem with the Arlequin method
    /// @return fluid model elements information
    std::vector<Element *> &ElementVec(){return elements_;}

    ProblemParameters &getProblemParameters(){
        return fProbParameters;
    }

    DomainIntegration* getNumericalIntegration(){return numIntegration;}

    int &Dimension() {return fDimension;}
    int &GetDefaultOrder() {return fOrder;}
     int &NumDOF(){return numDOF;}
};

#endif