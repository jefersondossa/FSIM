#ifndef CompMesh_H
#define CompMesh_H


#include "CompMeshTools.h"
#include "DataTypes.h"
#include "Element.h"
#include "ProblemParameters.h"
#include "DomainIntegration.h"
#include "Analysis.h"
#include "GmshTools.h"

class Element;
class Analysis;

class CompMesh{
private:
    int fDimension = 2;
    int fNState = 0;
    int fOrder = 1;
    int nBdNodes = 0;
    // int nElNodes = 0;
    int numDOF = 0;

    /// Defines the vector of fluid nodes
    std::vector<Node *>       fNodeVector;

    /// Defines the vector of fluid boundaries mesh nodes
    std::vector<Element *>    fBoundaryVector;

    /// Defines the vector of fluid elements
    std::vector<Element *>    fElementVector;

public:
    CompMesh() = default;
    CompMesh(ProblemParameters &pparam, int DIM, int order){
        fDimension = DIM;
        fProbParameters = pparam;
        fOrder = order;
        nBdNodes = 3*(1-fOrder)+DIM*(2*fOrder-1);
        if (fProbParameters.ProbType() == ProblemType::ENavierStokes || fProbParameters.ProbType() == ProblemType::EStokes){
            fNState = (DIM+1);
        } else if (fProbParameters.ProbType() == ProblemType::EPoisson) {
            fNState = 1;
        } else if (fProbParameters.ProbType() == ProblemType::EElastic || fProbParameters.ProbType() == ProblemType::ESolidPositional){
            fNState = DIM;
        } else {
            PanicButton();
        }
    };
    
    int* part_elem;      //Fluid Domain Decomposition - Elements
    int* part_nodes;     //Fluid Domain Decomposition - Nodes
    

    ProblemParameters fProbParameters;
    
    /// Gets the fluid model nodes and export for solving the overlapping
    /// mesh problem with the Arlequin method
    /// @return fluid model nodes information
    std::vector<Node *> &NodeVec(){return fNodeVector;}
    int64_t NNodes(){return fNodeVector.size();}
    void SetNumNodes(int64_t nnodes){fNodeVector.resize(nnodes);}

    /// Gets the fluid model elements and export for solving the overlapping
    /// mesh problem with the Arlequin method
    /// @return fluid model elements information
    std::vector<Element *> &ElementVec(){return fElementVector;}
    int64_t NElements(){return fElementVector.size();}
    void SetNumElements(int64_t nelements){fElementVector.resize(nelements);}

    /// Gets the fluid model elements and export for solving the overlapping
    /// mesh problem with the Arlequin method
    /// @return fluid model elements information
    std::vector<Element *> &BoundaryVec(){return fBoundaryVector;}
    int64_t NBoundElements(){return fBoundaryVector.size();}

    ProblemParameters &getProblemParameters(){
        return fProbParameters;
    }

    void SetDefaultOrder(int order){
        fOrder = order;
    }
    void SetDimension(int dim){
        fDimension = dim;
    }
    void SetNStateVariables(int nstate){
        fNState = nstate;
    }

    int &Dimension() {return fDimension;}
    int &GetDefaultOrder() {return fOrder;}
    int &NState(){return fNState;}
    // int &NElNodes() {return nElNodes;}
    int &NBdNodes() {return nBdNodes;}
    int &NGlobalDOF() {return numDOF;}
};

#endif