#ifndef CompMesh_H
#define CompMesh_H


#include "CompMeshTools.h"
#include "DataTypes.h"
#include "Element.h"
#include "Boundary.h"
#include "ProblemParameters.h"
#include "DomainIntegration.h"
#include "Analysis.h"
#include "GmshTools.h"

class Element;
class Boundary;
class Analysis;
// namespace GmshTools{};

class CompMesh{
private:
    int fDimension = 2;
    int fNState = 0;
    int fOrder = 1;
    int nBdNodes = 0;
    int nElNodes = 0;
    int nLocDOF = 0;
    int numDOF = 0;

    /// Defines the vector of fluid nodes
    std::vector<Node *>       fNodeVector;

    /// Defines the vector of fluid boundaries mesh nodes
    std::vector<Boundary *>   fBoundaryVector;

    /// Defines the vector of fluid elements
    std::vector<Element *>    fElementVector;

public:
    CompMesh() = default;
    CompMesh(ProblemParameters &pparam, int DIM, int order){
        fDimension = DIM;
        fProbParameters = pparam;
        fOrder = order;
        nElNodes = (3+(DIM-2)*fOrder)*(2+3*fOrder+fOrder*fOrder)/6;
        nBdNodes = 3*(1-fOrder)+DIM*(2*fOrder-1);
        numIntegration = new DomainIntegration(DIM,order);
        if (fProbParameters.ProbType() == ProblemType::ENavierStokes || fProbParameters.ProbType() == ProblemType::EStokes){
            fNState = (DIM+1);
            nLocDOF = nElNodes*(DIM+1);
        } else if (fProbParameters.ProbType() == ProblemType::EPoisson) {
            fNState = 1;
            nLocDOF = nElNodes;
        } else if (fProbParameters.ProbType() == ProblemType::EElastic || fProbParameters.ProbType() == ProblemType::ESolidPositional){
            fNState = DIM;
            nLocDOF = nElNodes*DIM;
        } else {
            PanicButton();
        }

        
        
    };
    
    int* part_elem;      //Fluid Domain Decomposition - Elements
    int* part_nodes;     //Fluid Domain Decomposition - Nodes
    

    ProblemParameters fProbParameters;
    DomainIntegration* numIntegration; //Numerical integration
    
    void SetNLocDOF(int nlocdof){
        nLocDOF = nlocdof;
    }

    /// Gets the fluid model nodes and export for solving the overlapping
    /// mesh problem with the Arlequin method
    /// @return fluid model nodes information
    std::vector<Node *> &NodeVec(){return fNodeVector;}
    int64_t NNodes(){return fNodeVector.size();}

    /// Gets the fluid model elements and export for solving the overlapping
    /// mesh problem with the Arlequin method
    /// @return fluid model elements information
    std::vector<Element *> &ElementVec(){return fElementVector;}
    int64_t NElements(){return fElementVector.size();}

    /// Gets the fluid model elements and export for solving the overlapping
    /// mesh problem with the Arlequin method
    /// @return fluid model elements information
    std::vector<Boundary *> &BoundaryVec(){return fBoundaryVector;}
    int64_t NBoundElements(){return fBoundaryVector.size();}

    ProblemParameters &getProblemParameters(){
        return fProbParameters;
    }

    DomainIntegration* getNumericalIntegration(){return numIntegration;}

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
    int &NElNodes() {return nElNodes;}
    int &NBdNodes() {return nBdNodes;}
    int &NLocDOF() {return nLocDOF;}
    int &NGlobalDOF() {return numDOF;}
};

#endif