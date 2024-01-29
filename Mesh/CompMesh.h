#ifndef CompMesh_H
#define CompMesh_H


#include "CompMeshTools.h"
#include "DataTypes.h"
#include "Element.h"
#include "DomainIntegration.h"
#include "Analysis.h"
#include "GmshTools.h"
#include "WeakForm.h"
#include <vector>
#include "GraphMesh.h"

class Element;
class Analysis;
class GraphMesh;

class CompMesh{
private:
    int fDimension = 2;
    int fNState = 0;
    int fOrder = 1;
    int nBdNodes = 0;

    /// Defines the vector of fluid nodes
    std::vector<Node *>       fNodeVector;

    /// Defines the vector of fluid elements
    std::vector<Element *>    fElementVector;

    std::map<int,WeakForm *> fMaterialVector;

    GraphMesh *fGraphMesh;
    
public:
    int* part_elem;      //Domain Decomposition - Elements
    int* part_nodes;     //Domain Decomposition - Nodes

    /// @brief Default constructor
    CompMesh() = default;

    /// @brief Default destructor
    ~CompMesh() = default;

    /// @brief Inserts a weak form to the computational mesh
    /// @param wf weak form
    void InsertMaterial(WeakForm *wf);

    /// @brief Returns the weak for associated with a given material id
    /// @param matid material id
    /// @return weak form object
    WeakForm* Material(int matid){
        return fMaterialVector[matid];
    }

    std::map<int,WeakForm *> &MaterialVector(){
        return fMaterialVector;
    }
    
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

    void SetDefaultOrder(int order){
        fOrder = order;
    }
    void SetDimension(int dim){
        fDimension = dim;
    }
    void SetNStateVariables(int nstate){
        fNState = nstate;
    }

    /// @brief Returns the problem dimension
    /// @return dimension
    int &Dimension() {
        return fDimension;
    }
    
    /// @brief returns the polynomial order
    /// @return polynominal order
    int &GetDefaultOrder() {
        return fOrder;
    }

    /// @brief returns the number of state variables
    /// @return state variables
    int &NState(){
        return fNState;
    }

    /// @brief Returns the number of boundary nodes
    /// @return element boundary nodes
    int &NBdNodes() {
        return nBdNodes;
    }

    /// @brief returns the number of global DOF's
    /// @return number of DOF's
    int64_t NGlobalDOF() {
        return fNodeVector.size()*fNState;
    }

    /// @brief Returns the graphycal mesh to print results.
    /// @return graphycal mesh 
    GraphMesh* GetGraphMesh();
};

#endif