#ifndef CompMesh_H
#define CompMesh_H


#include "CompMeshTools.h"
#include "DataTypes.h"
#include "Analysis.h"
#include "GmshTools.h"
#include "WeakForm.h"
#include <vector>
#include <set>
#include "GraphMesh.h"
#include "Connect.h"
#include "GeoMesh.h"

enum class ApproxType
{
    EIsoparametric = 0,
    EHierarquic = 1,
    EIncompleteHierarquic = 2,
    EIsogeometric = 3
};

class Element;
class Analysis;
class GraphMesh;
class GeoMesh;
class WeakForm;

class CompMesh{
protected:
    int fDimension = 2;
    int fNState = 0;

    /// @brief Default polynomial order
    int fOrder = 1;

    ///Number of DOFS
    int fNGlobalDOF = 0;

    // /// Defines the vector of fluid nodes
    // std::vector<Node *>       fNodeVector;

    /// Defines the vector of fluid elements
    std::vector<Element *>    fElementVector;

    /// Connect vector
    std::vector<Connect *>    fConnectVector;

    std::map<int,WeakForm *> fMaterialVector;

    GraphMesh *fGraphMesh = nullptr;

    ApproxType fApproxType = ApproxType::EIsoparametric;

    GeoMesh* fReference = nullptr;

    bool fIsDisconnected = false;
    
public:
    int* part_elem;      //Domain Decomposition - Elements
    int* part_nodes;     //Domain Decomposition - Nodes

    /// @brief Default constructor
    CompMesh() = default;

    /// @brief Constructor
    CompMesh(GeoMesh *gmesh, ApproxType approxType);

    /// @brief Default destructor
    ~CompMesh() = default;

    GeoMesh* Reference() {return fReference;}

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
    

    std::vector<Connect *> &ConnectVec(){return fConnectVector;}
    int NConnects(){return fConnectVector.size();}
    void SetNumConnects(int nconnects){fConnectVector.resize(nconnects);}

    /// Gets the fluid model elements and export for solving the overlapping
    /// mesh problem with the Arlequin method
    /// @return fluid model elements information
    std::vector<Element *> &ElementVec(){return fElementVector;}
    int NElements(){return fElementVector.size();}
    void SetNumElements(int nelements){fElementVector.resize(nelements);}

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

    /// @brief returns the number of global DOF's
    /// @return number of DOF's
    int &NGlobalDOF() {
        return fNGlobalDOF;
    }

    /// @brief Returns the graphycal mesh to print results.
    /// @return graphycal mesh 
    GraphMesh* GetGraphMesh();
    GraphMesh* GetGraphMesh(std::set<int> &elsToPrint);
    void ClearGraphMesh(){
        if (fGraphMesh){
            delete fGraphMesh;
            fGraphMesh = nullptr;
        }
    }

    void Integrate(std::set<int> &matIds, std::vector<std::string> &varNames, std::map<std::string,VecDouble> &result);

    void SetSolution(VecDouble &sol);

    void SetApproxType(ApproxType type){
        fApproxType = type;
    }

    ApproxType &GetApproxType(){
        return fApproxType;
    }

    virtual void AutoBuild();
    void BuildElements();
    void BuildConnectivity();
    void BuildConnects();
    void BuildHierarquicConnects();
    void BuildIsoparametricConnects();

    void Print(std::string filename);

    void CreateDisconnectedElements(){
        fIsDisconnected = true;
    }
    
};


#endif