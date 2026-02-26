#ifndef GeoMesh_h
#define GeoMesh_h

#include "CompMesh.h"
#include "GraphMesh.h"
#include "Node.h"
#include "GeoElement.h"
#include "Connect.h"

class CompMesh;
class GraphMesh;
class GeoElement;

class GeoMesh{
private:

    int fDimension = -1;

    /// Defines the vector of fluid nodes
    std::vector<Node *>       fNodeVector;

    /// Defines the vector of fluid elements
    std::vector<GeoElement *> fElementVector;

    GraphMesh *fGraphMesh = nullptr;



public:
    GeoMesh() = default;
    ~GeoMesh() = default;

    void SetDimension(int dim){
        fDimension = dim;
    }

    int Dimension() {return fDimension;}

    /// Gets the mesh nodes
    /// @return fluid model nodes information
    std::vector<Node *> &NodeVec(){return fNodeVector;}
   
    int64_t NNodes(){return fNodeVector.size();}
    void SetNumNodes(int64_t nnodes){fNodeVector.resize(nnodes);}

    /// Gets the fluid model elements and export for solving the overlapping
    /// mesh problem with the Arlequin method
    /// @return fluid model elements information
    std::vector<GeoElement *> &ElementVec(){return fElementVector;}
    int64_t NElements(){return fElementVector.size();}
    void SetNumElements(int64_t nelements){fElementVector.resize(nelements);}

    void Print(std::string filename);


};

#endif

