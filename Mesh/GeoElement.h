#ifndef GeoElement_h
#define GeoElement_h

#include "GeoMesh.h"

class GeoMesh;

class GeoElement
{
protected:
    GeoMesh *fMesh = nullptr;

    VecInt        fGeoNodes; // geometric nodes

    int64_t       fIndex;    //Element index
    
    VecDouble     xK, XK;    // Bounding box coordinates

    std::vector<int64_t> fNeighborElements;

    int fMaterial;

    int fPrintType;


public:
    GeoElement() = default;
    ~GeoElement() = default;

    /// Sets the element geometric nodes
    /// @param int* elenodes geometric Nodes
    void setGeometricNodes(VecInt &geonodes){fGeoNodes = geonodes;};

    /// GetGeometric nodes geometric nodes
    /// @return element geonodes
    VecInt &getGeometricNodes(){return fGeoNodes;};

    int &PrintType() {return fPrintType;}
    int &Material() {return fMaterial;}

    /// Pushs back a term of the inverse incidence, i.e., an element which
    /// contains the node
    /// @param int element
    void pushNeighborElement(int el) {
        fNeighborElements.push_back(el);
    }


    /// Gets the number of elements which contains the node
    /// @return int number of elements which contains the node
    int getNumberOfNeighborElements(){
        return fNeighborElements.size();
    }

    /// Gets an specific member of the inverse incidence
    /// @param int index @return int element of the inverse incidence
    int64_t &getNeighborElement(int i){
        return fNeighborElements[i];
    }

    void clearInverseIncidence(){
        fNeighborElements.clear();
        fNeighborElements.shrink_to_fit();
    }


    void sortEraseNeighborElements(){
        std::sort(fNeighborElements.begin(), fNeighborElements.end());
        fNeighborElements.erase(std::unique(fNeighborElements.begin(),fNeighborElements.end()), fNeighborElements.end());
    }

    virtual int Dimension() = 0;

    VecDouble GetGeometricCenter();


    virtual int NCornerNodes() = 0;
    virtual int NEdges() = 0;
    virtual int NFaces() = 0;
    virtual int NVolumes() = 0;

    virtual VecDouble NodeCoord(int inode) = 0;

    virtual ElementType Type() = 0;


    /// Compute and store the spatial jacobian matrix
    /// @param bounded_vector integration point adimensional coordinates
    virtual void ComputeJacobian() = 0;
    virtual void ComputeCurrentJacobian() = 0;
    virtual void ComputeJacobianSearch() = 0;

    /// Compute and store the shape function spatial derivatives
    /// @param bounded_vector integration point adimensional coordinates
    virtual void ComputeSpatialDerivatives() = 0;
    virtual void ComputeCurrentSpatialDerivatives() = 0;
    virtual void ComputeHighOrderSpatialDerivatives() = 0;

    virtual double getJacobian() = 0;
};

#endif