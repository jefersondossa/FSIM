#ifndef GeoElementT_h
#define GeoElementT_h

#include "GeoElement.h"

template<class geoshape>
class GeoElementT:public GeoElement{

protected:

public:
    GeoElementT();
    ~GeoElementT() = default;

    GeoElementT(int64_t index, VecInt &geonodes, GeoMesh* mesh, int &material);


    int Dimension() override{
        return geoshape::Dimension;
    }

    int NCornerNodes() override {return geoshape::NCornerNodes;}
    int NEdges() override {return geoshape::NEdges;}
    int NFaces() override {return geoshape::NFaces;}
    int NVolumes() override {return geoshape::NVolumes;}


    VecDouble NodeCoord(int inode) override{
        VecDouble xnode(3);
        MatrixDouble coords(3,geoshape::NShape);
        geoshape::getCoordinates(coords);
        for (int i = 0; i < 3; i++){
            xnode[i] = coords(i,inode);
        }
        return xnode;
    }

    ElementType Type() override{
        return geoshape::ElType;
    }

    //........................Element basic information.........................
    /// Compute and store the spatial jacobian matrix
    /// @param bounded_vector integration point adimensional coordinates
    void ComputeJacobian() override;
    void ComputeCurrentJacobian() override;
    void ComputeJacobianSearch() override;

    /// Compute and store the shape function spatial derivatives
    /// @param bounded_vector integration point adimensional coordinates
    void ComputeSpatialDerivatives() override;
    void ComputeCurrentSpatialDerivatives() override;
    void ComputeHighOrderSpatialDerivatives() override;

};

#endif