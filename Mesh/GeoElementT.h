#ifndef GeoElementT_h
#define GeoElementT_h

#include "GeoElement.h"
#include "IntPointData.h"
#include "Element.h"


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
    void ComputeJacobian(IntPointData &data) override;
    void ComputeCurrentJacobian(IntPointData &data, Element* compel) override;
    void ComputeJacobianSearch(IntPointData &data) override;

    


    //.................Element intersection and correspondence..................
    /// Gets the element intersection parameters 
    /// @param minimum coordinates @param maximum coordinates 
    /// @param minimum inner product @param maximum inner product
    /// @param side lenght
    void setIntersectionParameters(VecDouble &x, VecDouble &X) override ;


};

#endif