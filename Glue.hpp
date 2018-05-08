//------------------------------------------------------------------------------
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2017 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------
 
//------------------------------------------------------------------------------
//----------------------------------GLUE ZONE-----------------------------------
//------------------------------------------------------------------------------

#ifndef GLUE_H
#define GLUE_H

#include "Node.hpp"

/// Defines the gluing zone mesh

template<int DIM>
class Glue{

public:
    //Defines de class Node
    typedef Node<DIM>                                           Nodes;
    //Defines the element mesh connectivity - velocity
    typedef ublas::bounded_vector<int, 4*DIM-2>                 Connectivity;
    //Defines a blas-type vector with dimension = DIM
    typedef ublas::bounded_vector<double, DIM>                  DimVector;
    //Defines a blas-type matrix with dimension = DIM x DIM
    typedef ublas::bounded_matrix<double, DIM, DIM>             DimMatrix;
    //Defines the vector which contains the element nodal coordinates
    typedef ublas::bounded_matrix<double, 4*DIM-2, DIM>         LocalNodes;
    //Defines the local vector type with dimension 15 for DIM=2 and 34 for DIM=3
    typedef ublas::bounded_vector<double, 22*DIM-26>            LocalVector;
    //Defines the local matrix type with dimension 15x15
    //for DIM=2 and 34x34 for DIM=3
    typedef ublas::bounded_matrix<double, 22*DIM-26, 22*DIM-26> LocalMatrix;

    typedef typename Nodes::VecLocD                             VecLoc;
    typedef ublas::bounded_vector<double, 4*DIM-2>              NodalValue;

private:
    QuadShapeFunction<DIM> shapeQuad; //Quadratic shape function
    std::vector<Nodes *>   nodes_;    //Velocity nodes
    Connectivity  connect_;           //Velocity mesh connectivity 
    int           index_;             //Element index
    LocalNodes    localNodes_;        //Nodal coordinates - velocity
 
    LocalMatrix   jacobianNRMatrix;   //Newton's method jacobian
    LocalVector   rhsVector;          //RHS vector of Newton's method
    int elemCorrespondent_;
 
    static const double k1;
    static const double k2;
     
    //Second derivatives of velocity shape functions
    typename QuadShapeFunction<DIM>::ValueDDeriv ddphi_dx; 
    //First derivatives of velocity shape functions
    typename QuadShapeFunction<DIM>::ValueDeriv  dphi_dx, dphiL_dx;
    //Values of velocity shape functins
    typename QuadShapeFunction<DIM>::Values      phi_;     


public:
    /// Glue Zone Element constructor
    /// @param int glue zone element index 
    /// @param int fluid element correspondent
    Glue(int index, int elemCorrespondent){
        index_ = index;
        elemCorrespondent_ = elemCorrespondent;
    };

    /// Returns the fluid element correspondence
    /// @return fluid element correspondence
    int getElemCorrespondent(){return elemCorrespondent_;};

    /// Sets the element connectivity
    /// @param Connectivity element connectivity
    void setConnectivity(Connectivity& connect){connect_ = connect;};

    /// Sets the vector of nodes of the gluing zone mesh
    /// @param vector<Nodes> gluing zone mesh nodes
    void setNodes(std::vector<Nodes *> nodes){nodes_ = nodes;};

    /// Gets the element connectivity
    /// @return gluing zone element connectivity
    Connectivity getConnectivity(){return connect_;};
    
};

template<>
double const Glue<2>::k1 = 1.0;
template<>
double const Glue<2>::k2 = 0.0;



//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------


#endif

