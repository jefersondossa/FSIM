//------------------------------------------------------------------------------
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2017 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-----------------------------------ELEMENT------------------------------------
//------------------------------------------------------------------------------

#ifndef GLUE_H
#define GLUE_H

#include "Node.hpp"

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
    //Element definition
    Glue(int index, int elemCorrespondent){
        index_ = index;
        elemCorrespondent_ = elemCorrespondent;
        
        // setLocalNodes();
        // getIntegPointCoordinates();
    };

    int getElemCorrespondent(){return elemCorrespondent_;};
    void setConnectivity(Connectivity& connect){connect_ = connect;};
    void setNodes(std::vector<Nodes *> nodes){nodes_ = nodes;};
    Connectivity getConnectivity(){return connect_;};
    
};

template<>
double const Glue<2>::k1 = 1.0;
template<>
double const Glue<2>::k2 = 0.0;



//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//----------------------SET ELEMENT INTERSECTION PARAMETERS---------------------
//------------------------------------------------------------------------------
// template<>
// void Element<2>::setIntersectionParameters(DimVector x, DimVector X, 
//                                            double *Dk, double *dk,
//                            std::vector<ublas::bounded_vector<double,2> > dii) {
//     xK = x; XK = X; 
//     dCk[0] = Dk[0];
//     dCk[1] = Dk[1]; 
//     dCk[2] = Dk[2];  
    
//     dck[0] = dk[0];
//     dck[1] = dk[1];
//     dck[2] = dk[2];   
    
//     di = dii;
    
//     return;
// };

#endif

