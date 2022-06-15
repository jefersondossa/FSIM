//------------------------------------------------------------------------------
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2017 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//--------------------------------FLUID BOUNDARY--------------------------------
//------------------------------------------------------------------------------

#ifndef BOUNDARY_H
#define BOUNDARY_H

/// Defines the fluid boundary object and its properties

template<int DIM, int DEG>
class Boundary{

public: 
    /// Boundary element constructor
    /// @param BoundConnect boundary element connectivity
    /// @param int boundary element index
    /// @param int boundary element constrain type 
    /// @param int boundary element constrain value @see Node::setConstrains()
    Boundary(int *connec, int index, int *constrain, double *values, int gr){
        
        for(int i = 0; i<nBdNodes; i++) connectB_[i] = connec[i];
        index_ = index;
        group_ = gr;

        for(int i = 0; i<DIM; i++){
            constrainType[i] = constrain[i];
            constrainValue[i] = values[i];
        };

        element_ = -1;
        elementSide_ = 1000;
    };     

    /// Returns the boundary element constrain component type
    /// @param int direction @return boundary element constrain component type
    int getConstrain(int dir){return constrainType[dir];}

    /// Returns the boundary element constrain component value
    /// @param int direction @return boundary element constrain component value
    double getConstrainValue(int dir){return constrainValue[dir];}

    /// Returns the boundary element connectivity
    /// @return boundary element connectivity
    int* getBoundaryConnectivity(){return connectB_;}
    void setBoundaryConnectivity(int* connec){for(int i = 0; i<3*(DIM-1); i++) connectB_[i] = connec[i];}

    /// Sets the boundary element group
    /// @param int boundary element group
    void setBoundaryGroup(int gr){group_ = gr;}

    /// Gets the boundary element group
    /// @return boundary element group
    int getBoundaryGroup(){return group_;};

    /// Sets the fluid element correspondence
    /// @param int fluid element correspondence
    void setElement(int el){element_ = el;}

    /// Sets the fluid element correspondence side at the boundary
    /// @param int fluid element correspondence side
    void setElementSide(int el){elementSide_ = el;}

    /// Gets the fluid element correspondence
    /// @return fluid element correspondence
    int getElement(){return element_;}

    /// Gets the fluid element correspondence side at the boundary
    /// @return fluid element correspondence side
    int getElementSide(){return elementSide_;}

private:
    const int nBdNodes = 3*(1-DEG)+DIM*(2*DEG-1);
    int       connectB_[3*(1-DEG)+DIM*(2*DEG-1)];//Boundary element connectivity
    int          index_;            //Boundary element index
    int          constrainType[DIM];  //Element type of constrain
    double       constrainValue[DIM]; //Element constrain value
    int          element_;          //Fluid Element
    int          elementSide_;      //Fluid Element Side
    int          group_;            //Element boundary group
};






















#endif
