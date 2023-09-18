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

#include "Node.h"

/// Defines the gluing zone mesh
class Glue{
private:
    VecInt connect_;           //Velocity mesh connectivity 
    int index_;             //Element index
    int elemCorrespondent_;

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
    void setConnectivity(VecInt &connect){
        connect_.resize(connect.size());
        connect_ = connect;
    };

    /// Gets the element connectivity
    /// @return gluing zone element connectivity
    VecInt &getConnectivity(){return connect_;};
    
};


//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------


#endif

