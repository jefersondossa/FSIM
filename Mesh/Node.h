//------------------------------------------------------------------------------
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2017 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------------NODES-------------------------------------
//------------------------------------------------------------------------------

#ifndef NODE_H
#define NODE_H

// #include <boost/numeric/ublas/io.hpp>

// using namespace boost::numeric;

#include <math.h> 
#include <stdio.h> 
#include <vector>
#include <algorithm>
#include <functional>
#include <array>
#include <iomanip>
#include <cmath>
#include <limits>
#include <ctime>
#include <numeric>
#include <iostream>
#include<fstream>

#include "DataTypes.h"

/// Defines the node object and stores all nodal variables information
class Node{
private:
    /// @brief Nodal coordinates
    VecDouble fCoord;
    /// @brief Previous time/load step nodal coordinates 
    VecDouble fPrevCoord;
    /// @brief Initial coordinates
    VecDouble fInitialCoord;
    /// @brief Node index
    int64_t fIndex;
    /// @brief Inverse incidence, i.e., elements connected by the node
    std::vector<int64_t> fInverseIncidence;
    
public:
    /// @brief Class constructor
    /// @param coor Nodal coordinates
    /// @param index Node index
    Node(VecDouble &coor, int64_t index);

    /// @brief Default destructor
    ~Node() = default;

    /// @brief Returns the node index
    /// @return index
    int64_t &Index() {
        return fIndex;
    }

    /// Returns the node coordinate vector
    /// @return node coordinate vector
    VecDouble &getCoordinates() {
        return fCoord;
    };

    /// Returns the node coordinate component value
    /// @return node coordinate component value
    double getCoordinateValue(int dir) const {
        return fCoord[dir];
    };

    /// @brief Returns the node component value for a given direction
    /// @param dir direction
    /// @return previous time/load step coordinate
    double getPreviousCoordinateValue(int dir) const {
        return fPrevCoord[dir];
    };

    /// Returns the node initial coordinate vector
    /// @return node initial coordinate vector
    VecDouble &getInitialCoordinates() {
        return fInitialCoord;
    };

    /// Returns the node coordinate vector at the previous time step
    /// @return node coordinate vector at the previous time step
    VecDouble &getPreviousCoordinates() {
        return fPrevCoord;
    }

    /// Increment the coordinate vector
    /// @param int direction @param double increment value
    void incrementCoordinate(int dir, double u);

    /// Sets the previous coordinate vector
    /// @param int direction @param double value
    void setPreviousCoordinates(int dir, double u);

    /// Sets the node coordinate vector
    /// @param VecLocD Coordinate
    void setCoordinates(VecDouble &coor);

    /// Pushs back a term of the inverse incidence, i.e., an element which
    /// contains the node
    /// @param int element
    void pushInverseIncidence(int el);

    /// Gets the number of elements which contains the node
    /// @return int number of elements which contains the node
    int64_t getNumberOfElements() const {
        return fInverseIncidence.size();
    }

    /// Gets an specific member of the inverse incidence
    /// @param int index @return int element of the inverse incidence
    int getInverseIncidenceElement(int i) const {
        return fInverseIncidence[i];
    }

    /// @brief Clears the inverse incidence vector
    void clearInverseIncidence();
    
    
};


#endif

