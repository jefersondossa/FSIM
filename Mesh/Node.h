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
    /// @brief Number of state variables
    int fNStateVariables;
    /// @brief Solution
    VecDouble fSolution;
    /// @brief Solution for the previous time/load step
    VecDouble fPrevSolution;
    /// @brief Nodal coordinates
    VecDouble fCoord;
    /// @brief Previous time/load step nodal coordinates 
    VecDouble fPrevCoord;
    /// @brief Initial coordinates
    VecDouble fInitialCoord;
    /// @brief Solution first time derivative
    VecDouble fDTimeSolution;
    /// @brief Solution second time derivative
    VecDouble fDDTimeSolution;
    /// @brief Node index
    int64_t fIndex;
    /// @brief Node dimension
    int fDimension;
    /// @brief Nodal mesh velocity
    VecDouble        fMeshVelocity;
    /// @brief Nodal Arlequin weight function
    double           fWeightFunction;
    /// @brief Nodal Arlequin previous time step weight funtion
    double           fPrevWeightFunction;
    /// @brief Inverse incidence, i.e., elements connected by the node
    std::vector<int64_t> fInverseIncidence;

    bool fHasBC = false;
    
public:
    /// @brief Class constructor
    /// @param coor Nodal coordinates
    /// @param index Node index
    /// @param nState Number of state variables
    Node(VecDouble &coor, int64_t index, int nState = 1);

    /// @brief Default destructor
    ~Node() = default;

    void SetHasBC(){
        fHasBC = true;
    }

    bool HasBC(){
        return fHasBC;
    }

    /// @brief Sets the number of state variables
    /// @param nstate number of state variables
    void SetNStateVariables(int nstate);

    /// @brief Gets the number of state variables
    /// @return number of state variables
    int &GetNStateVariables() {
        return fNStateVariables;
    }

    /// @brief Sets the nodal solution
    /// @param istate state variable index
    /// @param sol solution
    void SetSolution(int istate, double sol);

    /// @brief Increments the solution of a given state variable
    /// @param istate state variable index
    /// @param sol increment
    void IncrementSolution(int istate, double sol);

    /// @brief Returns a vector with all nodal solutions
    /// @return Solution vector
    VecDouble &Solution() {
        return fSolution;
    }
    
    /// @brief Returns the previous time/load step solution 
    /// @return Previous time/load step solution vector
    VecDouble &PrevSolution() {
        return fPrevSolution;
    }

    /// @brief Returns the nodal solution of a given state variable
    /// @param istate state variable index
    /// @return solution
    double GetSolution(int istate) const {
        return fSolution[istate];
    }

    /// @brief Returns the nodal solution of a given state variable in the previous time/load step
    /// @param istate state variable index
    /// @return solution
    double GetPreviousSolution(int istate) const {
        return fPrevSolution[istate];
    }
    void SetPreviousSolution(int istate, double val);

    /// @brief Returns the solution first time derivative of a given state variable
    /// @param istate state variable index
    /// @return solution first time derivative
    double GetDSolutionDTime(int istate) const {
        return fDTimeSolution[istate];
    }
    VecDouble &SolutionDTime(){
        return fDTimeSolution;
    }
    void SetDSolutionDTime(int istate, double val);

    /// @brief Returns the solution second time derivative of a given state variable
    /// @param istate state variable index
    /// @return solution second time derivative
    double GetDSolutionDDTime(int istate) const {
        return fDDTimeSolution[istate];
    }
    VecDouble &SolutionDDTime(){
        return fDDTimeSolution;
    }
    void SetDSolutionDDTime(int istate, double val);

    void AllocateTimeDerivatives();

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
    
    //.........................Mesh Velocity functions..........................
    /// Sets the node mesh velocity
    /// @param double* mesh velocity
    void setMeshVelocity(VecDouble &u);
    void setMeshVelocityComponent(int dir,double u);

    /// Gets the node mesh velocity
    /// @param int direction @return mesh velocity component
    double getMeshVelocity(int dir) const {
        return fMeshVelocity[dir];
    }
   
    //............................Arlequin functions............................
    /// Sets the nodal energy weight function value
    /// @param double weight function value
    void setWeightFunction(double val);

    /// Gets the nodal energy weight function value
    /// @return weight function value
    double getWeightFunction() const {
        return fWeightFunction;
    };
    double getPreviousWeightFunction() const {
        return fPrevWeightFunction;
    };
};


#endif

