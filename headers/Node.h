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


#include "DataTypes.h"

/// Defines the node object and stores all nodal variables information

// template<int spaceDim>

class Node{
public:
    static const int spaceDim = 2;

private:
    int fNStateVariables;
    VecDouble fSolution;
    VecDouble fPrevSolution;
    VecDouble fCoord;
    VecDouble fPrevCoord;
    int fIndex, fDimension;

    //Main variables
    VecDouble        coord_;                  //Nodal coordinate vector
    VecDouble        initialCoord_;           //Initial nodal coordinate vector
    VecDouble        coordUpdated_;           //Updated nodal coordinate vector
    int              index_;                       //Node index
    VecDouble        previousCoord_;          //Previous nodal coordinate vector
    VecDouble        nNodal_;                 //Nodal normal vector
    
    //Fluid
    VecInt           constrainType;        //Constrain direction
    VecDouble        constrainValue;       //Nodal prescribed velocity

    VecInt           constrainTypeLaplace; //Constrain direction Laplace
    VecDouble        constrainValueLaplace;//Nodal prescribed value

    VecDouble        velocity_;               //Nodal velocity
    VecDouble        previousVelocity_;       //Previous time step velocity

    VecDouble        acceleration_;           //Nodal acceleration
    VecDouble        previousAcceleration_;   //Previous time step acceleration

    double           pressure_;               //Nodal pressure 

    double           vorticity_;              //Vorticity

    VecDouble        meshVelocity_;           //Nodal mesh velocity
    VecDouble        previousMeshVelocity_;   //Previous time step mesh velocity
    VecDouble        meshAcceleration_;

    //Arlequin
    int              elemCorresp;             //Element correspondence
    VecDouble        xsiCorresp;              //Adimensional coord correspond

    double           presArlequin_;           //Glue zone pressure    
    VecDouble        velArlequin_;            //Glue zone velocity


    VecDouble        lagMultiplier_;          //Nodal Lagrange Multiplier value
    double           weightFunction_;         //Nodal Energy Weight Function
    double           previousWeightFunction_;
    double           distGlueZone;            //Signaled distance to glue zone

    std::vector<int> invIncidence;
    
public:
    Node(VecDouble &coor, int index, int nState){
        fNStateVariables = nState;
        fIndex = index;
        fDimension = coor.size();

        fSolution.resize(fNStateVariables);
        fPrevSolution.resize(fNStateVariables);
        fCoord.resize(fDimension);
        fPrevCoord.resize(fDimension);
        fPrevSolution.setZero();
        fSolution.setZero();

        fCoord = coor;
        fPrevCoord = coor;

        index_ = index; 

        pressure_ = 0.;          
        vorticity_ = 0.;
        elemCorresp = 0;   
        presArlequin_=0.;
        weightFunction_ = 0.;   previousWeightFunction_ = 0.;
        distGlueZone = 0.;
        
        invIncidence.clear();

        coordUpdated_.resize(spaceDim);
        initialCoord_.resize(spaceDim);
        constrainValue.resize(spaceDim);
        constrainValueLaplace.resize(spaceDim);
        constrainType.resize(spaceDim);
        constrainTypeLaplace.resize(spaceDim);
        meshVelocity_.resize(spaceDim);
        previousMeshVelocity_.resize(spaceDim);
        meshAcceleration_.resize(spaceDim);
        xsiCorresp.resize(spaceDim);
        velArlequin_.resize(spaceDim);
        lagMultiplier_.resize(spaceDim);
        acceleration_.resize(spaceDim);
        previousAcceleration_.resize(spaceDim);
        velocity_.resize(spaceDim);
        previousVelocity_.resize(spaceDim);

        for (int i = 0; i < spaceDim; ++i){
            lagMultiplier_[i]= 0.;
            meshVelocity_[i] = 0.;
            previousMeshVelocity_[i] = 0.;
            meshAcceleration_[i] = 0.;
            velArlequin_[i] = 0.;
            velocity_[i] = 0.;
            previousVelocity_[i] = 0.;
            acceleration_[i] = 0.;
            previousAcceleration_[i] = 0.;
            xsiCorresp[i] = 0.;

            constrainTypeLaplace[i] = 0;
            constrainValueLaplace[i] = 0.;
            constrainType[i] = 0;
            constrainValue[i] = 0.;

            // previousCoord_[i] = coor[i];
            coordUpdated_[i] = coor[i];
            initialCoord_[i] = coor[i];
            // coord_[i] = coor[i];
        }
        nNodal_.resize(spaceDim);
    }

    ///Constructor - Defines a node with index and coordinates
    Node(VecDouble &coor, int index){
        index_ = index; 

        pressure_ = 0.;          
        vorticity_ = 0.;
        elemCorresp = 0;   
        presArlequin_=0.;
        weightFunction_ = 0.;   previousWeightFunction_ = 0.;
        distGlueZone = 0.;
        
        invIncidence.clear();

        coord_.resize(3);
        previousCoord_.resize(spaceDim);
        coordUpdated_.resize(spaceDim);
        initialCoord_.resize(spaceDim);
        constrainValue.resize(spaceDim);
        constrainValueLaplace.resize(spaceDim);
        constrainType.resize(spaceDim);
        constrainTypeLaplace.resize(spaceDim);
        meshVelocity_.resize(spaceDim);
        previousMeshVelocity_.resize(spaceDim);
        meshAcceleration_.resize(spaceDim);
        xsiCorresp.resize(spaceDim);
        velArlequin_.resize(spaceDim);
        lagMultiplier_.resize(spaceDim);
        acceleration_.resize(spaceDim);
        previousAcceleration_.resize(spaceDim);
        velocity_.resize(spaceDim);
        previousVelocity_.resize(spaceDim);

        for (int i = 0; i < spaceDim; ++i){
            lagMultiplier_[i]= 0.;
            meshVelocity_[i] = 0.;
            previousMeshVelocity_[i] = 0.;
            meshAcceleration_[i] = 0.;
            velArlequin_[i] = 0.;
            velocity_[i] = 0.;
            previousVelocity_[i] = 0.;
            acceleration_[i] = 0.;
            previousAcceleration_[i] = 0.;
            xsiCorresp[i] = 0.;

            constrainTypeLaplace[i] = 0;
            constrainValueLaplace[i] = 0.;
            constrainType[i] = 0;
            constrainValue[i] = 0.;

            previousCoord_[i] = coor[i];
            coordUpdated_[i] = coor[i];
            initialCoord_[i] = coor[i];
            coord_[i] = coor[i];
        }
        nNodal_.resize(spaceDim);

    };

    int &GetNStateVariables() {return fNStateVariables;}
    void SetSolution(int istate, double sol){
        fPrevSolution[istate] = fSolution[istate];
        fSolution[istate] = sol;
    }
    void IncrementSolution(int istate, double sol){
        fSolution[istate] += sol;
    }
    double GetSolution(int istate){
        return fSolution[istate];
    }

    /// Clear all node object variables
    void clearVariables();

    /// Returns the node coordinate vector
    /// @return node coordinate vector
    VecDouble &getCoordinates() {return fCoord;};

    /// Returns the node coordinate component value
    /// @return node coordinate component value
    double getCoordinateValue(int dir) const {return fCoord[dir];};
    double getPreviousCoordinateValue(int dir) const {return fPrevCoord[dir];};

    /// Returns the node initial coordinate vector
    /// @return node initial coordinate vector
    VecDouble &getInitialCoordinates() {return initialCoord_;};

    /// Returns the node coordinate vector at the previous time step
    /// @return node coordinate vector at the previous time step
    VecDouble &getPreviousCoordinates() {return previousCoord_;}

    /// Returns the node updated coordinate vector
    /// @return node coordinate updated vector
    VecDouble &getUpdatedCoordinates() {return coordUpdated_;};

    /// Increment the coordinate vector
    /// @param int direction @param double increment value
    void incrementCoordinate(int dir, double u){coord_[dir] += u;};

    /// Sets the previous coordinate vector
    /// @param int direction @param double value
    void setPreviousCoordinates(int dir, double u){previousCoord_[dir] = u;};

    /// Sets the node coordinate vector
    /// @param VecLocD Coordinate
    void setCoordinates(VecDouble &coor){for (int i=0; i<spaceDim; i++) coord_[i] = coor[i];};

    /// Sets the updated coordinate vector
    /// @param VecLocD Updated Coordinate
    void setUpdatedCoordinates(VecDouble &coor){for (int i=0; i<spaceDim; i++) coordUpdated_[i] = coor[i];};

    /// Updates node coordinate vector
    /// @param int direction @param double updated value
    void updateCoordinate(int dir, double val){coord_[dir] = val;};
    
    /// Sets nodal correspondence of overlapped mesh
    /// @param double element @param VecLocD adimensional coordinates
    void setNodalCorrespondence(double elem, VecDouble &xsi){
        elemCorresp = elem;
        for (int i=0; i<spaceDim; i++) xsiCorresp[i] = xsi[i];};

    /// Gets the nodal normal vector
    /// @return nodal normal vector
    VecDouble &getInnerNormal() {return nNodal_;};

    /// Sets nodal normal vector
    /// @param VecLocD nodal normal vector
    void setInnerNormal(VecDouble &n){for (int i=0; i<spaceDim; i++) nNodal_[i] = n[i];};

    /// Clears nodal normal vector
    /// @param VecLocD nodal normal vector
    void clearInnerNormal() {nNodal_.setZero();};

    /// Gets nodal correspondence of overlapped mesh - element
    /// @return Element correspondence of overlapped mesh
    int getNodalElemCorrespondence() const {return elemCorresp;}

    /// Gets nodal correspondence of overlapped mesh - adim. coordinate
    /// @return Adim. coordinate correspondence of overlapped mesh
    VecDouble &getNodalXsiCorrespondence() {return xsiCorresp;}

    /// Pushs back a term of the inverse incidence, i.e., an element which
    /// contains the node
    /// @param int element
    void pushInverseIncidence(int el) {invIncidence.push_back(el);}

    /// Gets the number of elements which contains the node
    /// @return int number of elements which contains the node
    int getNumberOfElements() const {return invIncidence.size();}

    /// Gets an specific member of the inverse incidence
    /// @param int index @return int element of the inverse incidence
    int getInverseIncidenceElement(int i) const {return invIncidence[i];}
    void clearInverseIncidence(){
        invIncidence.clear();
        invIncidence.shrink_to_fit();
    }
    
    //............................Velocity functions............................
    /// Sets the velocity vector
    /// @param double* velocity vector
    void setVelocity(VecDouble &u){for (int i=spaceDim; i--; ) velocity_[i] = u[i];};
    void setVelocityComponent(int dir, double val){velocity_[dir] = val;};

    void setPreviousVelocityComponent(int dir, double val){previousVelocity_[dir] = val;};


    /// Sets the velocity vector
    /// @param double* velocity vector
    void setVelocityGlobal(VecDouble &u);

    /// Sets the previous velocity vector
    /// @param double* previous time step velocity vector
    void setPreviousVelocity(VecDouble &u){for (int i=spaceDim; i--; ) previousVelocity_[i] = u[i];};

    /// Increment the velocity vector
    /// @param int direction @param double increment value
    void incrementVelocity(int dir, double u){velocity_[dir] += u;};

    /// Returns the node velocity vector
    /// @return node velocity vector
    double getVelocity(int dir) const {return velocity_[dir];}

    /// Returns the node previous time step velocity vector
    /// @return node previous time step velocity vector
    double getPreviousVelocity(int dir) const {return previousVelocity_[dir];}

    /// Sets the vorticity at the node
    /// @param double vorticity
    void setVorticity(double div) {vorticity_ = div;}

    /// Returns the nodal vorticity
    /// @return node vorticity
    double getVorticity() const {return vorticity_;}

    /// Clears the nodal vorticity
    void clearVorticity() {vorticity_ = 0.;}

    /// Increment the value of the nodal vorticity
    /// @param double increment 
    void incrementVorticity(double val) {vorticity_ += val/invIncidence.size();}

    //..........................Acceleration functions..........................
    /// Sets the acceleration vector
    /// @param double* acceleration vector
    void setAcceleration(VecDouble &u){for (int i=spaceDim; i--; ) acceleration_[i] = u[i];};
    void setAccelerationComponent(int dir, double val){acceleration_[dir] = val;};


    /// Sets the acceleration vector
    /// @param double* acceleration vector
    void setAccelerationGlobal(VecDouble &u);

    /// Increment the velocity vector
    /// @param int direction @param double increment value
    void incrementAcceleration(int dir, double u){acceleration_[dir] += u;};

    /// Sets the previous time step acceleration vector
    /// @param double* previous time step acceleration vector
    void setPreviousAcceleration(VecDouble &u){for (int i=spaceDim; i--; ) previousAcceleration_[i] = u[i];};
    void setPreviousAccelerationComponent(int dir, double val){previousAcceleration_[dir] = val;};

    /// Gets the acceleration vector
    /// @return acceleration vector
    double getAcceleration(int dir) const {return acceleration_[dir];}

    /// Gets the previous time step acceleration vector
    /// @return previous time step acceleration vector
    double getPreviousAcceleration(int dir) const {return previousAcceleration_[dir];}

    //............................Pressure functions............................
    /// Sets the nodal pressure
    /// @param double pressure
    void setPressure(double p){pressure_ = p;};

    /// Increments the nodal pressure
    /// @param double pressure increment
    void incrementPressure(double p){pressure_ += p; };

    /// Gets the nodal pressure value
    /// @return nodal pressure value
    double getPressure() const {return pressure_;};

    //.........................Mesh Velocity functions..........................
    /// Sets the node mesh velocity
    /// @param double* mesh velocity
    void setMeshVelocity(VecDouble &u){
        for (int i=spaceDim; i--; ){
            previousMeshVelocity_[i] = meshVelocity_[i];
            meshVelocity_[i] = u[i];          
        };
    };
    void setMeshVelocityComponent(int dir,double u){meshVelocity_[dir] = u;} ;
    void setMeshAccelerationComponent(int dir,double u){meshAcceleration_[dir] = u;} ;

    /// Sets the previous time step mesh velocity
    /// @param int direction @param double previous time step mesh velocity valu
    void setPreviousMeshVelocity(int dir, double u){previousMeshVelocity_[dir] = u;};

    /// Gets the node mesh velocity
    /// @param int direction @return mesh velocity component
    double getMeshVelocity(int dir) const {return meshVelocity_[dir];}
    double getMeshAcceleration(int dir) const {return meshAcceleration_[dir];}

    /// Gets the previous time step mesh velocity
    /// @param int direction @return previous time step mesh velocity component
    double getPreviousMeshVelocity(int dir) {return previousMeshVelocity_[dir];}
    void setPreviousMeshVelocityComponent(int dir, double val){previousMeshVelocity_[dir] = val;}


    //...........................Constrains functions...........................

    /// Sets all node constrains     
    /// @param int direction 
    /// @param int type: 0 - free, 1 - constrained, 2 - glue zone, 
    /// 3 - fluid-structure interface @param double constrain value
    void setConstrains(int dir, int type, double value){
        constrainType[dir] = type;
        constrainValue[dir] = value;
        velocity_[dir] = value; 
        // fSolution[0] = value;
        // previousVelocity_(dir) = value;
    };
    void SetBoundaryCondition(int dir, int type, double value){
        constrainType[dir] = type;
        constrainValue[dir] = value;
        fSolution[dir] = value;
        // previousVelocity_(dir) = value;
    };

    void setConstrainValue(int dir, double value){
        constrainValue[dir] = value;
    }

    /// Gets node constrain type
    /// @return constrain type
    int getConstrains(int dir) const {return constrainType[dir];}

    /// Gets node constrain value
    /// @return constrain value
    double getConstrainValue(int dir) const {return constrainValue[dir];}

    /// Sets constrains for solving the mesh moving problem
    /// @param int direction 
    /// @param int constrain type: 0 - free, 1 - constrained
    /// @param double constrain value
    void setConstrainsLaplace(int dir, int type, double value){
        constrainTypeLaplace[dir] = type;
        constrainValueLaplace[dir] = value;
    };

    /// Gets constrains of mesh moving problem
    /// @return constrain type
    int getConstrainsLaplace(int dir) const {return constrainTypeLaplace[dir];}

    /// Gets constrain value of mesh moving problem
    /// return constrain value
    double getConstrainValueLaplace(int dir) const {return constrainValueLaplace[dir];}

    //............................Arlequin functions............................
    /// Sets Lagrange Multiplier value
    /// @param int direction @param double component value
    void setLagrangeMultiplier(int dir, double lMult){
        lagMultiplier_[dir] = lMult;};

    /// Increment the Lagrange Multiplier vector
    /// @param int direction @param double increment value
    void incrementLagrangeMultiplier(int dir, double lMult){
        lagMultiplier_[dir] += lMult;};

    /// Gets Lagrange Multiplier component value
    /// @param int direction @return component value
    double getLagrangeMultiplier(int dir) const {return lagMultiplier_[dir];};

    /// Sets the nodal energy weight function value
    /// @param double weight function value
    void setWeightFunction(double val) {previousWeightFunction_ = weightFunction_; weightFunction_ = val;};

    /// Gets the nodal energy weight function value
    /// @return weight function value
    double getWeightFunction() const {return weightFunction_;};
    double getPreviousWeightFunction() const {return previousWeightFunction_;};

    /// Sets the interpolated Arlequin pressure 
    /// @param double interpolated pressure value
    void setPressureArlequin(double p) {presArlequin_ = p;};

    /// Gets interpolated Arlequin pressure
    /// @return interpolated Arlequin pressure
    double getPressureArlequin() const {return presArlequin_;};

    /// Sets the interpolated Arlequin velocity
    /// @param int direction @param double interpolated velocity component value
    void setVelocityArlequin(int dir, double v) {velArlequin_[dir] = v;};

    /// Gets the interpolated Arlequin velocity component
    /// @param int direction @return interpolated Arlequin velocity component
    double getVelocityArlequin(int dir) const {return velArlequin_[dir];};

    //.......................Signaled distance functions........................
    /// Sets the Signaled distance function value
    /// @param double signaled distance value
    void setDistFunction(double dist) {distGlueZone = dist;};

    /// Gets the Signaled distance function value
    /// @return signaled distance value
    double getDistFunction() const {return distGlueZone;};

};


#endif

