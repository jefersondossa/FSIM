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

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric;

template<int DIM>

class Node{
public:
    //store doubles to the dimension size
    typedef ublas::bounded_vector<double,DIM>        VecLocD;
    //poir integer double (help to set constraints)
    typedef std::vector<std::pair<int, double> >     VecConstrD;
    //help to build loops in order to enforce bcs
    typedef VecConstrD::iterator                     VecConstrDIt;

    static const int spaceDim = DIM;

private:
    //Main variables
    VecLocD          coord_;                  //Nodal coordinate vector
    VecLocD          initialCoord_;           //Initial nodal coordinate vector
    VecLocD          coordUpdated_;           //Updated nodal coordinate vector
    int              index_;                  //Node index
    VecLocD          previousCoord_;          //Previous nodal coordinate vector


    //Fluid
    int              constrainType[3];        //Constrain direction
    double           constrainValue[3];       //Nodal prescribed velocity

    int              constrainTypeLaplace[3]; //Constrain direction Laplace
    double           constrainValueLaplace[3];//Nodal prescribed value

    VecLocD          velocity_;               //Nodal velocity
    VecLocD          previousVelocity_;       //Previous time step velocity

    VecLocD          acceleration_;           //Nodal acceleration
    VecLocD          previousAcceleration_;   //Previous time step acceleration

    double           pressure_;               //Nodal pressure 
    double           previousPressure_;       //Previous time step pressure

    double           divergent_;              //Velocity divergent

    VecLocD          meshVelocity_;           //Nodal mesh velocity
    VecLocD          previousMeshVelocity_;   //Previous time step mesh velocity

    //Arlequin
    int              elemCorresp;             //Element correspondence
    VecLocD          xsiCorresp;              //Adimensional coord correspond

    double           presArlequin_;           //Glue zone pressure    
    VecLocD          velArlequin_;            //Glue zone velocity

    VecLocD          lagMultiplier_;          //Nodal Lagrange Multiplier value
    double           weightFunction_;         //Nodal Energy Weight Function
    double           distGlueZone;            //Signaled distance to glue zone

    //Potential problem
    VecLocD          gradient_;               //Potential gradient
    double           potential_;              //Potential value
    
public:
    //!Constructor - Defines a node with index and coordinates
    Node(VecLocD& coor, int index){
        coord_ = coor;
        initialCoord_ = coor;
        coordUpdated_ = coor;
        previousCoord_ = coor;
        index_ = index; 

        constrainType[0] = 0;    constrainType[1] = 0;    constrainType[2] = 0;
        constrainValue[0] = 0;   constrainValue[1] = 0;   constrainValue[2] = 0;
        pressure_ = 0.;          previousPressure_ = 0.;  divergent_ = 0.;
        elemCorresp = 0;         
        velocity_.clear();   previousVelocity_.clear();   acceleration_.clear();
        previousAcceleration_.clear(); xsiCorresp.clear(); gradient_.clear();
        potential_ = 0.;     lagMultiplier_.clear();      weightFunction_ = 0.;
        distGlueZone = 0.;

        constrainTypeLaplace[0] = 0;    constrainTypeLaplace[1] = 0;
        constrainTypeLaplace[2] = 0;    constrainValueLaplace[0] = 0;
        constrainValueLaplace[1] = 0;   constrainValueLaplace[2] = 0;

    };

    //Clear all nodal variables
    void clearVariables();

    //Returns the nodal coordinate vector
    VecLocD getCoordinates() {return coord_;};
    VecLocD getInitialCoordinates() {return initialCoord_;};
    VecLocD getPreviousCoordinates() {return previousCoord_;;}
    VecLocD getUpdatedCoordinates() {return coordUpdated_;};

    //Increment coordinate vector
    void incrementCoordinate(int dir, double u);

    //Sets previous coordinate vector
    void setPreviousCoordinates(int dir, double u);

    //Set the nodal coordinates
    void setCoordinates(VecLocD& coor){
        coord_ = coor;};
    void setUpdatedCoordinates(VecLocD& coor){
        coordUpdated_ = coor;};

    //Updates nodal coordinates
    void updateCoordinate(int dir, double val){
        coord_(dir) = val;};
    
    //Sets nodal correspondence of overlapped mesh
    void setNodalCorrespondence(double elem, const VecLocD& xsi){
        elemCorresp = elem;
        xsiCorresp = xsi;
    };

    //Gets nodal correspondence of overlapped mesh
    int getNodalElemCorrespondence() {return elemCorresp;}
    VecLocD getNodalXsiCorrespondence() {return xsiCorresp;}
    
    //Velocity functions
    void setVelocity(double *u);
    void setPreviousVelocity(double *u);
    void incrementVelocity(int dir, double u);
    double getVelocity(int dir) {return velocity_(dir);}
    double getPreviousVelocity(int dir) {return previousVelocity_(dir);}
    void setVelocityDivergent(double div) {divergent_ = div;}
    double getVelocityDivergent() {return divergent_;}

    //Acceleration functions
    void setAcceleration(double *u);
    void setPreviousAcceleration(double *u);
    double getAcceleration(int dir) {return acceleration_(dir);}
    double getPreviousAcceleration(int dir) {return previousAcceleration_(dir);}

    //Pressure functions
    void setPressure(double p);
    void incrementPressure(double p);
    double getPressure() {return pressure_;};
    void setPreviousPressure(double p);
    double getPreviousPressure() {return previousPressure_;}

    //Mesh Velocity functions
    void setMeshVelocity(double *u);
    void setPreviousMeshVelocity(int dir, double u);
    double getMeshVelocity(int dir) {return meshVelocity_(dir);}
    double getPreviousMeshVelocity(int dir) {return previousMeshVelocity_(dir);}

    //Constrain functions
    void setConstrains(int dir, int type, double value){
        constrainType[dir] = type;
        constrainValue[dir] = value;
        velocity_(dir) = value;
        previousVelocity_(dir) = value;
    };
    int getConstrains(int dir) {return constrainType[dir];}
    double getConstrainValue(int dir) {return constrainValue[dir];}

    void setConstrainsLaplace(int dir, int type, double value){
        constrainTypeLaplace[dir] = type;
        constrainValueLaplace[dir] = value;
    };
    int getConstrainsLaplace(int dir) {return constrainTypeLaplace[dir];}
    double getConstrainValueLaplace(int dir){return constrainValueLaplace[dir];}

    //Arlequin functions
    void setLagrangeMultiplier(int dir, double lMult){
        lagMultiplier_(dir) = lMult;};
    void incrementLagrangeMultiplier(int dir, double lMult){
        lagMultiplier_(dir) += lMult;};
    double getLagrangeMultiplier(int dir) {return lagMultiplier_(dir);};
    void setWeightFunction(double val) {weightFunction_ = val;};
    double getWeightFunction() {return weightFunction_;};
    void setPressureArlequin(double p) {presArlequin_ = p;};
    double getPressureArlequin() {return presArlequin_;};
    void setVelocityArlequin(int dir, double v) {velArlequin_(dir) = v;};
    double getVelocityArlequin(int dir) {return velArlequin_(dir);};

    //Signaled distance
    void setDistFunction(double dist) {distGlueZone = dist;};
    double getDistFunction(){return distGlueZone;};


    //Potential problem functions
    void setGradientComponent(double val, int dir) {gradient_(dir) = val;};
    double getGradientComponent(int dir) {return gradient_(dir);};
    void setPotential(double p) {potential_ = p;};
    void incrementPotential(double p) {potential_ += p;};
    double getPotential() {return potential_;};

};

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-------------------------------CLEAR VARIABLES--------------------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::clearVariables(){
    pressure_ = 0.;          previousPressure_ = 0.;  divergent_ = 0.;
    elemCorresp = 0;

    if (constrainType[0] != 1){
        velocity_(0) = 0.;   
        previousVelocity_(0) = 0.;   
        acceleration_(0) = 0.;
        previousAcceleration_(0) = 0.;
    };

    if (constrainType[1] != 1){
        velocity_(1) = 0.;   
        previousVelocity_(1) = 0.;   
        acceleration_(1) = 0.;
        previousAcceleration_(1) = 0.;
    };
    
    xsiCorresp.clear();    gradient_.clear();
    potential_ = 0.;       lagMultiplier_.clear();      weightFunction_ = 0.;
    
    return;
};

template<>
void Node<3>::clearVariables(){

    pressure_ = 0.;          previousPressure_ = 0.;  divergent_ = 0.;
    elemCorresp = 0;

    if (constrainType[0] != 1){
        velocity_(0) = 0.;   
        previousVelocity_(0) = 0.;   
        acceleration_(0) = 0.;
        previousAcceleration_(0) = 0.;
    };

    if (constrainType[1] != 1){
        velocity_(1) = 0.;   
        previousVelocity_(1) = 0.;   
        acceleration_(1) = 0.;
        previousAcceleration_(1) = 0.;
    };

    if (constrainType[2] != 1){
        velocity_(2) = 0.;   
        previousVelocity_(2) = 0.;   
        acceleration_(2) = 0.;
        previousAcceleration_(2) = 0.;
    };
    
    xsiCorresp.clear();    gradient_.clear();
    potential_ = 0.;       lagMultiplier_.clear();     weightFunction_ = 0.;
    
    return;
};

//------------------------------------------------------------------------------
//----------------------------NODAL VELOCITY VALUES-----------------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::setVelocity(double *u){
    //Sets Velocity value
    velocity_(0) = u[0];
    velocity_(1) = u[1]; 
    return;
};

template<>
void Node<3>::setVelocity(double *u){
    //Sets Velocity value
    velocity_(0) = u[0];
    velocity_(1) = u[1]; 
    velocity_(2) = u[2]; 
    return;
};

//------------------------------------------------------------------------------
//-----------------------INCREMENT NODAL VELOCITY VALUES------------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::incrementVelocity(int dir, double u){
    //All element nodes
    velocity_(dir) += u;
    return;
};

template<>
void Node<3>::incrementVelocity(int dir, double u){
    //All element nodes
    velocity_(dir) += u;

    return;
};

//------------------------------------------------------------------------------
//---------------------INCREMENT NODAL COORDINATE VALUES------------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::incrementCoordinate(int dir, double u){
    //All element nodes
    coord_(dir) += u;
    return;
};

template<>
void Node<3>::incrementCoordinate(int dir, double u){
    //All element nodes
    coord_(dir) += u;

    return;
};

//------------------------------------------------------------------------------
//---------------------SETS PREVIOUS NODAL VELOCITY VALUES----------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::setPreviousCoordinates(int dir, double u){
    //All element nodes
    previousCoord_(dir) = u;
    return;
};

//------------------------------------------------------------------------------
//----------------------------NODAL PRESSURE VALUES-----------------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::setPressure(double p){
    //Only element nodes 0, 1 and 2  
    pressure_ = p; 
    return;
};

template<>
void Node<3>::setPressure(double p){
    //Only element nodes 0, 1, 2 and 3
    pressure_ = p; 
    return;
};

//------------------------------------------------------------------------------
//----------------------------NODAL PRESSURE VALUES-----------------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::setPreviousPressure(double p){
    //Only element nodes 0, 1 and 2  
    previousPressure_ = p; 
    return;
};

template<>
void Node<3>::setPreviousPressure(double p){
    //Only element nodes 0, 1, 2 and 3
    previousPressure_ = p; 
    return;
};

//------------------------------------------------------------------------------
//-----------------------INCREMENT NODAL PRESSURE VALUES------------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::incrementPressure(double p){
    //Only element nodes 0, 1 and 2  
    pressure_ += p; 
    return;
};

template<>
void Node<3>::incrementPressure(double p){
    //Only element nodes 0, 1, 2 and 3
    pressure_ += p; 
    return;
};

//------------------------------------------------------------------------------
//---------------------SETS PREVIOUS NODAL VELOCITY VALUES----------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::setPreviousVelocity(double *u){
    //All element nodes
    previousVelocity_(0) = u[0];
    previousVelocity_(1) = u[1]; 
    return;
};

template<>
void Node<3>::setPreviousVelocity(double *u){
    //All element nodes
    previousVelocity_(0) = u[0];
    previousVelocity_(1) = u[1]; 
    previousVelocity_(2) = u[2]; 
    return;
};

//------------------------------------------------------------------------------
//--------------------------NODAL ACCELERATION VALUES---------------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::setAcceleration(double *u){
    //All element nodes
    acceleration_(0) = u[0];
    acceleration_(1) = u[1]; 
    return;
};

template<>
void Node<3>::setAcceleration(double *u){
    //All element nodes
    acceleration_(0) = u[0];
    acceleration_(1) = u[1]; 
    acceleration_(2) = u[2]; 
    return;
};

//------------------------------------------------------------------------------
//-------------------SETS PREVIOUS NODAL ACCELERATION VALUES--------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::setPreviousAcceleration(double *u){
    //All element nodes
    previousAcceleration_(0) = u[0];
    previousAcceleration_(1) = u[1]; 
    return;
};

template<>
void Node<3>::setPreviousAcceleration(double *u){
    //All element nodes
    previousAcceleration_(0) = u[0];
    previousAcceleration_(1) = u[1]; 
    previousAcceleration_(2) = u[2]; 
    return;
};


//------------------------------------------------------------------------------
//----------------------------NODAL VELOCITY VALUES-----------------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::setMeshVelocity(double *u){
    //All element nodes
    previousMeshVelocity_(0) = meshVelocity_(0);
    previousMeshVelocity_(1) = meshVelocity_(1); 

    meshVelocity_(0) = u[0];
    meshVelocity_(1) = u[1]; 
    return;
};

template<>
void Node<3>::setMeshVelocity(double *u){
    //All element nodes
    meshVelocity_(0) = u[0];
    meshVelocity_(1) = u[1]; 
    meshVelocity_(2) = u[2]; 
    return;
};

//------------------------------------------------------------------------------
//----------------------------NODAL VELOCITY VALUES-----------------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::setPreviousMeshVelocity(int dir, double u){
    //All element nodes
    previousMeshVelocity_(dir) = u;

    return;
};


#endif
