#include "Node.h"


//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

Node::Node(VecDouble &coor, int64_t index){
    fIndex = index;
    fDimension = coor.size();
    fCoord.resize(fDimension);

    fCoord = coor;
    // fPrevCoord = coor;
    fInitialCoord = coor;

    fWeightFunction = 0.;   fPrevWeightFunction = 0.;
    fInverseIncidence.clear();
}

void Node::incrementCoordinate(int dir, double u){
    fCoord[dir] += u;
};

void Node::setPreviousCoordinates(int dir, double u){
    fPrevCoord[dir] = u;
};

void Node::setCoordinates(VecDouble &coor){
    for (int i=0; i<fDimension; i++){
        fCoord[i] = coor[i];
    } 
};

void Node::pushInverseIncidence(int el) {
    fInverseIncidence.push_back(el);
}

void Node::clearInverseIncidence(){
    fInverseIncidence.clear();
    fInverseIncidence.shrink_to_fit();
}

void Node::setMeshVelocity(VecDouble &u){
    for (int i=fDimension; i--; ){
        fMeshVelocity[i] = u[i];          
    };
};

void Node::setMeshVelocityComponent(int dir,double u){
    fMeshVelocity[dir] = u;
} ;

void Node::setWeightFunction(double val) {
    fPrevWeightFunction = fWeightFunction; 
    fWeightFunction = val;
};
