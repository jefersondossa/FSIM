#include "Node.h"


//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

Node::Node(VecDouble &coor, int64_t index){
    fIndex = index;
    fCoord.resize(coor.size());

    fCoord = coor;
    // fPrevCoord = coor;
    fInitialCoord = coor;

    fInverseIncidence.clear();
}

void Node::incrementCoordinate(int dir, double u){
    fCoord[dir] += u;
};

void Node::setPreviousCoordinates(int dir, double u){
    fPrevCoord[dir] = u;
};

void Node::setCoordinates(VecDouble &coor){
    for (int i=0; i<coor.size(); i++){
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
