#include "Node.h"


//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

Node::Node(VecDouble &coor, int64_t index, int nState){
    fNStateVariables = nState;
    fIndex = index;
    fDimension = coor.size();

    fSolution.resize(fNStateVariables);
    fPrevSolution.resize(fNStateVariables);
    fCoord.resize(fDimension);
    fPrevSolution.setZero();
    fSolution.setZero();

    fCoord = coor;
    fPrevCoord = coor;
    fInitialCoord = coor;

    fWeightFunction = 0.;   fPrevWeightFunction = 0.;
    
    fInverseIncidence.clear();
}


void Node::SetNStateVariables(int nstate){
    fNStateVariables = nstate;
    fSolution.resize(fNStateVariables);
    fDTimeSolution.resize(fNStateVariables);
    fDDTimeSolution.resize(fNStateVariables);
    fPrevSolution.resize(fNStateVariables);
}


void Node::SetSolution(int istate, double sol){
    fPrevSolution[istate] = fSolution[istate];
    fSolution[istate] = sol;
}

void Node::IncrementSolution(int istate, double sol){
    fSolution[istate] += sol;
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
