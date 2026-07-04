#include "Connect.h"

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

Connect::Connect(int nstate, int nshape, int order, int index, int seqnum){
    fNState = nstate;
    fNShape = nshape;
    fOrder = order;
    fIndex = index;
    fSequenceNumber = seqnum;

    fSolution.resize(fNShape*fNState);
    fSolution.setZero();
    fPrevSolution.resize(fNShape*fNState);
    fPrevSolution.setZero();
    // fDTimeSolution.resize(fNShape);
    // fDDTimeSolution.resize(fNShape);
}

void Connect::AllocateTimeDerivatives(){
    fDTimeSolution.resize(fNShape*fNState);
    fDTimeSolution.setZero();
    fDDTimeSolution.resize(fNShape*fNState);
    fDDTimeSolution.setZero();
}

void Connect::SetSolution(int istate, REAL sol){
    // fPrevSolution[istate] = fSolution[istate];
    fSolution[istate] = sol;
}

void Connect::IncrementSolution(int istate, REAL sol){
    fSolution[istate] += sol;
}

void Connect::SetDSolutionDTime(int istate, REAL val){
    fDTimeSolution[istate] = val;
};

void Connect::SetDSolutionDDTime(int istate, REAL val){
    fDDTimeSolution[istate] = val;
};

void Connect::SetPreviousSolution(int istate, REAL val){
    fPrevSolution[istate] = val;
};


void Connect::setMeshVelocity(VecDouble &u){
    for (int i=3; i--; ){
        fMeshVelocity[i] = u[i];          
    };
};

void Connect::setMeshVelocityComponent(int dir,REAL u){
    fMeshVelocity[dir] = u;
} ;

void Connect::setWeightFunction(REAL val) {
    fPrevWeightFunction = fWeightFunction; 
    fWeightFunction = val;
};
