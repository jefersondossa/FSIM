#include "Connect.h"

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

Connect::Connect(int nstate, int nshape, int order, int64_t index, int64_t seqnum){
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

void Connect::SetSolution(int istate, double sol){
    // fPrevSolution[istate] = fSolution[istate];
    fSolution[istate] = sol;
}

void Connect::IncrementSolution(int istate, double sol){
    fSolution[istate] += sol;
}

void Connect::SetDSolutionDTime(int istate, double val){
    fDTimeSolution[istate] = val;
};

void Connect::SetDSolutionDDTime(int istate, double val){
    fDDTimeSolution[istate] = val;
};

void Connect::SetPreviousSolution(int istate, double val){
    fPrevSolution[istate] = val;
};