#include "Connect.h"

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

Connect::Connect(int nstate, int nshape, int order, int64_t index){
    fNState = nstate;
    fNShape = nshape*fNState;
    fOrder = order;
    fIndex = index; 

    fSolution.resize(fNShape);
    fSolution.setZero();
    fPrevSolution.resize(fNShape);
    fPrevSolution.setZero();
    // fDTimeSolution.resize(fNShape);
    // fDDTimeSolution.resize(fNShape);
}

void Connect::AllocateTimeDerivatives(){
    fDTimeSolution.resize(fNShape);
    fDTimeSolution.setZero();
    fDDTimeSolution.resize(fNShape);
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