#ifndef INTPOINTDATA
#define INTPOINTDATA

#include "DataTypes.h"


class IntPointData {

public:
    MatrixDouble fA0;
    MatrixDouble fA0Inv;
    MatrixDouble fAxes0;
    MatrixDouble fAxes1;

    MatrixDouble fA1;

    double fJacA0;
    double fJacA1;

    double fWeight;

    VecDouble fX;
    VecDouble fX1;
    VecDouble fAdimCoord;

    MatrixDouble fDPhiX0;
    MatrixDouble fDPhiX1;

    VecDouble    fPhi;
    MatrixDouble fDPhi;

    VecDouble fWeightFunction;
    VecDouble fPrevWeightFunction;
    VecDouble fDistFunction;

    VecDouble fSol;
    MatrixDouble fDSolDx;

    bool fNeedsSol = false;
    bool fNeedsDSol = false;
    

};


#endif