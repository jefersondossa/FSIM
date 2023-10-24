#ifndef INTPOINTDATA
#define INTPOINTDATA

#include "DataTypes.h"


class IntPointData {

public:
    MatrixDouble fA0;
    MatrixDouble fA0Inv;

    MatrixDouble fA1;
    MatrixDouble fA1Inv;

    double fJacA0;
    double fJacA1;

    double fWeight;

    VecDouble fX;
    VecDouble fAdimCoord;

    MatrixDouble fDPhiX0;
    MatrixDouble fDPhiX1;
};


#endif