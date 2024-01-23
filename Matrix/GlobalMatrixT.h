#ifndef STRUCMATRIX
#define STRUCMATRIX
#include "SparseEigen.h"
#include "PETScDataTypes.h"
#include "GlobalMatrix.h"

template<class MatType>
class GlobalMatrixT : public GlobalMatrix
{
protected:
    MatType *fGlobalMatrix;

public:
    GlobalMatrixT(int64_t rows, int64_t cols);
    
    ~GlobalMatrixT();
};



#endif