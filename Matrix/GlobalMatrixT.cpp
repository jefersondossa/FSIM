#include "GlobalMatrixT.h"

template<class MatType>
GlobalMatrixT<MatType>::GlobalMatrixT(int64_t rows, int64_t cols){
    fGlobalMatrix = new MatType(rows,cols);

}

template<class MatType>
GlobalMatrixT<MatType>::~GlobalMatrixT()
{
}

#include "EigenSpMatrix.h"
template class GlobalMatrixT<EigenSpMatrix>;

#ifdef HAS_PETSC
#include "PETScMatrix.h"
template class GlobalMatrixT<PETScMatrix>;
#endif