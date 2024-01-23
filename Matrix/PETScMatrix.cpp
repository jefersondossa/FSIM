#include "PETScMatrix.h"

#ifdef HAS_PETSC

PETScMatrix::PETScMatrix(int64_t rows, int64_t cols) : MatrixType(rows,cols){


}

PETScMatrix::~PETScMatrix(){
}


#endif