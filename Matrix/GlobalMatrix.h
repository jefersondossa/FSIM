#ifndef GLOBALMATRIX_H
#define GLOBALMATRIX_H
#include "PanicButton.h"

class GlobalMatrix
{
protected:

public:
    GlobalMatrix() = default;

    GlobalMatrix(int64_t rows, int64_t cols);

    ~GlobalMatrix();
};


#endif

