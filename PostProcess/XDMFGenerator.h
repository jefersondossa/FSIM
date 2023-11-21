#ifndef XDMFGENERATOR_H
#define XDMFGENERATOR_H

#include <string>
#include <iostream>

/// @brief Implements the results print in a XDMF file, using HDF5 library
namespace XDMFGenerator{
    //Need refactor
    void ReadInitialSolution(const std::string& inputPrev, const std::string& inputCurr);

};

#endif