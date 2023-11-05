#ifndef VTUGENERATOR_H
#define VTUGENERATOR_H

#include<fstream>
#include<iostream>
#include<string>
#include "Element.h"
#include "CompMesh.h"

namespace VTUGenerator{

    void PrintResults(CompMesh *cmesh, std::string filename, int step = 0);

};

#endif