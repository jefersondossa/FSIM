#ifndef VTUGENERATOR_H
#define VTUGENERATOR_H

#include<fstream>
#include<iostream>
#include<string>
#include "Element.h"
#include "GraphMesh.h"

namespace VTUGenerator{

    void PrintResults(CompMesh *cmesh, std::string filename, std::vector<std::string> &scalnames, std::vector<std::string> &vecnames, int step = 0);

};

#endif