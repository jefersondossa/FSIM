#ifndef VTUGENERATOR_H
#define VTUGENERATOR_H

#include<fstream>
#include<iostream>
#include<string>
#include "Element.h"
#include "GraphMesh.h"

namespace VTUGenerator{
    struct CustomCellField {
        std::string name;
        std::vector<double> data;
    };

    /// @brief Prints results in a .vtu file, readable in ParaView
    /// @param cmesh Computational mesh
    /// @param filename name of the output file
    /// @param scalnames scalar variables to be printed
    /// @param vecnames vector variables to be printed
    /// @param step (optional) step, for incremental and transient analyses 
    void PrintResults(CompMesh *cmesh, std::string filename, std::vector<std::string> &scalnames, std::vector<std::string> &vecnames, const std::vector<CustomCellField>& custom_cel_fields = {}, int step = 0);
    void PrintResults(Arlequin *arl, std::string filename);
    void PrintResultsGraph(CompMesh *cmesh, std::string filename, std::vector<std::string> &scalnames, std::vector<std::string> &vecnames, int step = 0);

};

#endif