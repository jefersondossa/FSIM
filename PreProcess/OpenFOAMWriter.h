#ifndef OPENFOAM_WRITER_H
#define OPENFOAM_WRITER_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include "DataTypes.h"

// Struct to store vertex data
struct Vertex {
    double x, y, z;
};

// Struct to store face data (triangular faces assumed)
struct Face {
    int v1, v2, v3;
};


class OpenFOAMWriter {
protected:
    std::string fInputMsh;
    std::string fOutputObj;


public:
    //Default constructor and destructor
    OpenFOAMWriter() = default;
    ~OpenFOAMWriter() = default;

    //Constructor
    OpenFOAMWriter(std::string inputFile);


    //Read the msh file
    bool ParseMSH(const std::string &inputFile, std::vector<Vertex> &vertices, std::vector<Face> &faces);

    //Write the obj file
    bool WriteOBJ(const std::string &outputFile, const std::vector<Vertex> &vertices, const std::vector<Face> &faces);
    
    //Find Max and Min coordinates of the building
    void FindMaxMin(const std::vector<Vertex> &vertices, double &maxX, double &maxY, double &maxZ, double &minX, double &minY, double &minZ);

    bool WriteBlockMeshDict(double dInlet, double dOutlet, double cellSize);

    bool WriteInitialConditions(VecDouble &internalField, double puniform, double nutuniform, double kuniform, double omegauniform);
    bool UInitial(VecDouble &internalField);
    bool pInitial(double puniform);
    bool nutInitial(double nutuniform);
    bool kInitial(double kuniform);
    bool omegaInitial(double omegauniform);

    bool WriteConstant();
    
};

#endif