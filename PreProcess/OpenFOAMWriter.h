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
    double fMaxX, fMaxY, fMaxZ;
    double fMinX, fMinY, fMinZ;
    double fDInlet, fDOutlet;
    double fRefProportion;
    int fNSubdomains;
    std::string fClass;//NBR 6123 - Building Class
    int fNRefinements;

public:
    //Default constructor and destructor
    OpenFOAMWriter() = default;
    ~OpenFOAMWriter() = default;

    //Constructor
    OpenFOAMWriter(std::string inputFile);

    void ClearAllFiles();

    void StartFromPreviousResults(double dt, double endTime, double writeInterval);

    //Read the msh file
    bool ParseMSH(const std::string &inputFile, std::vector<Vertex> &vertices, std::vector<Face> &faces);

    //Write the obj file
    bool WriteOBJ(const std::string &outputFile, const std::vector<Vertex> &vertices, const std::vector<Face> &faces);
    
    //Find Max and Min coordinates of the building
    void FindMaxMin(const std::vector<Vertex> &vertices, double &maxX, double &maxY, double &maxZ, double &minX, double &minY, double &minZ);

    //folder 0
    bool WriteInitialConditions(VecDouble &internalField, double puniform, double nutuniform, double kuniform, double omegauniform);
    bool UInitial(VecDouble &internalField);
    bool pInitial(double puniform);
    bool nutInitial(double nutuniform);
    bool kInitial(double kuniform);
    bool omegaInitial(double omegauniform);
    void VelocityNBR6123();

    //folder constant
    bool WriteConstant(double density, double viscosity);
    bool WriteMomentumTransport();
    bool WritePhysicalProperties(double density, double viscosity); 

    //folder system
    bool WriteSystem(double dInlet, double dOutlet, double cellSizeX, double cellSizeY, double cellSizeZ, double refproportion, int nref, double dt, double endTime, double writeInterval);
    bool WriteBlockMeshDict(double dInlet, double dOutlet, double cellSizeX, double cellSizeY, double cellSizeZ);
    bool WriteSurfaceFeaturesDict();
    bool WriteControlDict(double dt, double endTime, double writeInterval);
    bool WriteFvSchemes();
    bool WriteFvSolution();
    bool WriteSnappyHexMeshDict();
    bool WriteMeshQualityDict();
    bool WriteDecomposeParDict();

    int GetNSubdomains() {return fNSubdomains;};

};

#endif