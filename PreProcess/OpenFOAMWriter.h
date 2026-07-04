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
    REAL x, y, z;
};

// Struct to store face data (triangular faces assumed)
struct Face {
    int v1, v2, v3;
};


class OpenFOAMWriter {
protected:
    std::string fInputMsh;
    std::string fOutputObj;
    REAL fMaxX, fMaxY, fMaxZ;
    REAL fMinX, fMinY, fMinZ;
    REAL fDInlet, fDOutlet;
    REAL fRefProportion;
    int fNSubdomains;
    std::string fClass;//NBR 6123 - Building Class
    int fNRefinements;
    VecDouble fInternalField;

public:
    //Default constructor and destructor
    OpenFOAMWriter() = default;
    ~OpenFOAMWriter() = default;

    //Constructor
    OpenFOAMWriter(std::string inputFile);

    void ClearAllFiles();

    void StartFromPreviousResults(REAL dt, REAL endTime, REAL writeInterval);

    //Read the msh file
    bool ParseMSH(const std::string &inputFile, std::vector<Vertex> &vertices, std::vector<Face> &faces);

    //Write the obj file
    bool WriteOBJ(const std::string &outputFile, const std::vector<Vertex> &vertices, const std::vector<Face> &faces);
    
    //Find Max and Min coordinates of the building
    void FindMaxMin(const std::vector<Vertex> &vertices, REAL &maxX, REAL &maxY, REAL &maxZ, REAL &minX, REAL &minY, REAL &minZ);

    //folder 0
    bool WriteInitialConditions(VecDouble &internalField, REAL puniform, REAL nutuniform, REAL kuniform, REAL omegauniform);
    bool UInitial();
    bool pInitial(REAL puniform);
    bool nutInitial(REAL nutuniform);
    bool kInitial(REAL kuniform);
    bool omegaInitial(REAL omegauniform);
    void VelocityNBR6123();

    //folder constant
    bool WriteConstant(REAL density, REAL viscosity);
    bool WriteMomentumTransport();
    bool WritePhysicalProperties(REAL density, REAL viscosity); 

    //folder system
    bool WriteSystem(REAL dInlet, REAL dOutlet, REAL cellSizeX, REAL cellSizeY, REAL cellSizeZ, REAL refproportion, int nref, REAL dt, REAL endTime, REAL writeInterval);
    bool WriteBlockMeshDict(REAL dInlet, REAL dOutlet, REAL cellSizeX, REAL cellSizeY, REAL cellSizeZ);
    bool WriteSurfaceFeaturesDict();
    bool WriteControlDict(REAL dt, REAL endTime, REAL writeInterval);
    bool WriteFvSchemes();
    bool WriteFvSolution();
    bool WriteSnappyHexMeshDict();
    bool WriteMeshQualityDict();
    bool WriteDecomposeParDict();
    bool WriteInitialFields();
    
    int GetNSubdomains() {return fNSubdomains;};

};

#endif