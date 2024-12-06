#include "OpenFOAMWriter.h"
#include <filesystem> 

OpenFOAMWriter::OpenFOAMWriter(std::string inputFile){
    fInputMsh = inputFile;

    std::filesystem::create_directory("OpenFOAMRun");
    std::filesystem::create_directory("OpenFOAMRun/constant");
    std::filesystem::create_directory("OpenFOAMRun/constant/geometry");
    std::filesystem::create_directory("OpenFOAMRun/0");
    std::filesystem::create_directory("OpenFOAMRun/system");

    fOutputObj = "OpenFOAMRun/constant/geometry/buindings.obj";
    //Delete the output file if it exists
    {
        std::string command = "rm -f " + fOutputObj + ".gz"; 
        system(command.c_str());
    }
}

bool OpenFOAMWriter::ParseMSH(const std::string &inputFile, std::vector<Vertex> &vertices, std::vector<Face> &faces) {
    std::ifstream infile(inputFile);
    if (!infile.is_open()) {
        std::cerr << "Error: Unable to open the .msh file." << std::endl;
        return false;
    }

    std::string line;
    bool readingNodes = false, readingElements = false;

    while (std::getline(infile, line)) {
        // Detect sections of the .msh file
        if (line.find("$Nodes") != std::string::npos) {
            readingNodes = true;
            continue;
        }
        if (line.find("$EndNodes") != std::string::npos) {
            readingNodes = false;
            continue;
        }
        if (line.find("$Elements") != std::string::npos) {
            readingElements = true;
            continue;
        }
        if (line.find("$EndElements") != std::string::npos) {
            readingElements = false;
            continue;
        }

        // Read vertices
        if (readingNodes) {
            std::istringstream iss(line);
            int id;
            double x, y, z;
            if (iss >> id >> x >> y >> z) {
                vertices.push_back({x, y, z});
            }
        }

        // Read faces (triangular faces assumed)
        if (readingElements) {
            std::istringstream iss(line);
            int id, elementType, numTags;
            if (iss >> id >> elementType >> numTags) {
                // Skip tags
                for (int i = 0; i < numTags; i++) {
                    int tag;
                    iss >> tag;
                }
                // If it's a triangular element (type 2 in Gmsh)
                if (elementType == 2) {
                    int v1, v2, v3;
                    if (iss >> v1 >> v2 >> v3) {
                        faces.push_back({v1, v2, v3});
                    }
                }
            }
        }
    }

    infile.close();
    return true;
}


// Function to write to an .obj file
bool OpenFOAMWriter::WriteOBJ(const std::string &outputFile, const std::vector<Vertex> &vertices, const std::vector<Face> &faces) {
    std::ofstream outfile(outputFile);
    if (!outfile.is_open()) {
        std::cerr << "Error: Unable to open the .obj file." << std::endl;
        return false;
    }

    // Write vertices
    for (const auto &vertex : vertices) {
        outfile << "v " << vertex.x << " " << vertex.y << " " << vertex.z << "\n";
    }

    // Write faces
    for (const auto &face : faces) {
        outfile << "f " << face.v1 << " " << face.v2 << " " << face.v3 << "\n";
    }

    outfile.close();
    return true;
}


// Find max and min values of x, y, z
void OpenFOAMWriter::FindMaxMin(const std::vector<Vertex> &vertices, double &maxX, double &maxY, double &maxZ, double &minX, double &minY, double &minZ) {
    maxX = maxY = maxZ = -1e9;
    minX = minY = minZ = 1e9;
    for (const auto &vertex : vertices) {
        if (vertex.x > maxX) maxX = vertex.x;
        if (vertex.y > maxY) maxY = vertex.y;
        if (vertex.z > maxZ) maxZ = vertex.z;
        if (vertex.x < minX) minX = vertex.x;
        if (vertex.y < minY) minY = vertex.y;
        if (vertex.z < minZ) minZ = vertex.z;
    }
}

// Function to write blockMeshDict file
bool OpenFOAMWriter::WriteBlockMeshDict(double dInlet, double dOutlet, double cellSize) {
    
    std::vector<Vertex> vertices;
    std::vector<Face> faces;

    if (!ParseMSH(fInputMsh, vertices, faces)) {
        std::cerr << "Failed to parse the .msh file." << std::endl;
        return 1;
    }

    if (!WriteOBJ(fOutputObj, vertices, faces)) {
        std::cerr << "Failed to write the .obj file." << std::endl;
        return 1;
    }

    //Compress obj file to .gz
    {
        std::string command = "gzip " + fOutputObj; 
        system(command.c_str());
    }   
    
    std::string outputFile = "OpenFOAMRun/system/blockMeshDict";
    std::ofstream outfile(outputFile);
    if (!outfile.is_open()) {
        std::cerr << "Error: Unable to open the .obj file." << std::endl;
        return false;
    }

    double maxX, maxY, maxZ;
    double minX, minY, minZ;
    FindMaxMin(vertices, maxX, maxY, maxZ, minX, minY, minZ);

    minX -= dInlet;
    maxX += dInlet;
    minY -= dInlet;
    maxY += dOutlet;
    maxZ += dInlet;

    int xCells = (maxX-minX)/cellSize;
    int yCells = (maxY-minY)/cellSize;
    int zCells = (maxZ-minZ)/cellSize;

    outfile << "/*--------------------------------*- C++ -*----------------------------------*\n";
    outfile << "  =========                 |\n";
    outfile << "  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n";
    outfile << "   \\\\    /   O peration     | Website:  https://openfoam.org\n";
    outfile << "    \\\\  /    A nd           | Version:  8\n";
    outfile << "     \\\\/     M anipulation  |\n";
    outfile << "\\*---------------------------------------------------------------------------*/\n";
    outfile << "FoamFile\n";
    outfile << "{\n";
    outfile << "    format      ascii;\n";
    outfile << "    class       dictionary;\n";
    outfile << "    object      blockMeshDict;\n";
    outfile << "}\n";
    outfile << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n";
    outfile << "\n\n";
    outfile << "backgroundMesh\n";
    outfile << "{\n";
    outfile << "   xMin    " << minX <<";\n";
    outfile << "   xMax    " << maxX <<";\n";
    outfile << "   yMin    " << minY <<";\n";   
    outfile << "   yMax    " << maxY <<";\n";
    outfile << "   zMin    " << minZ <<";\n";
    outfile << "   zMax    " << maxZ <<";\n";
    outfile << "   xCells  " << xCells <<";\n";
    outfile << "   yCells  " << yCells <<";\n";
    outfile << "   zCells  " << zCells <<";\n";
    outfile << "}\n";
    outfile << "\n\n";
    outfile << "convertToMeters 1;\n";
    outfile << "\n\n";
    outfile << "vertices\n";
    outfile << "(\n";
    outfile << "    ($!backgroundMesh/xMin $!backgroundMesh/yMin $!backgroundMesh/zMin)\n";
    outfile << "    ($!backgroundMesh/xMax $!backgroundMesh/yMin $!backgroundMesh/zMin)\n"; 
    outfile << "    ($!backgroundMesh/xMax $!backgroundMesh/yMax $!backgroundMesh/zMin)\n"; 
    outfile << "    ($!backgroundMesh/xMin $!backgroundMesh/yMax $!backgroundMesh/zMin)\n"; 
    outfile << "\n";
    outfile << "    ($!backgroundMesh/xMin $!backgroundMesh/yMin $!backgroundMesh/zMax)\n"; 
    outfile << "    ($!backgroundMesh/xMax $!backgroundMesh/yMin $!backgroundMesh/zMax)\n"; 
    outfile << "    ($!backgroundMesh/xMax $!backgroundMesh/yMax $!backgroundMesh/zMax)\n"; 
    outfile << "    ($!backgroundMesh/xMin $!backgroundMesh/yMax $!backgroundMesh/zMax)\n"; 
    outfile << ");\n";
    outfile << "\n\n";
    outfile << "blocks\n";
    outfile << "(\n";
    outfile << "    hex (0 1 2 3 4 5 6 7) \n";
    outfile << "    (\n";
    outfile << "    $!backgroundMesh/xCells\n";
    outfile << "    $!backgroundMesh/yCells\n";
    outfile << "    $!backgroundMesh/zCells\n";
    outfile << "    )\n";
    outfile << "    simpleGrading (1 1 1)\n";
    outfile << ");\n";
    outfile << "\n\n";
    outfile << "boundary\n";
    outfile << "(\n";
    outfile << "    inlet\n";
    outfile << "    {\n";
    outfile << "        type patch;\n";
    outfile << "        faces\n";
    outfile << "        (\n";
    outfile << "            (0 3 7 4)\n";
    outfile << "        );\n";
    outfile << "    }\n";
    outfile << "\n";
    outfile << "    outlet\n";
    outfile << "    {\n";
    outfile << "        type patch;\n";
    outfile << "        faces\n";
    outfile << "        (\n";
    outfile << "            (1 5 6 2)\n";
    outfile << "        );\n";
    outfile << "    }\n";
    outfile << "\n";
    outfile << "    ground\n";
    outfile << "    {\n";
    outfile << "        type wall;\n";
    outfile << "        faces\n";
    outfile << "        (\n";
    outfile << "            (0 1 2 3)\n";
    outfile << "        );\n";
    outfile << "    }\n";
    outfile << "\n";
    outfile << "    frontAndBack\n";
    outfile << "    {\n";
    outfile << "        type symmetry;\n";
    outfile << "        faces\n";
    outfile << "        (\n";
    outfile << "            (0 4 5 1)\n";
    outfile << "            (3 2 6 7)\n";
    outfile << "            (4 7 6 5)\n";
    outfile << "        );\n";
    outfile << "    }\n";
    outfile << ");\n";
    outfile << "\n\n";
    outfile << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n";

    return true;
}

bool OpenFOAMWriter::WriteInitialConditions(VecDouble &internalField, double puniform, double nutuniform, double kuniform, double omegauniform){

    UInitial(internalField);
    pInitial(puniform);
    nutInitial(nutuniform);
    kInitial(kuniform);
    omegaInitial(omegauniform);
    return true;
}

bool OpenFOAMWriter::UInitial(VecDouble &internalField){

    std::ofstream UFile("OpenFOAMRun/0/U");
    if (!UFile.is_open()) {
        std::cerr << "Error: Unable to open 0/U file." << std::endl;
        return false;
    }

    UFile << "/*--------------------------------*- C++ -*----------------------------------*\n";
    UFile << "  =========                 |\n";
    UFile << "  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n";
    UFile << "   \\\\    /   O peration     | Website:  https://openfoam.org\n";
    UFile << "    \\\\  /    A nd           | Version:  8\n";
    UFile << "     \\\\/     M anipulation  |\n";
    UFile << "\\*---------------------------------------------------------------------------*/\n";
    UFile << "FoamFile\n";
    UFile << "{\n";
    UFile << "    format      ascii;\n";
    UFile << "    class       volVectorField;\n";
    UFile << "    object      U;\n";
    UFile << "}\n";
    UFile << "\\*---------------------------------------------------------------------------*/\n";
    UFile << "\n\n";
    UFile << "dimensions      [0 1 -1 0 0 0 0];\n";
    UFile << "internalField   uniform (" << internalField[0] << " " << internalField[1] << " " << internalField[2] << ");\n";
    UFile << "boundaryField\n";
    UFile << "{\n";
    UFile << "    inlet\n";
    UFile << "    {\n";
    UFile << "        type            fixedValue;\n";
    UFile << "        value           uniform uniform (" << internalField[0] << " " << internalField[1] << " " << internalField[2] << ");\n";
    UFile << "    }\n";
    UFile << "\n";
    UFile << "    outlet\n";
    UFile << "    {\n";
    UFile << "        type            zeroGradient;\n";
    UFile << "    }\n";
    UFile << "\n";
    UFile << "    wall\n";
    UFile << "    {\n";
    UFile << "        type            noSlip;\n";
    UFile << "    }\n";
    UFile << "\n";
    UFile << "    #includeEtc \"caseDicts/setConstraintTypes\"\n";
    UFile << "}\n";
    UFile << "\n\n";
    UFile << "// ************************************************************************* //\n";

    return true;
}

bool OpenFOAMWriter::pInitial(double puniform){

    std::ofstream PFile("OpenFOAMRun/0/p");
    if (!PFile.is_open()) {
        std::cerr << "Error: Unable to open 0/p file." << std::endl;
        return false;
    }

    PFile << "/*--------------------------------*- C++ -*----------------------------------*\n";
    PFile << "  =========                 |\n";
    PFile << "  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n";
    PFile << "   \\\\    /   O peration     | Website:  https://openfoam.org\n";
    PFile << "    \\\\  /    A nd           | Version:  8\n";
    PFile << "     \\\\/     M anipulation  |\n";
    PFile << "\\*---------------------------------------------------------------------------*/\n";
    PFile << "FoamFile\n";
    PFile << "{\n";
    PFile << "    format      ascii;\n";
    PFile << "    class       volScalarField;\n";
    PFile << "    object      p;\n";
    PFile << "}\n";
    PFile << "\\*---------------------------------------------------------------------------*/\n";
    PFile << "\n\n";
    PFile << "dimensions      [0 2 -2 0 0 0 0];\n";
    PFile << "internalField   uniform " << puniform << ";\n";
    PFile << "boundaryField\n";
    PFile << "{\n";
    PFile << "    inlet\n";
    PFile << "    {\n";
    PFile << "        type            zeroGradient;\n";
    PFile << "    }\n";
    PFile << "\n";
    PFile << "    outlet\n";
    PFile << "    {\n";
    PFile << "        type            totalPressure;\n";
    PFile << "        p0              uniform 0;\n";
    PFile << "        value           uniform 0;\n";
    PFile << "    }\n";
    PFile << "\n";
    PFile << "    wall\n";
    PFile << "    {\n";
    PFile << "        type            zeroGradient;\n";
    PFile << "    }\n";
    PFile << "\n";
    PFile << "    #includeEtc \"caseDicts/setConstraintTypes\"\n";
    PFile << "}\n";
    PFile << "\n\n";
    PFile << "// ************************************************************************* //\n";

    return true;
}

bool OpenFOAMWriter::nutInitial(double  nutuniform){

    std::ofstream nutFile("OpenFOAMRun/0/nut");
    if (!nutFile.is_open()) {
        std::cerr << "Error: Unable to open 0/nut file." << std::endl;
        return false;
    }

    nutFile << "/*--------------------------------*- C++ -*----------------------------------*\n";
    nutFile << "  =========                 |\n";
    nutFile << "  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n";
    nutFile << "   \\\\    /   O peration     | Website:  https://openfoam.org\n";
    nutFile << "    \\\\  /    A nd           | Version:  8\n";
    nutFile << "     \\\\/     M anipulation  |\n";
    nutFile << "\\*---------------------------------------------------------------------------*/\n";
    nutFile << "FoamFile\n";
    nutFile << "{\n";
    nutFile << "    format      ascii;\n";
    nutFile << "    class       volScalarField;\n";
    nutFile << "    object      nut;\n";
    nutFile << "}\n";
    nutFile << "\\*---------------------------------------------------------------------------*/\n";
    nutFile << "\n\n";
    nutFile << "dimensions      [0 2 -1 0 0 0 0];\n";
    nutFile << "internalField   uniform " << nutuniform << ";\n";
    nutFile << "boundaryField\n";
    nutFile << "{\n";
    nutFile << "    inlet\n";
    nutFile << "    {\n";
    nutFile << "        type            zeroGradient;\n";
    nutFile << "    }\n";
    nutFile << "\n";
    nutFile << "    outlet\n";
    nutFile << "    {\n";
    nutFile << "        type            zeroGradient;\n";
    nutFile << "    }\n";
    nutFile << "\n";
    nutFile << "    wall\n";
    nutFile << "    {\n";
    nutFile << "        type            nutUSpaldingWallFunction;\n";
    nutFile << "        value           $internalField;\n";
    nutFile << "    }\n";
    nutFile << "\n";
    nutFile << "    #includeEtc \"caseDicts/setConstraintTypes\"\n";
    nutFile << "}\n";
    nutFile << "\n\n";
    nutFile << "// ************************************************************************* //\n";

    return true;
}

bool OpenFOAMWriter::kInitial(double kuniform){
    std::ofstream kFile("OpenFOAMRun/0/k");
    if (!kFile.is_open()) {
        std::cerr << "Error: Unable to open 0/k file." << std::endl;
        return false;
    }

    kFile << "/*--------------------------------*- C++ -*----------------------------------*\n";
    kFile << "  =========                 |\n";
    kFile << "  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n";
    kFile << "   \\\\    /   O peration     | Website:  https://openfoam.org\n";
    kFile << "    \\\\  /    A nd           | Version:  8\n";
    kFile << "     \\\\/     M anipulation  |\n";
    kFile << "\\*---------------------------------------------------------------------------*/\n";
    kFile << "FoamFile\n";
    kFile << "{\n";
    kFile << "    format      ascii;\n";
    kFile << "    class       volScalarField;\n";
    kFile << "    object      k;\n";
    kFile << "}\n";
    kFile << "\\*---------------------------------------------------------------------------*/\n";
    kFile << "\n\n";
    kFile << "dimensions      [0 2 -2 0 0 0 0];\n";
    kFile << "internalField   uniform " << kuniform << ";\n";
    kFile << "boundaryField\n";
    kFile << "{\n";
    kFile << "    inlet\n";
    kFile << "    {\n";
    kFile << "        type            fixedValue;\n";
    kFile << "        value           $internalField;\n";
    kFile << "    }\n";
    kFile << "\n";
    kFile << "    outlet\n";
    kFile << "    {\n";
    kFile << "        type            zeroGradient;\n";
    kFile << "    }\n";
    kFile << "\n";
    kFile << "    wall\n";
    kFile << "    {\n";
    kFile << "        type            kqRWallFunction;\n";
    kFile << "        value           $internalField;\n";
    kFile << "    }\n";
    kFile << "\n";
    kFile << "    #includeEtc \"caseDicts/setConstraintTypes\"\n";
    kFile << "}\n";
    kFile << "\n\n";
    kFile << "// ************************************************************************* //\n";


    return true;
}

bool OpenFOAMWriter::omegaInitial(double omegauniform){
    std::ofstream omegaFile("OpenFOAMRun/0/omega");
    if (!omegaFile.is_open()) {
        std::cerr << "Error: Unable to open 0/omega file." << std::endl;
        return false;
    }

    omegaFile << "/*--------------------------------*- C++ -*----------------------------------*\n";
    omegaFile << "  =========                 |\n";
    omegaFile << "  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n";
    omegaFile << "   \\\\    /   O peration     | Website:  https://openfoam.org\n";
    omegaFile << "    \\\\  /    A nd           | Version:  8\n";
    omegaFile << "     \\\\/     M anipulation  |\n";
    omegaFile << "\\*---------------------------------------------------------------------------*/\n";
    omegaFile << "FoamFile\n";
    omegaFile << "{\n";
    omegaFile << "    format      ascii;\n";
    omegaFile << "    class       volScalarField;\n";
    omegaFile << "    object      omega;\n";
    omegaFile << "}\n";
    omegaFile << "\\*---------------------------------------------------------------------------*/\n";
    omegaFile << "\n\n";
    omegaFile << "dimensions      [0 0 -1 0 0 0 0];\n";
    omegaFile << "internalField   uniform " << omegauniform << ";\n";
    omegaFile << "boundaryField\n";
    omegaFile << "{\n";
    omegaFile << "    inlet\n";
    omegaFile << "    {\n";
    omegaFile << "        type            fixedValue;\n";
    omegaFile << "        value           $internalField;\n";
    omegaFile << "    }\n";
    omegaFile << "\n";
    omegaFile << "    outlet\n";
    omegaFile << "    {\n";
    omegaFile << "        type            zeroGradient;\n";
    omegaFile << "    }\n";
    omegaFile << "\n";
    omegaFile << "    wall\n";
    omegaFile << "    {\n";
    omegaFile << "        type            omegaWallFunction;\n";
    omegaFile << "        value           $internalField;\n";
    omegaFile << "    }\n";
    omegaFile << "\n";
    omegaFile << "    #includeEtc \"caseDicts/setConstraintTypes\"\n";
    omegaFile << "}\n";
    omegaFile << "\n\n";
    omegaFile << "// ************************************************************************* //\n";


    return true;
}

bool OpenFOAMWriter::WriteConstant(){


    return true;
}