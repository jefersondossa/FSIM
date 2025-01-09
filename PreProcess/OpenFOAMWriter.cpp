#include "OpenFOAMWriter.h"
#include <filesystem> 
#include <thread>

OpenFOAMWriter::OpenFOAMWriter(std::string inputFile){
    fInputMsh = inputFile;
}

void OpenFOAMWriter::ClearAllFiles(){
    system("cd OpenFOAMRun && foamCleanTutorials");

    if (std::filesystem::exists("OpenFOAMRun")) {
        // Delete the folder and its contents
        std::filesystem::remove_all("OpenFOAMRun");
        // std::cout << "Folder deleted successfully: " << folderPath << std::endl;
    }
    std::filesystem::create_directory("OpenFOAMRun");
    std::filesystem::create_directory("OpenFOAMRun/constant");
    std::filesystem::create_directory("OpenFOAMRun/constant/geometry");
    // std::filesystem::create_directory("OpenFOAMRun/constant/geometry/buildings.obj");
    std::filesystem::create_directory("OpenFOAMRun/0");
    std::filesystem::create_directory("OpenFOAMRun/system");

    fOutputObj = "OpenFOAMRun/constant/geometry/buildings.obj";
    //Delete the output file if it exists
    {
        std::string command = "rm -f OpenFOAMRun/constant/geometry/buildings.obj.gz"; 
        system(command.c_str());
    }

}

void OpenFOAMWriter::VelocityNBR6123(){

    double dx = fMaxX-fMinX;
    double dy = fMaxY-fMinY;
    double dz = fMaxZ-fMinZ;

    double dmax = std::max({dx, dy, dz});
    double bm;
    double Fr;
    double p;
    double zg = fMaxZ + fDInlet;
    
    //Table 1 and 2 - NBR 6123
    if (dmax <= 20.){
        fClass = "A";
        Fr = 1.;
        if (zg<=250.){
            bm = 1.1;
            p = 0.06;
        } else if (zg<=300.){
            bm = 1.0;
            p = 0.085;
        } else if (zg<=350.){
            bm = 0.94;
            p = 0.10;
        } else if (zg<=420.){
            bm = 0.86;
            p = 0.12;
        } else if (zg<=500.){
            bm = 0.74;
            p = 0.15;
        } else{
            std::cout << "Not provided by NBR 6123\n";
            PanicButton();
        }
    } else if (dmax <= 50.){
        fClass = "B";
        Fr = 0.98;
        if (zg<=250.){
            bm = 1.11;
            p = 0.065;
        } else if (zg<=300.){
            bm = 1.0;
            p = 0.09;
        } else if (zg<=350.){
            bm = 0.94;
            p = 0.105;
        } else if (zg<=420.){
            bm = 0.85;
            p = 0.125;
        } else if (zg<=500.){
            bm = 0.73;
            p = 0.16;
        } else{
            std::cout << "Not provided by NBR 6123\n";
            PanicButton();
        }
    } else {
        fClass = "C";
        Fr = 0.95;
        if (zg<=250.){
            bm = 1.12;
            p = 0.07;
        } else if (zg<=300.){
            bm = 1.0;
            p = 0.1;
        } else if (zg<=350.){
            bm = 0.93;
            p = 0.115;
        } else if (zg<=420.){
            bm = 0.84;
            p = 0.135;
        } else if (zg<=500.){
            bm = 0.71;
            p = 0.175;
        } else{
            std::cout << "Not provided by NBR 6123\n";
            PanicButton();
        }
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
bool OpenFOAMWriter::WriteBlockMeshDict(double dInlet, double dOutlet, double cellSizeX, double cellSizeY, double cellSizeZ) {
    
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
        std::string command = "gzip -r OpenFOAMRun/constant/geometry/buildings.obj"; 
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
    fMaxX = maxX;
    fMaxY = maxY;
    fMaxZ = maxZ;
    fMinX = minX;
    fMinY = minY;
    fMinZ = minZ;
    fDInlet = dInlet;
    fDOutlet = dOutlet;
    
    minX -= dInlet;
    maxX += dOutlet;
    minY -= dInlet;
    maxY += dInlet;
    maxZ += dInlet;
    

    int xCells = (maxX-minX)/cellSizeX;
    int yCells = (maxY-minY)/cellSizeY;
    int zCells = (maxZ-minZ)/cellSizeZ;

    outfile << "/*--------------------------------*- C++ -*----------------------------------*\n";
    outfile << "  =========                 |\n";
    outfile << "  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n";
    outfile << "   \\\\    /   O peration     | Website:  https://openfoam.org\n";
    outfile << "    \\\\  /    A nd           | Version:  11\n";
    outfile << "     \\\\/     M anipulation  |\n";
    outfile << "//---------------------------------------------------------------------------*/\n";
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
    UFile << "    \\\\  /    A nd           | Version:  11\n";
    UFile << "     \\\\/     M anipulation  |\n";
    UFile << "//---------------------------------------------------------------------------*/\n";
    UFile << "FoamFile\n";
    UFile << "{\n";
    UFile << "    format      ascii;\n";
    UFile << "    class       volVectorField;\n";
    UFile << "    object      U;\n";
    UFile << "}\n";
    UFile << "//---------------------------------------------------------------------------*/\n";
    UFile << "\n\n";
   
    UFile << "dimensions      [0 1 -1 0 0 0 0];\n";
    UFile << "internalField   uniform (" << internalField[0] << " " << internalField[1] << " " << internalField[2] << ");\n";
    UFile << "boundaryField\n";
    UFile << "{\n";
    UFile << "    inlet\n";
    UFile << "    {\n";
    UFile << "        type            codedFixedValue;\n";
    UFile << "        value           uniform ( 0 0 0 );\n";
    UFile << "        redirectType    codedFixedValue;\n";
    UFile << "        name            parabolicVelocity;\n";
    UFile << "        code\n";
    UFile << "        #{\n";
    UFile << "              const vectorField& Cf = patch().Cf();\n";
    UFile << "              vectorField& field = *this;\n\n";
    UFile << "              const scalar V0 = " << internalField[0] << ";\n";
    UFile << "              const scalar S1 = 1.;\n";
    UFile << "              const scalar S3 = 1.;\n";
    UFile << "              const scalar bm = 1.;\n";
    UFile << "              const scalar Fr = 0.16;\n";
    UFile << "              const scalar p = 0.65;\n\n";
    UFile << "              forAll(Cf, faceI)\n";
    UFile << "              {\n";
    UFile << "                  const scalar x = Cf[faceI][0];\n";
    UFile << "                  const scalar y = Cf[faceI][1];\n";
    UFile << "                  const scalar z = Cf[faceI][2];\n";
    UFile << "                  field[faceI] = vector(V0*S1*S3*bm*Fr*pow(z/10,p),0,0);\n";
    UFile << "              }\n";
    UFile << "        #};\n";
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
    PFile << "    \\\\  /    A nd           | Version:  11\n";
    PFile << "     \\\\/     M anipulation  |\n";
    PFile << "//---------------------------------------------------------------------------*/\n";
    PFile << "FoamFile\n";
    PFile << "{\n";
    PFile << "    format      ascii;\n";
    PFile << "    class       volScalarField;\n";
    PFile << "    object      p;\n";
    PFile << "}\n";
    PFile << "//---------------------------------------------------------------------------*/\n";
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
    nutFile << "    \\\\  /    A nd           | Version:  11\n";
    nutFile << "     \\\\/     M anipulation  |\n";
    nutFile << "//---------------------------------------------------------------------------*/\n";
    nutFile << "FoamFile\n";
    nutFile << "{\n";
    nutFile << "    format      ascii;\n";
    nutFile << "    class       volScalarField;\n";
    nutFile << "    object      nut;\n";
    nutFile << "}\n";
    nutFile << "//---------------------------------------------------------------------------*/\n";
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
    kFile << "    \\\\  /    A nd           | Version:  11\n";
    kFile << "     \\\\/     M anipulation  |\n";
    kFile << "//---------------------------------------------------------------------------*/\n";
    kFile << "FoamFile\n";
    kFile << "{\n";
    kFile << "    format      ascii;\n";
    kFile << "    class       volScalarField;\n";
    kFile << "    object      k;\n";
    kFile << "}\n";
    kFile << "//---------------------------------------------------------------------------*/\n";
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
    omegaFile << "    \\\\  /    A nd           | Version:  11\n";
    omegaFile << "     \\\\/     M anipulation  |\n";
    omegaFile << "//---------------------------------------------------------------------------*/\n";
    omegaFile << "FoamFile\n";
    omegaFile << "{\n";
    omegaFile << "    format      ascii;\n";
    omegaFile << "    class       volScalarField;\n";
    omegaFile << "    object      omega;\n";
    omegaFile << "}\n";
    omegaFile << "//---------------------------------------------------------------------------*/\n";
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

bool OpenFOAMWriter::WriteConstant(double density, double viscosity){

    WriteMomentumTransport();
    WritePhysicalProperties(density, viscosity);

    return true;
}

bool OpenFOAMWriter::WriteMomentumTransport(){
    std::ofstream outFile("OpenFOAMRun/constant/momentumTransport");
    if (!outFile.is_open()) {
        std::cerr << "Error: Unable to open constant/momentumTransport file." << std::endl;
        return false;
    }

    outFile << "/*--------------------------------*- C++ -*----------------------------------*\n";
    outFile << "  =========                 |\n";
    outFile << "  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n";
    outFile << "   \\\\    /   O peration     | Website:  https://openfoam.org\n";
    outFile << "    \\\\  /    A nd           | Version:  11\n";
    outFile << "     \\\\/     M anipulation  |\n";
    outFile << "//---------------------------------------------------------------------------*/\n";
    outFile << "FoamFile\n";
    outFile << "{\n";
    outFile << "    format      ascii;\n";
    outFile << "    class       dictionary;\n";
    outFile << "    object      momentumTransport;\n";
    outFile << "}\n";
    outFile << "//---------------------------------------------------------------------------//\n";
    outFile << "\n\n";
    outFile << "simulationType  RAS;\n";
    outFile << "RAS\n";
    outFile << "{\n";
    outFile << "    RASModel        kOmegaSST;\n";
    outFile << "    turbulence      on;\n";
    outFile << "    printCoeffs     on;\n";
    outFile << "}\n";
    outFile << "\n\n";
    outFile << "// ************************************************************************* //\n";


    return true;
}

bool OpenFOAMWriter::WritePhysicalProperties(double density, double viscosity){

    std::ofstream outFile("OpenFOAMRun/constant/physicalProperties");
    if (!outFile.is_open()) {
        std::cerr << "Error: Unable to open constant/physicalProperties file." << std::endl;
        return false;
    }

    outFile << "/*--------------------------------*- C++ -*----------------------------------*\n";
    outFile << "  =========                 |\n";
    outFile << "  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n";
    outFile << "   \\\\    /   O peration     | Website:  https://openfoam.org\n";
    outFile << "    \\\\  /    A nd           | Version:  11\n";
    outFile << "     \\\\/     M anipulation  |\n";
    outFile << "//---------------------------------------------------------------------------*/\n";
    outFile << "FoamFile\n";
    outFile << "{\n";
    outFile << "    format      ascii;\n";
    outFile << "    class       dictionary;\n";
    outFile << "    object      physicalProperties;\n";
    outFile << "}\n";
    outFile << "//---------------------------------------------------------------------------//\n";
    outFile << "\n\n";
    outFile << "transportModel  Newtonian;\n";
    outFile << "nu              nu [ 0 2 -1 0 0 0 0 ] " << viscosity << ";\n";
    outFile << "rho             rho [ 1 -3 0 0 0 0 0 ] " << density << ";\n";
    outFile << "\n\n";
    outFile << "// ************************************************************************* //\n";

    return true;
}

bool OpenFOAMWriter::WriteSurfaceFeaturesDict(){
    std::ofstream outFile("OpenFOAMRun/system/surfaceFeaturesDict");
    if (!outFile.is_open()) {
        std::cerr << "Error: Unable to open system/surfaceFeaturesDict file." << std::endl;
        return false;
    }

    outFile << "/*--------------------------------*- C++ -*----------------------------------*\n";
    outFile << "  =========                 |\n";
    outFile << "  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n";
    outFile << "   \\\\    /   O peration     | Website:  https://openfoam.org\n";
    outFile << "    \\\\  /    A nd           | Version:  11\n";
    outFile << "     \\\\/     M anipulation  |\n";
    outFile << "//---------------------------------------------------------------------------*/\n";
    outFile << "FoamFile\n";
    outFile << "{\n";
    outFile << "    format      ascii;\n";
    outFile << "    class       dictionary;\n";
    outFile << "    object      surfaceFeaturesDict;\n";
    outFile << "}\n";
    outFile << "//---------------------------------------------------------------------------//\n";
    outFile << "\n\n";
    // outFile << "surfaceFeatures\n";
    // outFile << "{\n";
    outFile << "    surfaces (\"buildings.obj\");\n";
    outFile << "    #includeEtc \"caseDicts/surface/surfaceFeaturesDict.cfg\"\n";
    // outFile << "}\n";
    outFile << "\n\n";
    outFile << "// ************************************************************************* //\n";


    return true;
}

bool OpenFOAMWriter::WriteControlDict(double dt, double endTime, double writeInterval){
    std::ofstream outFile("OpenFOAMRun/system/controlDict");
    if (!outFile.is_open()) {
        std::cerr << "Error: Unable to open system/controlDict file." << std::endl;
        return false;
    }

    outFile << "/*--------------------------------*- C++ -*----------------------------------*\n";
    outFile << "  =========                 |\n";
    outFile << "  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n";
    outFile << "   \\\\    /   O peration     | Website:  https://openfoam.org\n";
    outFile << "    \\\\  /    A nd           | Version:  11\n";
    outFile << "     \\\\/     M anipulation  |\n";
    outFile << "//---------------------------------------------------------------------------*/\n";
    outFile << "FoamFile\n";
    outFile << "{\n";
    outFile << "    format      ascii;\n";
    outFile << "    class       dictionary;\n";
    outFile << "    object      controlDict;\n";
    outFile << "}\n";
    outFile << "//---------------------------------------------------------------------------//\n";
    outFile << "\n\n";
    outFile << "application     foamRun;\n";
    outFile << "solver          incompressibleFluid;\n";

    outFile << "startFrom       latestTime;\n";
    outFile << "startTime       0;\n";
    outFile << "stopAt          endTime;\n";
    outFile << "endTime         " << endTime <<";\n";
    outFile << "deltaT          " << dt << ";\n";
    outFile << "writeControl    adjustableRunTime;\n";
    outFile << "writeInterval   "<< writeInterval << ";\n";
    outFile << "purgeWrite      0;\n";
    outFile << "writeFormat     ascii;\n";
    outFile << "writePrecision  8;\n";
    outFile << "writeCompression off;\n";
    outFile << "timeFormat      general;\n";
    outFile << "timePrecision   6;\n";
    outFile << "runTimeModifiable true;\n";
    outFile << "adjustTimeStep  yes;\n";
    outFile << "maxCo           20;\n";
    
    outFile << "\n\n";
    outFile << "// ************************************************************************* //\n";

    outFile << "///////////////////////////////////////////////////////////////////////////\n";
    outFile << "functions\n";
    outFile << "{\n";
    outFile << "minmaxdomain_scalar\n";
    outFile << "{\n";
    outFile << "    type            volFieldValue;\n";
    outFile << "    libs            (\"libfieldFunctionObjects.so\");\n";
    outFile << " \n";
    outFile << "    enabled 	    true; 	//true or false\n";
    outFile << "    log             true;	//write to screen\n";
    outFile << "    writeControl    timeStep; \n";
    outFile << "    writeInterval   1; \n\n";
    outFile << "    writeFields     false;	//write solution to field value - Not needed when only reporting value to screen\n";
    outFile << "    writeLocation   true;	//write location in the output file \n";
    outFile << "    select      all; \n";
    outFile << "    operation       none; \n\n";
    outFile << "fields \n";
    outFile << "( \n";
    outFile << " p k omega nut\n";
    outFile << " );\n";
    outFile << "} \n\n";
    outFile << "minmaxdomain_vector \n";
    outFile << "{ \n";
    outFile << "    type            volFieldValue; \n";
    outFile << "    libs            (\"libfieldFunctionObjects.so\"); \n";
    outFile << "    enabled 	    true; 	//true or false \n";
    outFile << "    log             true;	//write to screen \n";
    outFile << "    writeControl    timeStep; \n";
    outFile << "    writeInterval   1; \n";
    outFile << "    writeFields     false;	//write solution to field value - Not needed when only reporting value to screen \n";
    outFile << "    writeLocation   true;	//write location in the output file \n";
    outFile << "    select      all; \n";
    outFile << "    operation       none; \n";
    outFile << "    fields \n";
    outFile << "    (U); \n";
    outFile << "} \n";
    outFile << "mindomain_scalar \n";
    outFile << "{ \n";
    outFile << "    $minmaxdomain_scalar \n";
    outFile << "    operation       min; \n";
    outFile << "} \n";
    outFile << "mindomain_vector \n";
    outFile << "{ \n";
    outFile << "    $minmaxdomain_vector \n";
    outFile << "    operation       minMag; \n";
    outFile << "} \n";
    outFile << "maxdomain_scalar \n";
    outFile << "{ \n";
    outFile << "    $minmaxdomain_scalar \n";
    outFile << "    operation       max; \n";
    outFile << "} \n";
    outFile << "maxdomain_vector \n";
    outFile << "{ \n";
    outFile << "    $minmaxdomain_vector \n";
    outFile << "    operation       maxMag; \n";
    outFile << "} \n";
    outFile << "/////////////////////////////////////////////////////////////////////////// \n";
    outFile << "forces_object \n";
    outFile << "{ \n";
    outFile << "    type forces; \n";
    outFile << "    libs (\"libforces.so\"); \n\n";
    outFile << "    enabled true;\n";
    outFile << "    writeControl   timeStep; \n";
    outFile << "    writeInterval  1; \n";
    outFile << "    patches (\"square\"); \n";
    outFile << "    pName p; \n";
    outFile << "    Uname U; \n";
    outFile << "    ////Density only for incompressible flows \n";
    outFile << "    rho rhoInf; \n";
    outFile << "    rhoInf 1; \n";
    outFile << "    CofR (0 0 0);//Centre of rotation  \n";
    outFile << "/////////////////////////////////////////////////////////////////////////// \n";
    outFile << "probes1 \n";
    outFile << "{ \n";
    outFile << "    type probes; \n";
    outFile << "    libs (\"libsampling.so\"); \n";
    outFile << "    probeLocations \n";
    outFile << "    (\n        (0.025 0 0)\n        (0.05 0 0)\n        (0.075 0 0)\n        (0.1 0 0)\n        (0.15 0 0)\n        (0.3 0 0)\n        (0.5 0 0)\n); \n";
    outFile << "    fields \n";
    outFile << "    (p U); \n";
    outFile << "    writeControl       timeStep; \n";
    outFile << "    writeInterval      1; \n";
    outFile << "} \n";
    outFile << "/////////////////////////////////////////////////////////////////////////// \n";
    outFile << "//Saves yplus field at the given interval \n";
    outFile << "    yplus_field \n";
    outFile << "    { \n";
    outFile << "         type yPlus;\n";
    outFile << "        libs (\"libutilityFunctionObjects.so\"); \n";
    outFile << "        enabled true; \n";
    outFile << "        writeControl outputTime; \n";
    outFile << "    } \n";
    outFile << "///////////////////////////////////////////////////////////////////////////\n";
    outFile << "//Print y+ stats every iteration without saving the field \n";
    outFile << "    yplus_stats \n";
    outFile << "    { \n";
    outFile << "        type yPlus; \n";
    outFile << "        libs (\"libutilityFunctionObjects.so\"); \n";
    outFile << "        enabled true; \n";
    outFile << "        writeControl timeStep; \n";
    outFile << "        writeInterval 1; \n";
    outFile << "        objects (); \n";
    outFile << "        log true; \n";
    outFile << "    } \n";
    outFile << "/////////////////////////////////////////////////////////////////////////// \n";
    outFile << "Q\n{ \n";
    outFile << "    type                Q; \n";
    outFile << "    libs  (\"libutilityFunctionObjects.so\"); \n";
    outFile << "    enabled true; \n";
    outFile << "    writeControl outputTime; \n";
    outFile << "} \n";
    outFile << "/////////////////////////////////////////////////////////////////////////// \n";
    outFile << "cellVol\n{ \n";
    outFile << "    type                writeCellVolumes; \n";
    outFile << "    libs  (\"libutilityFunctionObjects.so\"); \n";
    outFile << "    enabled true; \n";
    outFile << "    writeControl outputTime; \n";
    outFile << "} \n";
    outFile << "/////////////////////////////////////////////////////////////////////////// \n";
    outFile << "wallShearStress1\n{ \n";
    outFile << "    type                wallShearStress; \n";
    outFile << "    libs  (\"libutilityFunctionObjects.so\"); \n";
    outFile << "    enabled true; \n";
    outFile << "    writeControl outputTime; \n";
    outFile << "} \n";
    outFile << "/////////////////////////////////////////////////////////////////////////// \n";
    outFile << "    fieldAverage \n";
    outFile << "    { \n";
    outFile << "        type            fieldAverage; \n";
    outFile << "        libs (\"libfieldFunctionObjects.so\"); \n";
    outFile << "        enabled         true; \n";
    outFile << "        writeControl   outputTime; \n";
    outFile << "	    timeStart       2; \n";
    outFile << "        fields \n";
    outFile << "        ( \n";
    outFile << "            U \n";
    outFile << "            { \n";
    outFile << "                mean        on; \n";
    outFile << "                prime2Mean  on; \n";
    outFile << "                base        time; \n";
    outFile << "            } \n";
    outFile << "            p \n";
    outFile << "            { \n";
    outFile << "                mean        on; \n";
    outFile << "                prime2Mean  on; \n";
    outFile << "                base        time; \n";
    outFile << "            } \n";
    outFile << "            nut \n";
    outFile << "            { \n";
    outFile << "                mean        on; \n";
    outFile << "                prime2Mean  on; \n";
    outFile << "                base        time; \n";
    outFile << "            } \n";
    outFile << "            k \n";
    outFile << "            { \n";
    outFile << "                mean        on; \n";
    outFile << "                prime2Mean  off; \n";
    outFile << "                base        time; \n";
    outFile << "            } \n";
    outFile << "            omega \n";
    outFile << "            { \n";
    outFile << "                mean        on; \n";
    outFile << "                prime2Mean  off; \n";
    outFile << "                base        time; \n";
    outFile << "            } \n";
    outFile << "            Q \n";
    outFile << "            { \n";
    outFile << "                mean        on; \n";
    outFile << "                prime2Mean  off; \n";
    outFile << "                base        time; \n";
    outFile << "            } \n";
    outFile << "            yPlus \n";
    outFile << "            { \n";
    outFile << "                mean        on; \n";
    outFile << "                prime2Mean  off; \n";
    outFile << "                base        time; \n";
    outFile << "            } \n";
    outFile << "            wallShearStress \n";
    outFile << "            { \n";
    outFile << "                mean        on; \n";
    outFile << "                prime2Mean  off; \n";
    outFile << "                base        time; \n";
    outFile << "            } \n";
    outFile << "        ); \n";
    outFile << "    } \n";
    outFile << "}; \n";

    return true;
}

bool OpenFOAMWriter::WriteFvSchemes(){
    std::ofstream outFile("OpenFOAMRun/system/fvSchemes");
    if (!outFile.is_open()) {
        std::cerr << "Error: Unable to open system/fvSchemes file." << std::endl;
        return false;
    }

    outFile << "/*--------------------------------*- C++ -*----------------------------------*\n";
    outFile << "  =========                 |\n";
    outFile << "  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n";
    outFile << "   \\\\    /   O peration     | Website:  https://openfoam.org\n";
    outFile << "    \\\\  /    A nd           | Version:  11\n";
    outFile << "     \\\\/     M anipulation  |\n";
    outFile << "//---------------------------------------------------------------------------*/\n";
    outFile << "FoamFile\n";
    outFile << "{\n";
    outFile << "    format      ascii;\n";
    outFile << "    class       dictionary;\n";
    outFile << "    object      fvSchemes;\n";
    outFile << "}\n";
    outFile << "//---------------------------------------------------------------------------*/\n";
    outFile << "\n\n";
    outFile << "ddtSchemes\n";
    outFile << "{\n";
    outFile << "    default         Euler;\n";
    outFile << "}\n";
    outFile << "\n";
    outFile << "gradSchemes\n";
    outFile << "{\n";
    outFile << "    default         cellLimited Gauss linear 1;\n";
    outFile << "    grad(p)         cellLimited Gauss linear 0.5;\n";
    outFile << "}\n";
    outFile << "\n";
    outFile << "divSchemes\n";
    outFile << "{\n";
    outFile << "    default         none;\n";
    outFile << "    div(phi,U)      Gauss linearUpwindV default;\n";
    outFile << "    div(phi,k)      Gauss upwind;\n";
    outFile << "    div(phi,omega)  Gauss upwind;\n";
    outFile << "    div((nuEff*dev(T(grad(U))))) Gauss linear;\n";
    outFile << "}\n";
    outFile << "\n";
    outFile << "laplacianSchemes\n";
    outFile << "{\n";
    outFile << "    default         Gauss linear limited 1;\n";
    outFile << "}\n";
    outFile << "\n";
    outFile << "interpolationSchemes\n";
    outFile << "{\n";
    outFile << "    default         linear;\n";
    outFile << "}\n";
    outFile << "\n";
    outFile << "snGradSchemes\n";
    outFile << "{\n";
    outFile << "    default         limited 1;\n";
    outFile << "}\n";
    outFile << "\n";
    outFile << "wallDist\n";
    outFile << "{\n";
    outFile << "    method          meshWave;\n";
    outFile << "}\n";
    outFile << "\n\n";
    outFile << "// ************************************************************************* //\n";


    return true;
}

bool OpenFOAMWriter::WriteFvSolution(){
    std::ofstream outFile("OpenFOAMRun/system/fvSolution");
    if (!outFile.is_open()) {
        std::cerr << "Error: Unable to open system/fvSolution file." << std::endl;
        return false;
    }

    outFile << "/*--------------------------------*- C++ -*----------------------------------*\n";
    outFile << "  =========                 |\n";
    outFile << "  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n";
    outFile << "   \\\\    /   O peration     | Website:  https://openfoam.org\n";
    outFile << "    \\\\  /    A nd           | Version:  11\n";
    outFile << "     \\\\/     M anipulation  |\n";
    outFile << "//---------------------------------------------------------------------------*/\n";
    outFile << "FoamFile\n";
    outFile << "{\n";
    outFile << "    format      ascii;\n";
    outFile << "    class       dictionary;\n";
    outFile << "    object      fvSolution;\n";
    outFile << "}\n";
    outFile << "//---------------------------------------------------------------------------*/\n";
    outFile << "\n\n";
    outFile << "solvers\n";
    outFile << "{\n";
    outFile << "    p\n";
    outFile << "    {\n";
    outFile << "        solver          GAMG;\n";
    outFile << "        tolerance       1e-6;\n";
    outFile << "        relTol          0.01;\n";
    outFile << "        smoother        GaussSeidel;\n";
    outFile << "        nPreSweeps      0;\n";
    outFile << "        nPostSweeps     2;\n";
    outFile << "        cacheAgglomeration on;\n";
    outFile << "        agglomerator    faceAreaPair;\n";
    outFile << "        nCellsInCoarsestLevel 100;\n";
    outFile << "        mergeLevels     1;\n";
    outFile << "        minIter		 3;\n";
    outFile << "    }\n";
    outFile << "\n";
    outFile << "    pFinal\n";
    outFile << "    {\n";
    outFile << "        $p;\n";
    outFile << "        tolerance       1e-6;\n";
    outFile << "        relTol          0;\n";
    outFile << "        minIter	        1;\n";
    outFile << "    }\n";
    outFile << "\n";
    outFile << "    \"(U|UFinal)\"\n";
    outFile << "    {\n";
    outFile << "        solver          PBiCGStab;\n";
    outFile << "        preconditioner  DILU;\n";
    outFile << "        tolerance       1e-8;\n";
    outFile << "        relTol          0;\n";
    outFile << "        minIter	        3;\n";
    outFile << "    }\n";
    outFile << "\n";
    outFile << "    k\n";
    outFile << "    {\n";
    outFile << "        solver          PBiCGStab;\n";
    outFile << "        preconditioner  DILU;\n";
    outFile << "        tolerance       1e-8;\n";
    outFile << "        relTol          0.001;\n";
    outFile << "    }\n";
    outFile << "\n";
    outFile << "    kFinal\n";
    outFile << "    {\n";
    outFile << "        solver          PBiCGStab;\n";
    outFile << "        preconditioner  DILU;\n";
    outFile << "        tolerance       1e-8;\n";
    outFile << "        relTol          0;\n";
    outFile << "        minIter	        2;\n";
    outFile << "    }\n";
    outFile << "\n";
    outFile << "    omega\n";
    outFile << "    {\n";
    outFile << "        solver          PBiCG;\n";
    outFile << "        preconditioner  DILU;\n";
    outFile << "        tolerance       1e-8;\n";
    outFile << "        relTol          0;\n";
    outFile << "    }\n";
    outFile << "\n";
    outFile << "    omegaFinal\n";
    outFile << "    {\n";
    outFile << "        solver          PBiCGStab;\n";
    outFile << "        preconditioner  DILU;\n";
    outFile << "        tolerance       1e-8;\n";
    outFile << "        relTol          0;\n";
    outFile << "        minIter	        2;\n";
    outFile << "    }\n";
    outFile << "\n";
    outFile << "}\n";
    outFile << "\n\n";
    outFile << "PIMPLE\n";
    outFile << "{\n";
    outFile << "    momentumPredictor yes;\n";
    outFile << "    consistent yes;\n";
    outFile << "    nOuterCorrectors 20;\n";
    outFile << "    nCorrectors 2;\n";
    outFile << "    nNonOrthogonalCorrectors 1;\n";
    outFile << "    outerCorrectorResidualControl\n";
    outFile << "    {\n";
    outFile << "        \"(U|p)\"\n";
    outFile << "        {\n";
    outFile << "            tolerance 1e-3;\n";
    outFile << "            relTol 0.01;\n";
    outFile << "        }\n";
    outFile << "    }\n";
    outFile << "}\n";
    outFile << "\n\n";
    outFile << "relaxationFactors";
    outFile << "{\n";
    outFile << "    fields\n";
    outFile << "    {\n";
    outFile << "        \".*\"	0.9;\n";
    outFile << "    }\n";
    outFile << "    equations\n";
    outFile << "    {\n";
    outFile << "        \".*\"	0.9;\n";
    outFile << "    }\n";
    outFile << "}\n";
    outFile << "\n\n";
    outFile << "// ************************************************************************* //\n";





    return true;
}

bool OpenFOAMWriter::WriteSnappyHexMeshDict(){
    std::ofstream outFile("OpenFOAMRun/system/snappyHexMeshDict");
    if (!outFile.is_open()) {
        std::cerr << "Error: Unable to open system/snappyHexMeshDict file." << std::endl;
        return false;
    }

    outFile << "/*--------------------------------*- C++ -*----------------------------------*\n";
    outFile << "  =========                 |\n";
    outFile << "  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n";
    outFile << "   \\\\    /   O peration     | Website:  https://openfoam.org\n";
    outFile << "    \\\\  /    A nd           | Version:  11\n";
    outFile << "     \\\\/     M anipulation  |\n";
    outFile << "//---------------------------------------------------------------------------*/\n";
    outFile << "FoamFile\n";
    outFile << "{\n";
    outFile << "    format      ascii;\n";
    outFile << "    class       dictionary;\n";
    outFile << "    object      snappyHexMeshDict;\n";
    outFile << "}\n";
    outFile << "//---------------------------------------------------------------------------*/\n";
    outFile << "\n\n";
    outFile << "#includeEtc \"caseDicts/mesh/generation/snappyHexMeshDict.cfg\"\n\n\n";
    outFile << "castellatedMesh on;\n";
    outFile << "snap            on;\n";
    outFile << "addLayers       off;\n";
    outFile << "\n\n";
    outFile << "geometry\n";
    outFile << "{\n";
    outFile << "    buildings\n";
    outFile << "    {\n";
    outFile << "        type triSurfaceMesh;\n";
    outFile << "        file \"buildings.obj\";\n";
    outFile << "    }\n";
    outFile << "    refinementBox\n";
    outFile << "    {\n";
    outFile << "        type searchableBox;\n";
    outFile << "        min  (" << fMinX-fRefProportion*fDInlet << " " << fMinY-fRefProportion*fDInlet << " " << fMinZ-fRefProportion*fDInlet << ");\n";
    outFile << "        max  (" << fMaxX+fRefProportion*fDOutlet << " " << fMaxY+fRefProportion*fDInlet << " " << fMaxZ+fRefProportion*fDInlet << ");\n";
    outFile << "    }\n";
    outFile << "};\n";
    outFile << "\n\n";
    outFile << "castellatedMeshControls\n";
    outFile << "{\n";
    outFile << "    features\n";
    outFile << "    (\n";
    outFile << "        {file  \"buildings.eMesh\"; level 1; }\n";
    outFile << "    );\n";
    outFile << "    refinementSurfaces\n";
    outFile << "    {\n";
    outFile << "        buildings\n";
    outFile << "        {\n";
    outFile << "            level (3 3);\n";
    outFile << "            patchInfo { type wall; }\n";
    outFile << "        }\n";
    outFile << "    }\n";
    outFile << "    refinementRegions\n";
    outFile << "    {\n";
    outFile << "        refinementBox\n";
    outFile << "        {\n";
    outFile << "            mode inside;\n";
    outFile << "            level 3;\n";
    outFile << "        }\n";
    outFile << "    }\n";
    outFile << "    insidePoint (1 1 1);\n";
    outFile << "}\n";
    outFile << "\n\n";
    outFile << "snapControls\n";
    outFile << "{\n";
    outFile << "    explicitFeatureSnap    true;\n";
    outFile << "    implicitFeatureSnap    false;\n";
    outFile << "}\n";
    outFile << "\n\n";
    outFile << "addLayersControls\n";
    outFile << "{\n";
    outFile << "    layers\n";
    outFile << "    {\n";
    outFile << "        \"CAD.*\"\n";
    outFile << "        {\n";
    outFile << "            nSurfaceLayers 2;\n";
    outFile << "        }\n";
    outFile << "    }\n";
    outFile << "    relativeSizes       true;\n";
    outFile << "    expansionRatio      1.2;\n";
    outFile << "    finalLayerThickness 0.5;\n";
    outFile << "    minThickness        0.001;\n";
    outFile << "}\n";
    outFile << "\n\n";
    outFile << "meshQualityControls\n";
    outFile << "{}\n";
    outFile << "\n\n";
    outFile << "writeFlags\n";
    outFile << "();\n";
    outFile << "\n\n";
    outFile << "mergeTolerance 1e-6;\n";
    outFile << "\n\n";
    outFile << "// ************************************************************************* //\n";




    return true;
}

bool OpenFOAMWriter::WriteMeshQualityDict(){
    std::ofstream outFile("OpenFOAMRun/system/meshQualityDict");
    if (!outFile.is_open()) {
        std::cerr << "Error: Unable to open system/meshQualityDict file." << std::endl;
        return false;
    }

    outFile << "/*--------------------------------*- C++ -*----------------------------------*\n";
    outFile << "  =========                 |\n";
    outFile << "  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n";
    outFile << "   \\\\    /   O peration     | Website:  https://openfoam.org\n";
    outFile << "    \\\\  /    A nd           | Version:  11\n";
    outFile << "     \\\\/     M anipulation  |\n";
    outFile << "//---------------------------------------------------------------------------*/\n";
    outFile << "FoamFile\n";
    outFile << "{\n";
    outFile << "    format      ascii;\n";
    outFile << "    class       dictionary;\n";
    outFile << "    object      meshQualityDict;\n";
    outFile << "}\n";
    outFile << "//---------------------------------------------------------------------------*/\n";
    outFile << "\n\n";
    outFile << "#includeEtc \"caseDicts/mesh/generation/meshQualityDict.cfg\"\n";
    outFile << "\n\n";
    outFile << "// ************************************************************************* //\n";


    return true;
}

bool OpenFOAMWriter::WriteDecomposeParDict(){
    std::ofstream outFile("OpenFOAMRun/system/decomposeParDict");
    if (!outFile.is_open()) {
        std::cerr << "Error: Unable to open system/decomposeParDict file." << std::endl;
        return false;
    }
    fNSubdomains = std::thread::hardware_concurrency()/2;

    // fNSubdomains = sysconf(_SC_NPROCESSORS_ONLN);


    outFile << "/*--------------------------------*- C++ -*----------------------------------*\n";
    outFile << "  =========                 |\n";
    outFile << "  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n";
    outFile << "   \\\\    /   O peration     | Website:  https://openfoam.org\n";
    outFile << "    \\\\  /    A nd           | Version:  11\n";
    outFile << "     \\\\/     M anipulation  |\n";
    outFile << "//---------------------------------------------------------------------------*/\n";
    outFile << "FoamFile\n";
    outFile << "{\n";
    outFile << "    format      ascii;\n";
    outFile << "    class       dictionary;\n";
    outFile << "    object      decomposeParDict;\n";
    outFile << "}\n";
    outFile << "//---------------------------------------------------------------------------*/\n";
    outFile << "\n\n";
    outFile << "numberOfSubdomains " << fNSubdomains << ";\n";
    outFile << "\n\n";
    outFile << "method          scotch;\n";
    outFile << "\n\n";
    outFile << "// ************************************************************************* //\n";


    return true;
}

bool OpenFOAMWriter::WriteSystem(double dInlet, double dOutlet, double cellSizeX, double cellSizeY, double cellSizeZ, double refproportion, int nref, double dt, double endTime, double writeInterval){
    fRefProportion = refproportion;
    fNRefinements = nref;
    WriteBlockMeshDict(dInlet, dOutlet, cellSizeX, cellSizeY, cellSizeZ);
    WriteControlDict(dt, endTime, writeInterval);
    WriteFvSchemes();
    WriteFvSolution();
    WriteSurfaceFeaturesDict();
    WriteSnappyHexMeshDict();
    WriteMeshQualityDict();
    WriteDecomposeParDict();

    return true;
}


void OpenFOAMWriter::StartFromPreviousResults(double dt, double endTime, double writeInterval){
    WriteControlDict(dt, endTime, writeInterval);

}