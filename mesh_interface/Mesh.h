#ifndef MESH_H
#define MESH_H

#include <fstream>
#include <stdio.h>
#include "Geometry.h"
#include "ParametricElement.h"


#ifdef _WIN32
#include <direct.h>
#define getCurrentDir _getcwd
#define rm "del "
#else
#include <unistd.h>
#define getCurrentDir getcwd
#define rm "rm "
#endif

enum MeshAlgorithm
{
    AUTO,
    FRONT,
    DELAUNAY,
    ADAPT,
    PACK,
    QUAD,
    BAMG,
    HXT
};

// std::string getCurrentWorkingDir();

// std::pair<std::string, bool> createMesh(Geometry* geometry, const PartitionOfUnity& elementType, const MeshAlgorithm& algorithm = AUTO, std::string geofile = std::string(), const std::string& gmshPath = std::string(), const bool& plotMesh = true, const bool& showInfo = false);


// std::string getCurrentWorkingDir()
// {
// 	char buff[FILENAME_MAX];
// 	getCurrentDir(buff, FILENAME_MAX);
// 	std::string current_working_dir(buff);
// 	return current_working_dir;
// }

// std::pair<std::string, bool> createMesh(Geometry* geometry, const PartitionOfUnity& elementType, const MeshAlgorithm& algorithm, std::string geofile, const std::string& gmshPath, const bool& plotMesh, const bool& showInfo)
// {
// 	const std::unordered_map<PartitionOfUnity, std::string> dimension =
// 	{
// 		//lines
// 		{L2, "3"}, {L3, "3"}, {L4, "3"},
// 		//triangles
// 		{T3, "2"}, {T6, "2"}, {T10, "2"},
// 		//quadrilaterals
// 		{Q4, "2"}, {Q9, "2"}, {Q16, "2"},
// 		//tetrahedrons
// 		{TET4, "3"}, {TET10, "3"}, {TET20, "3"}
// 	};

// 	std::string mshfile;
// 	std::stringstream text;
// 	switch (elementType)
// 	{
// 		case Q4:
// 		case Q9:
// 		case Q16:
// 			text << "Recombine Surface {";
// 			int cont = 0;
// 			for (const auto& pair : geometry->getSurfaces())
// 			{
// 				text << pair.first;
// 				if (cont != (geometry->getNumberOfSurfaces() - 1))
// 					text << ", ";
// 				cont += 1;
// 			}
// 			text << "};";
// 			geometry->appendGmshCode(text.str());
// 			break;
// 	}
	
// 	bool deleteFiles = geofile.empty();
// 	if (deleteFiles)
// 	{
// 		mshfile = geofile + "temp.msh";
// 		geofile = geofile + "temp.geo";
// 	}
// 	else {
// 		mshfile = geofile + ".msh";
// 		geofile = geofile + ".geo";
// 	}
// 	std::ofstream file(geofile);
// 	file << geometry->getGmshCode();
// 	file.close();

// 	std::string gmshExe = (gmshPath.empty()) ? "gmsh" : gmshPath;
// 	std::string cmd = gmshExe;
// 	cmd += " -clscale 1.0 " + geofile + " -o " + mshfile + " -format msh2";

// 	switch (algorithm)
// 	{
// 		case FRONT:
// 			//cmd += " -algo front" + dimension.at(elementType) + "d";
// 			cmd += " -algo front2d";
// 			break;
// 		case DELAUNAY:
// 			cmd += " -algo del" + dimension.at(elementType) + "d";
// 			break;
// 		case ADAPT:
// 			cmd += " -algo meshadapt";
// 			break;
// 		case PACK:
// 			cmd += " -algo pack";
// 			break;
// 		case QUAD:
// 			cmd += " -algo delquad";
// 			break;
// 		case BAMG:
// 			cmd += " -algo bamg";
// 			break;
// 		case HXT:
// 			cmd += " -algo hxt";
// 		break;
// 	}

// 	switch (elementType)
// 	{
// 		case L2:
// 			cmd += " -1 -order 1";
// 			break;
// 		case L3:
// 			cmd += " -1 -order 2";
// 			break;
// 		case L4:
// 			cmd += " -1 -order 3";
// 			break;
// 		case T3:
// 		case Q4:
// 			cmd += " -2 -order 1";
// 			break;
// 		case T6:
// 		case Q9:
// 			cmd += " -2 -order 2";
// 			break;
// 		case T10:
// 		case Q16:
// 			cmd += " -2 -order 3";
// 			break;
// 		case TET4:
// 			cmd += " -3 -order 1";
// 			break;
// 		case TET10:
// 			cmd += " -3 -order 2";
// 			break;
// 		case TET20:
// 			cmd += " -3 -order 3";
// 			break;
// 	}

// 	if (!showInfo) cmd += " -v 0";

// 	system(cmd.c_str());

// 	if (deleteFiles)
// 		system((rm + geofile).c_str());

// 	if (plotMesh)
// 		system((gmshExe + " " + mshfile).c_str());

// 	return std::make_pair(mshfile, deleteFiles);
// }

#endif