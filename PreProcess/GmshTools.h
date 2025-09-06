#ifndef GMSHTOOLS_H
#define GMSHTOOLS_H

#include "CompMesh.h"
    
class CompMesh;
class Element;

#define GMSH_SHIFT 1

namespace GmshTools{

    /// @brief Renumbers the mesh nodes using METIS to reduce global stiffness bandwidth
    /// @param cmesh Computational mesh
    void RenumberConnectivity(CompMesh * cmesh);

    void BuildNeighbourInformation(CompMesh * cmesh);

    /** @brief Reads geometric mesh file from GMsh (.msh)
     * @param gmesh [output] Reference to a geometric mesh to be filled with elements from msh file
     * @param file_name Relative path to the .msh file you want to read
    */
    void Read(CompMesh& gmesh, const std::string& file_name);
    
    /** @brief Reads geometric mesh file from GMsh (.msh version 3)
     * @param gmesh [output] Reference to a geometric mesh to be filled with elements from msh file
     * @param file_name Relative path to the .msh file you want to read
    */
    void Read3(CompMesh &gmesh, const std::string &file_name);

    /** @brief Reads geometric mesh file from GMsh (.msh version 4.1)
     * @param gmesh [output] Reference to a geometric mesh to be filled with elements from msh file
     * @param file_name Relative path to the .msh file you want to read
    */
    void Read4(CompMesh &gmesh, const std::string &file_name);

};


#endif