#ifndef GMSHTOOLS_H
#define GMSHTOOLS_H

#include "Geometry.h"
#include "CompMesh.h"
    
class CompMesh;
class Element;

#define GMSH_SHIFT 1

namespace GmshTools{

    void MeshReading(Geometry* &geometry_, const std::string& mshfile, CompMesh* cmesh);
    std::vector<std::string> SplitLine(std::string str, std::string delim);
    void ReadNodes(std::ifstream &file, CompMesh * cmesh);
    void ReadElements(Geometry* &geometry_, std::ifstream &file, std::unordered_map<int, std::string> &physicalEntities, CompMesh * cmesh);
    void RenumberConnectivity(CompMesh * cmesh);
    void BoundaryConstrains(CompMesh * cmesh);
    void BoundarySides(CompMesh * cmesh);
    /** @brief Reads geometric mesh file from GMsh (.msh)
     * @param gmesh [output] Reference to a geometric mesh to be filled with elements from msh file
     * @param file_name Relative path to the .msh file you want to read
    */
    void Read(CompMesh& gmesh, const std::string& file_name);
    
    
// private:
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

    Element* CreateElement(CompMesh *cmesh, int64_t index, VecInt &connect);
};


#endif