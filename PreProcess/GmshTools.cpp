#include "GmshTools.h"
#include "ShapeHexahedron.h"
#include "ShapeOneDLin.h"
#include "ShapeOneDQua.h"
#include "ShapeOneDCub.h"
#include "ShapeQuadrilateralLin.h"
#include "ShapePoint.h"
#include "ShapeTetrahedronLin.h"
#include "ShapeTetrahedronQua.h"
#include "ShapeTetrahedronCub.h"
#include "ShapeTriangleLin.h"
#include "ShapeTriangleQua.h"
#include "ShapeTriangleCub.h"
#include <metis.h>
#include<cstdlib>
#include<fstream>
#include<iostream>
#include "ElementT.h"
#include <vector>

void GmshTools::RenumberConnectivity(CompMesh *cmesh){
    // Renumber nodes - start
    std::vector<int > neighborNodes;

    int* xadj;
    int numNd = cmesh->NNodes();
    std::vector<int> adjncy;
    xadj = new int[numNd+1]();
    adjncy.reserve(10*numNd);

    for (int iNode = 0; iNode < numNd; iNode++){
        neighborNodes.reserve(cmesh->NodeVec()[iNode] -> getNumberOfElements()*3);
        neighborNodes.push_back(iNode);

        for (int j = 0; j < cmesh->NodeVec()[iNode] -> getNumberOfElements(); j++){
            int elem = cmesh->NodeVec()[iNode] -> getInverseIncidenceElement(j);
            if (!cmesh->ElementVec()[elem]) continue;
            VecInt connec = cmesh->ElementVec()[elem] -> getConnectivity();

            // std::cout << "COMM " << connec[0] << " " << connec[4] << std::endl;
            bool flag = false;
            for (int i = 0; i < connec.size(); i++){
                for (int iNeig = 0; iNeig < neighborNodes.size(); iNeig++){
                    if (connec[i] == neighborNodes[iNeig]){
                        flag = true;
                        break;
                    }
                }
                if (flag == false) neighborNodes.push_back(connec[i]);
                flag = false;
            }
        }
        //Save nodal adjacency for domain partitioning
        xadj[iNode+1] = xadj[iNode] + neighborNodes.size() - 1;
        for (int i = 1; i < neighborNodes.size(); i++){
            adjncy.push_back(neighborNodes[i]);
        }
        neighborNodes.clear();
        neighborNodes.shrink_to_fit();
    }

    //Save a second adjacency vector in idx_t format
    int *adjncy2;
    int adj_size = adjncy.size();
    adjncy2 = new int[adj_size];
    for (int i = 0; i < adjncy.size(); i++) adjncy2[i] = adjncy[i];
    
    adjncy.clear(); adjncy.shrink_to_fit();

    int* perm;
    int* iperm;
    perm = new int[numNd];
    iperm = new int[numNd];

    // Call METIS for node renumbering

    // MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
    //                         numNodes, numNodes,
    //                         1,NULL,1,NULL,&A); 
    
    // // MatGetOrdering(A, MATORDERINGMETISND, IS *rperm, IS *cperm)

    // MatDestroy(&A);

    METIS_NodeND(&numNd, xadj, adjncy2, NULL, NULL, perm, iperm);

    //Reorder nodes
    for (int i = 0, j; i < cmesh->NNodes(); ++i) {
        for (j = iperm[i]; j < i; j = iperm[j]);
        if (j == i) while (j = iperm[j],j != i) std::swap(cmesh->NodeVec()[i],cmesh->NodeVec()[j]);
    }

    // Update connectivity
    for (int i = 0; i < cmesh->NElements(); i++){
        if (!cmesh->ElementVec()[i]) continue;
        VecInt connect = cmesh->ElementVec()[i] -> getConnectivity();

        //Reorder connectivity
        for (int k = 0; k < connect.size(); k++) connect[k] = iperm[connect[k]];
        cmesh->ElementVec()[i] -> setConnectivity(connect);
    }
    
    for (int i = 0; i < cmesh->NNodes(); i++) cmesh->NodeVec()[i] -> clearInverseIncidence();

    for (int i = 0; i < cmesh->NElements(); i++){
        if (!cmesh->ElementVec()[i]) continue;
        VecInt connect = cmesh->ElementVec()[i] -> getConnectivity();

        for (int k = 0; k < connect.size(); k++) cmesh->NodeVec()[connect[k]] -> pushInverseIncidence(i);
    }

    PetscLogDouble bytes = 0;
    PetscMemoryGetCurrentUsage(&bytes);
    PetscPrintf(PETSC_COMM_WORLD,"Memory used-1 %g M\n",bytes/(1024*1024));
    
    for (int i = 0; i < cmesh->NNodes(); i++){
        for (int j = 0; j < cmesh->NodeVec()[i] -> getNumberOfElements(); j++){
            int elJ = cmesh->NodeVec()[i] -> getInverseIncidenceElement(j);
            for (int k = 0; k <cmesh-> NodeVec()[i] -> getNumberOfElements(); k++)
                cmesh->ElementVec()[elJ] -> pushNeighborElement(cmesh->NodeVec()[i] -> getInverseIncidenceElement(k)); 
        }
    }

    PetscMemoryGetCurrentUsage(&bytes);
    PetscPrintf(PETSC_COMM_WORLD,"Memory used00 %g M\n",bytes/(1024*1024));

    for (int i = 0; i < cmesh->NElements(); i++) cmesh->ElementVec()[i] -> sortEraseNeighborElements();

    PetscMemoryGetCurrentUsage(&bytes);
    PetscPrintf(PETSC_COMM_WORLD,"Memory used11 %g M\n",bytes/(1024*1024));        

    // for (int i = 0; i < NElements(); i++){
    //     if (rank == 0) std::cout << "Neighbor " << i << " ";
    //     for (int j = 0; j < ElementVec()[i] -> getNumberOfNeighborElements(); j++){
    //         std::cout << ElementVec()[i] -> getNeighborElement(j) << " ";
    //     }
    //     std::cout << std::endl;
    // }

    delete [] perm;
    delete [] iperm;
    delete [] adjncy2;
    delete [] xadj;
    // Renumber nodes - end

}

static void InsertElement(CompMesh &gmesh, int elindex, std::ifstream & line);
static std::string GetFileVersion(const std::string& file_name);
Element* InsertElement(CompMesh * gmesh, int & physical_identifier, int & el_type, int  el_identifier, VecInt & node_identifiers);

int GetNumberofNodes(int & el_type){
    int n_nodes;
    switch (el_type) {
        case 1:
        {   // Line
            n_nodes = 2;
        }
            break;
        case 2:
        {
            // Triangle
            n_nodes = 3;
        }
            break;
        case 3:
        {
            // Quadrilateral
            n_nodes = 4;
        }
            break;
        case 4:
        {
            // Tetrahedron
            n_nodes = 4;
        }
            break;
        case 5:
        {
            // Hexahedra
            n_nodes = 8;
        }
            break;
        case 6:
        {
            // Prism
            n_nodes = 6;
            std::cout << "Element not implemented\n";
            PanicButton();
        }
            break;
        case 7:
        {
            // Pyramid
            n_nodes = 5;
            std::cout << "Element not implemented\n";
            PanicButton();
        }
            break;
        case 8:
        {
            // Quadratic Line
            n_nodes = 3;
        }
            break;
        case 9:
        {
            // Quadratic Triangle
            n_nodes = 6;
        }
            break;
        case 10:
        {
            // Quadratic Quadrilateral
            n_nodes = 9;
        }
            break;
        case 11:
        {
            // Quadratic Tetrahedron
            n_nodes = 10;
            
        }
            break;
        case 12:
        {
            // Quadratic Hexahedra
            n_nodes = 20;
        }
            break;
        case 13:
        {
            // Quadratic Prism
            n_nodes = 15;
            std::cout << "Element not implemented\n";
            PanicButton();
        }
            break;
        case 15:{
            // Point
            n_nodes = 1;
        }
            break;
        case 21:{
            // Cubic Triangle
            n_nodes = 10;
        }
            break;
        case 26:{
            // Cubic Line
            n_nodes = 4;  
        }
            break;
        default:
        {
            std::cout << "Element not impelemented." << std::endl;
            n_nodes = 0;
            PanicButton();
        }
            break;
    }
    
    return n_nodes;
}





void GmshTools::Read3(CompMesh &gmesh, const std::string &file_name)
{
     /** @brief MaterialVec */
    /** Structure of both, physical entities dimension and names */
    std::vector<std::map<int,std::string> > fMaterialDataVec(4);
    // std::string string_temp;
    
    // //  Mesh Creation
    
    // {
        
    //     // reading a general mesh information by filter
    //     std::ifstream read (file_name);
        
    //     while(read)
    //     {
    //         char buf[1024];
    //         read.getline(buf, 1024);
    //         std::string str(buf);
            
    //         if(str == "$MeshFormat" || str == "$MeshFormat\r")
    //         {
    //             read.getline(buf, 1024);
    //             std::string str(buf);
    //             std::cout << "Reading mesh format = " << str << std::endl;
                
    //         }
            
    //         if(str == "$PhysicalNames" || str == "$PhysicalNames\r" )
    //         {
                
    //             int64_t n_entities;
    //             read >> n_entities;
    //             int max_dimension = 0;
                
    //             int dimension, id;
    //             std::string name;
    //             std::pair<int, std::string> chunk;
                
    //             for (int64_t inode = 0; inode < n_entities; inode++) {
                    
    //                 read.getline(buf, 1024);
    //                 read >> dimension;
    //                 read >> id;
    //                 read >> name;
    //                 name.erase(0,1);
    //                 name.erase(name.end()-1,name.end());
    //                 fMaterialDataVec[dimension][id] = name;
                    
    //                 if (max_dimension < dimension) {
    //                     max_dimension = dimension;
    //                 }
    //             }
                
    //             char buf_end[1024];
    //             read.getline(buf_end, 1024);
    //             read.getline(buf_end, 1024);
    //             std::string str_end(buf_end);
    //             if(str_end == "$EndPhysicalNames" || str_end == "$EndPhysicalNames\r")
    //             {
    //                 std::cout << "Read mesh physical entities = " << n_entities << std::endl;
    //             }
    //             continue;
    //         }
            
    //         if(str == "$Nodes" || str == "$Nodes\r")
    //         {
                
    //             int64_t n_nodes;
    //             read >> n_nodes;
                
    //             gmesh.SetNumNodes(n_nodes);
                
    //             // needed for node insertion
                
    //             for (int64_t inode = 0; inode < n_nodes; inode++) {
                    
    //                 VecDouble co(3);
    //                 int64_t node_id;
    //                 read.getline(buf, 1024);
    //                 read >> node_id;
    //                 read >> co[0];
    //                 read >> co[1];
    //                 read >> co[2];
                    
    //                 gmesh.Node(node_id-GMSH_SHIFT).SetCo(co);
                    
    //             }
                
                
    //             char buf_end[1024];
    //             read.getline(buf_end, 1024);
    //             read.getline(buf_end, 1024);
    //             std::string str_end(buf_end);
    //             if(str_end == "$EndNodes" || str_end == "$EndNodes\r")
    //             {
    //                 std::cout << "Read mesh nodes = " <<  gmesh.NumNodes() << std::endl;
    //             }
    //             continue;
    //         }
            
    //         if(str == "$Elements" || str == "$Elements\r")
    //         {
                
    //             int64_t n_elements;
    //             read >> n_elements;
    //             gmesh.SetNumElements(n_elements);
                
    //             for (int64_t iel = 0; iel < n_elements; iel++) {
    //                 InsertElement(gmesh, iel, read);
    //             }
                
    //             char buf_end[1024];
    //             read.getline(buf_end, 1024);
    //             read.getline(buf_end, 1024);
    //             std::string str_end(buf_end);
    //             if(str_end == "$EndElements" || str_end == "$EndElements\r")
    //             {
    //                 std::cout << "Read mesh elements = " << gmesh.NumElements() << std::endl;
    //             }
    //             continue;
    //         }
            
    //     }
        
    // }
    
    // std::cout << "Read General Mesh Data -> done!" << std::endl;
    // gmesh.BuildConnectivity();
    
    
    // std::cout << "Geometric Mesh Connectivity -> done!" << std::endl;

    
}

static std::string GetFileVersion(const std::string& file_name){

    std::ifstream read(file_name.c_str());
    if(!read)
    {
        std::cout << "Couldn't open the file " << file_name << std::endl;
        std::cout << "Maybe you have a wrong relative path. Check the current work directory you're running this executable from." << std::endl;
        PanicButton();
    }

    while(read){
        char buf[1024];
        read.getline(buf, 1024);
        std::string str(buf);
        
        if(str == "$MeshFormat" || str == "$MeshFormat\r")
        {
            read.getline(buf, 1024);
            std::string str(buf);
            return str;
        } // Mesh Format
    }

    std::cerr << "Could not determine .msh file format version." << std::endl;
    PanicButton();
    return "void";
}



void GmshTools::Read4(CompMesh &gmesh, const std::string &file_name){
    int max_dimension = 0;
     /** @brief MaterialVec */
    /** Structure of both, physical entities dimension and names */
    std::vector<std::map<int,std::string> > fMaterialDataVec(4);
    // vector of 4 positions
    // first index is the dimension
    // for each dimension there is a map
    // I guess that the first key of the map is the entity tag as it apears in gmesh
    // I guess that the second key of the map is the number as it appears in .gmsh file
    std::array<std::map<int, std::vector<int>>,4> m_dim_entity_tag_and_physical_tag;
    std::array<std::map<int, int>,4> m_dim_physical_tag_and_physical_tag;

    std::ifstream read(file_name.c_str());
    if(!read)
    {
        std::cout << "Couldn't open the file " << file_name << std::endl;
        std::cout << "Maybe you have a wrong relative path. Check the current work directory you're running this executable from." << std::endl;
        PanicButton();
    }

    while(read){
        char buf[1024];
        read.getline(buf, 1024);
        std::string str(buf);
        
        if(str == "$MeshFormat" || str == "$MeshFormat\r")
        {
            read.getline(buf, 1024);
            std::string str(buf);
            std::cout << "Reading mesh format = " << str << std::endl;
            if(str[0] != '4'){
                std::cerr << __PRETTY_FUNCTION__ << '\n'
                        << "Was intended for .msh format 4.1" << std::endl;
                PanicButton();
            }
            read.getline(buf, 1024); //<- Skip $EndMeshFormat
            continue;
        } // Mesh Format

        if(str == "$PhysicalNames" || str == "$PhysicalNames\r" )
        {
            
            int n_physical_names;
            read >> n_physical_names;

            int dimension, matid;
            std::string name;
            // std::pair<int, std::string> chunk;


            for (int64_t i_name = 0; i_name < n_physical_names; i_name++) {
                
                read.getline(buf, 1024);
                read >> dimension;
                read >> matid;
                read >> name;
                name.erase(0,1);
                name.erase(name.end()-1,name.end());
                fMaterialDataVec[dimension][matid] = name;
                
                
                if (max_dimension < dimension) {
                    max_dimension = dimension;
                }
            }
            
            char buf_end[1024];
            read.getline(buf_end, 1024);
            read.getline(buf_end, 1024);
            std::string str_end(buf_end);
            if(str_end == "$EndPhysicalNames" || str_end == "$EndPhysicalNames\r")
            {
                std::cout << "Read mesh number of physical names = " << n_physical_names << std::endl;
            }            
            continue;
        } // Physical Names

        // if(str == "$Entities" || str == "$Entities\r")
        // {
        //     std::string line;
        //     // We're ignoring Geometrical Entities for this code
        //     while(getline(read, line) && line != "$EndEntities"){/*void*/}
        //     continue;
        // }// Entities
        if(str == "$Entities" || str == "$Entities\r")
        {
            int m_n_points,m_n_curves,m_n_surfaces,m_n_volumes;
            read >> m_n_points;
            read >> m_n_curves;
            read >> m_n_surfaces;
            read >> m_n_volumes;

            max_dimension = (m_n_curves > 0 ?   1 : max_dimension);
            max_dimension = (m_n_surfaces > 0 ? 2 : max_dimension);
            max_dimension = (m_n_volumes > 0 ?  3 : max_dimension);
            
            int n_physical_tag;
            std::pair<int, std::vector<int> > chunk;
            /// Entity bounding box data
            double x_min, y_min, z_min;
            double x_max, y_max, z_max;
            std::vector<int> n_entities = {m_n_points,m_n_curves,m_n_surfaces,m_n_volumes};
            std::vector<int> n_entities_with_physical_tag = {0,0,0,0};
            
            
            for (int i_dim = 0; i_dim <4; i_dim++) {
                for (int64_t i_entity = 0; i_entity < n_entities[i_dim]; i_entity++) {
                    
                    read.getline(buf, 1024);
                    read >> chunk.first;
                    read >> x_min;
                    read >> y_min;
                    read >> z_min;
                    if(i_dim > 0)
                    {
                        read >> x_max;
                        read >> y_max;
                        read >> z_max;
                    }
                    read >> n_physical_tag;
                    if(n_physical_tag == 0){
                        chunk.second.resize(1,-999);
                    }else{
                        chunk.second.resize(n_physical_tag);
                        for (int i_data = 0; i_data < n_physical_tag; i_data++) {
                            read >> chunk.second[i_data];
                        }
                    }
                    if(i_dim > 0)
                    {
                        size_t n_bounding_points;
                        read >> n_bounding_points;
                        for (int i_data = 0; i_data < n_bounding_points; i_data++) {
                            int point_tag;
                            read >> point_tag;
                        }
                    }
                    n_entities_with_physical_tag[i_dim] += n_physical_tag;
                    m_dim_entity_tag_and_physical_tag[i_dim].insert(chunk);
                }
            }

            int m_n_physical_points = n_entities_with_physical_tag[0];
            int m_n_physical_curves = n_entities_with_physical_tag[1];
            int m_n_physical_surfaces = n_entities_with_physical_tag[2];
            int m_n_physical_volumes = n_entities_with_physical_tag[3];
            
            char buf_end[1024];
            read.getline(buf_end, 1024);
            read.getline(buf_end, 1024);
            std::string str_end(buf_end);
            if(str_end == "$EndEntities" || str_end == "$EndEntities\r")
            {
                std::cout << "Read mesh entities = " <<  m_n_points + m_n_curves + m_n_surfaces + m_n_volumes << std::endl;
                std::cout << "Read mesh entities with physical tags = " <<  m_n_physical_points + m_n_physical_curves + m_n_physical_surfaces + m_n_physical_volumes << std::endl;
            }
            continue;
        }

        if(str == "$Nodes" || str == "$Nodes\r")
        {
            
            int64_t n_entity_blocks, n_nodes, min_node_tag, max_node_tag;
            read >> n_entity_blocks;
            read >> n_nodes;
            read >> min_node_tag;
            read >> max_node_tag;
            
            int64_t node_id;
            VecDouble coord(3);
            gmesh.SetNumNodes(max_node_tag);
            
            int entity_tag, entity_dim, entity_parametric, entity_nodes;
            for (int64_t i_block = 0; i_block < n_entity_blocks; i_block++)
            {
                read.getline(buf, 1024);
                read >> entity_dim;
                read >> entity_tag;
                read >> entity_parametric;
                read >> entity_nodes;
                
                if (entity_parametric != 0) {
                    std::cout << "GmshTools:: Characteristic not implemented." << std::endl;
                    PanicButton();
                }
                
                std::vector<int64_t> nodeids(entity_nodes,-1);
                for (int64_t inode = 0; inode < entity_nodes; inode++) {
                    read >> nodeids[inode];
                }
                for (int64_t inode = 0; inode < entity_nodes; inode++) {
                    read >> coord[0];
                    read >> coord[1];
                    read >> coord[2];
                    
                    Node *node = new Node(coord,nodeids[inode]-GMSH_SHIFT,gmesh.NState());
                    gmesh.NodeVec()[nodeids[inode]-GMSH_SHIFT] = node;
                    // gmesh.Node(nodeids[inode] - GMSH_SHIFT).SetCo(coord);                    
                }
            }
            
            char buf_end[1024];
            read.getline(buf_end, 1024);
            read.getline(buf_end, 1024);
            std::string str_end(buf_end);
            if(str_end == "$EndNodes" || str_end == "$EndNodes\r")
            {
                std::cout << "Read mesh nodes = " <<  gmesh.NNodes() << std::endl;
            }
            continue;
        } // Nodes


        if(str == "$Elements" || str == "$Elements\r")
        {
            
            int64_t n_entity_blocks, n_elements, min_element_tag, max_element_tag;
            read >> n_entity_blocks;
            read >> n_elements;
            read >> min_element_tag;
            read >> max_element_tag;
            gmesh.SetNumElements(n_elements);
            int64_t elcount = 0;
            
            int entity_tag, entity_dim, entity_el_type, entity_elements;
            for (int64_t i_block = 0; i_block < n_entity_blocks; i_block++)
            {
                read.getline(buf, 1024);
                read >> entity_dim;
                read >> entity_tag;
                read >> entity_el_type;
                read >> entity_elements;
                
                if(entity_elements == 0){
                    std::cout << "The entity with tag " << entity_tag << " does not have elements to insert" << std::endl;
                }
                
                for (int64_t iel = 0; iel < entity_elements; iel++) {
                    int physical_identifier;
                    int n_physical_identifier = 0;
                    if(m_dim_entity_tag_and_physical_tag[entity_dim].find(entity_tag) != m_dim_entity_tag_and_physical_tag[entity_dim].end())
                    {
                        n_physical_identifier = m_dim_entity_tag_and_physical_tag[entity_dim][entity_tag].size();
                    }
                    bool physical_identifier_Q = n_physical_identifier != 0;
                    if(physical_identifier_Q)
                    {
                        int gmsh_physical_identifier = m_dim_entity_tag_and_physical_tag[entity_dim][entity_tag][0];
                        physical_identifier = m_dim_physical_tag_and_physical_tag[entity_dim][gmsh_physical_identifier];
                        if(n_physical_identifier !=1){
                            std::cout << "The entity with tag " << entity_tag << std::endl;
                            std::cout << "Has associated the following physical tags : " << std::endl;
                            for (int i_data = 0; i_data < n_physical_identifier; i_data++) {
                                std::cout << m_dim_entity_tag_and_physical_tag[entity_dim][entity_tag][i_data] << std::endl;
                            }
                            
                            std::cout << "Automatically, the assigned pz physical tag = " << physical_identifier << " is used.  The other ones are dropped out." << std::endl;
                        }
                        
                        
                        read.getline(buf, 1024);
                        int el_identifier, n_el_nodes;
                        n_el_nodes = GetNumberofNodes(entity_el_type);
                        read >> el_identifier;
                        VecInt node_identifiers(n_el_nodes);
                        for (int i_node = 0; i_node < n_el_nodes; i_node++) {
                            read >> node_identifiers[i_node];
                        }
                        /// Internally the nodes index and element index is converted to zero based indexation
                        gmesh.ElementVec()[elcount] = InsertElement(&gmesh, gmsh_physical_identifier, entity_el_type, elcount, node_identifiers);
                        elcount++;
                        
                    }else{
                        read.getline(buf, 1024);
                        int el_identifier, n_el_nodes;
                        n_el_nodes = GetNumberofNodes(entity_el_type);
                        read >> el_identifier;
                        VecInt node_identifiers(n_el_nodes);
                        for (int i_node = 0; i_node < n_el_nodes; i_node++) {
                            read >> node_identifiers[i_node];
                        }
                        std::cout << "The entity with tag " << entity_tag << " does not have a physical tag, element " << el_identifier << " skipped " << std::endl;
                    }

                }
            }
            
            char buf_end[1024];
            read.getline(buf_end, 1024);
            read.getline(buf_end, 1024);
            std::string str_end(buf_end);
            if(str_end == "$EndElements" || str_end == "$EndElements\r")
            {
                std::cout << "Read mesh elements = " << gmesh.NElements() << std::endl;
            }
            continue;
        } // Elements

    }

    gmesh.SetDimension(max_dimension);
    // gmesh.BuildConnectivity();
}

/** @brief Insert elements following msh file format */
void InsertElement(CompMesh &gmesh, int elindex, std::ifstream & line){
    
    // first implementation based on linear elements: http://gmsh.info/doc/texinfo/gmsh.html#File-formats
    VecInt TopolPoint(1);
    VecInt TopolLine(2);
    VecInt TopolTriangle(3);
    VecInt TopolQuad(4);
    VecInt TopolTet(4);
    VecInt TopolPyr(5);
    VecInt TopolPrism(6);
    VecInt TopolHex(8);
    
    VecInt TopolLineQ(3);
    VecInt TopolTriangleQ(6);
    VecInt TopolQuadQ(8);
    VecInt TopolTetQ(10);
    VecInt TopolPyrQ(14);
    VecInt TopolPrismQ(15);
    VecInt TopolHexQ(20);
    
    
    int64_t element_id, type_id, div_id, physical_id, elementary_id;
    
    
    char buf[1024];
    buf[0] = 0;
    line.getline(buf, 1024);
    
    line >> element_id;
    line >> type_id;
    line >> div_id;
    line >> physical_id;
    line >> elementary_id;
    
    if (div_id != 2 || type_id > 15) {
        std::cout << "div_id " << div_id << " type_id " << type_id << std::endl;
        PanicButton();
    }
    int matid = physical_id;
    
    // switch (type_id) {
    //     case 1:
    //     {
    //         // Line
    //         line >> TopolLine[0]; //node 1
    //         line >> TopolLine[1]; //node 2
    //         element_id--;
    //         TopolLine[0]--;
    //         TopolLine[1]--;
    //         GeoElement *gel = new GeoElementTemplate< Geom1d> (TopolLine, matid, &gmesh,elindex);
    //         gmesh.SetElement(element_id, gel);
            
    //     }
    //         break;
    //     case 2:
    //     {
    //         // Triangle
    //         line >> TopolTriangle[0]; //node 1
    //         line >> TopolTriangle[1]; //node 2
    //         line >> TopolTriangle[2]; //node 3
    //         element_id--;
    //         TopolTriangle[0]--;
    //         TopolTriangle[1]--;
    //         TopolTriangle[2]--;
    //         GeoElement *gel = new GeoElementTemplate< GeomTriangle> (TopolTriangle, matid, &gmesh,elindex);
    //         gmesh.SetElement(element_id, gel);

    //     }
    //         break;
    //     case 3:
    //     {
    //         // Quadrilateral
    //         line >> TopolQuad[0]; //node 1
    //         line >> TopolQuad[1]; //node 2
    //         line >> TopolQuad[2]; //node 3
    //         line >> TopolQuad[3]; //node 4
    //         element_id--;
    //         TopolQuad[0]--;
    //         TopolQuad[1]--;
    //         TopolQuad[2]--;
    //         TopolQuad[3]--;
    //         GeoElement *gel = new GeoElementTemplate< GeomQuad> (TopolQuad, matid, &gmesh,elindex);
    //         gmesh.SetElement(element_id, gel);

    //     }
    //         break;
    //     case 4:
    //     {
    //         // Tetrahedron
    //         line >> TopolTet[0]; //node 1
    //         line >> TopolTet[1]; //node 2
    //         line >> TopolTet[2]; //node 3
    //         line >> TopolTet[3]; //node 4
    //         element_id--;
    //         TopolTet[0]--;
    //         TopolTet[1]--;
    //         TopolTet[2]--;
    //         TopolTet[3]--;
    //         GeoElement * gel = new GeoElementTemplate< GeomTetrahedron> (TopolTet, matid, &gmesh,elindex);
    //         gmesh.SetElement(element_id, gel);

    //     }
    //         break;
    //     case 5:
    //     {
    //         // Hexahedra
    //         line >> TopolHex[0]; //node 1
    //         line >> TopolHex[1]; //node 2
    //         line >> TopolHex[2]; //node 3
    //         line >> TopolHex[3]; //node 4
    //         line >> TopolHex[4]; //node 5
    //         line >> TopolHex[5]; //node 6
    //         line >> TopolHex[6]; //node 7
    //         line >> TopolHex[7]; //node 8
    //         element_id--;
    //         TopolHex[0]--;
    //         TopolHex[1]--;
    //         TopolHex[2]--;
    //         TopolHex[3]--;
    //         TopolHex[4]--;
    //         TopolHex[5]--;
    //         TopolHex[6]--;
    //         TopolHex[7]--;
    //         PanicButton();
    //     }
    //         break;
    //     case 6:
    //     {
    //         // Prism
    //         line >> TopolPrism[0]; //node 1
    //         line >> TopolPrism[1]; //node 2
    //         line >> TopolPrism[2]; //node 3
    //         line >> TopolPrism[3]; //node 4
    //         line >> TopolPrism[4]; //node 5
    //         line >> TopolPrism[5]; //node 6
    //         element_id--;
    //         TopolPrism[0]--;
    //         TopolPrism[1]--;
    //         TopolPrism[2]--;
    //         TopolPrism[3]--;
    //         TopolPrism[4]--;
    //         TopolPrism[5]--;
    //         PanicButton();
    //     }
    //         break;
    //     case 7:
    //     {
    //         // Pyramid
    //         line >> TopolPyr[0]; //node 1
    //         line >> TopolPyr[1]; //node 2
    //         line >> TopolPyr[2]; //node 3
    //         line >> TopolPyr[3]; //node 4
    //         line >> TopolPyr[4]; //node 5
    //         element_id--;
    //         TopolPyr[0]--;
    //         TopolPyr[1]--;
    //         TopolPyr[2]--;
    //         TopolPyr[3]--;
    //         TopolPyr[4]--;
    //         PanicButton();
    //     }
    //         break;
    //     case 8:
    //     {
    //         // Quadratic line
    //         line >> TopolLineQ[0]; //node 1
    //         line >> TopolLineQ[1]; //node 2
    //         line >> TopolLineQ[2]; //node 2
    //         element_id--;
    //         TopolLineQ[0]--;
    //         TopolLineQ[1]--;
    //         TopolLineQ[2]--;
    //         PanicButton();
    //     }
    //         break;
    //     case 9:
    //     {
    //         // Triangle
    //         line >> TopolTriangleQ[0]; //node 1
    //         line >> TopolTriangleQ[1]; //node 2
    //         line >> TopolTriangleQ[2]; //node 3
    //         line >> TopolTriangleQ[3]; //node 4
    //         line >> TopolTriangleQ[4]; //node 5
    //         line >> TopolTriangleQ[5]; //node 6
    //         element_id--;
    //         TopolTriangleQ[0]--;
    //         TopolTriangleQ[1]--;
    //         TopolTriangleQ[2]--;
    //         TopolTriangleQ[3]--;
    //         TopolTriangleQ[4]--;
    //         TopolTriangleQ[5]--;
    //         PanicButton();
    //     }
    //         break;
    //     case 10:
    //     {
    //         // Quadrilateral
    //         line >> TopolQuadQ[0]; //node 1
    //         line >> TopolQuadQ[1]; //node 2
    //         line >> TopolQuadQ[2]; //node 3
    //         line >> TopolQuadQ[3]; //node 4
    //         line >> TopolQuadQ[4]; //node 5
    //         line >> TopolQuadQ[5]; //node 6
    //         line >> TopolQuadQ[6]; //node 7
    //         line >> TopolQuadQ[7]; //node 8
    //         element_id--;
    //         TopolQuadQ[0]--;
    //         TopolQuadQ[1]--;
    //         TopolQuadQ[2]--;
    //         TopolQuadQ[3]--;
    //         TopolQuadQ[4]--;
    //         TopolQuadQ[5]--;
    //         TopolQuadQ[6]--;
    //         TopolQuadQ[7]--;
    //         PanicButton();
    //     }
    //         break;
    //     case 11:
    //     {
    //         // Tetrahedron
    //         line >> TopolTetQ[0]; //node 1
    //         line >> TopolTetQ[1]; //node 2
    //         line >> TopolTetQ[2]; //node 3
    //         line >> TopolTetQ[3]; //node 4
            
    //         line >> TopolTetQ[4]; //node 4 -> 4
    //         line >> TopolTetQ[5]; //node 5 -> 5
    //         line >> TopolTetQ[6]; //node 6 -> 6
    //         line >> TopolTetQ[7]; //node 7 -> 7
    //         line >> TopolTetQ[9]; //node 9 -> 8
    //         line >> TopolTetQ[8]; //node 8 -> 9
            
    //         element_id--;
    //         TopolTetQ[0]--;
    //         TopolTetQ[1]--;
    //         TopolTetQ[2]--;
    //         TopolTetQ[3]--;
            
    //         TopolTetQ[4]--;
    //         TopolTetQ[5]--;
    //         TopolTetQ[6]--;
    //         TopolTetQ[7]--;
    //         TopolTetQ[8]--;
    //         TopolTetQ[9]--;
    //         PanicButton();
            
    //     }
    //         break;
    //     case 12:
    //     {
    //         // Hexahedra
    //         line >> TopolHexQ[0]; //node 1
    //         line >> TopolHexQ[1]; //node 2
    //         line >> TopolHexQ[2]; //node 3
    //         line >> TopolHexQ[3]; //node 4
    //         line >> TopolHexQ[4]; //node 5
    //         line >> TopolHexQ[5]; //node 6
    //         line >> TopolHexQ[6]; //node 7
    //         line >> TopolHexQ[7]; //node 8
            
    //         line >> TopolHexQ[8];  //node 8  -> 8
    //         line >> TopolHexQ[11]; //node 11 -> 9
    //         line >> TopolHexQ[12]; //node 12 -> 10
    //         line >> TopolHexQ[9];  //node 9  -> 11
    //         line >> TopolHexQ[13]; //node 13 -> 12
    //         line >> TopolHexQ[10]; //node 10 -> 13
    //         line >> TopolHexQ[14]; //node 14 -> 14
    //         line >> TopolHexQ[15]; //node 15 -> 15
    //         line >> TopolHexQ[16]; //node 16 -> 16
    //         line >> TopolHexQ[19]; //node 19 -> 17
    //         line >> TopolHexQ[17]; //node 17 -> 18
    //         line >> TopolHexQ[18]; //node 18 -> 19
            
    //         element_id--;
    //         TopolHexQ[0]--;
    //         TopolHexQ[1]--;
    //         TopolHexQ[2]--;
    //         TopolHexQ[3]--;
    //         TopolHexQ[4]--;
    //         TopolHexQ[5]--;
    //         TopolHexQ[6]--;
    //         TopolHexQ[7]--;
            
    //         TopolHexQ[8]--;
    //         TopolHexQ[9]--;
    //         TopolHexQ[10]--;
    //         TopolHexQ[11]--;
    //         TopolHexQ[12]--;
    //         TopolHexQ[13]--;
    //         TopolHexQ[14]--;
    //         TopolHexQ[15]--;
    //         TopolHexQ[16]--;
    //         TopolHexQ[17]--;
    //         TopolHexQ[18]--;
    //         TopolHexQ[19]--;
    //         PanicButton();
    //     }
    //         break;
    //     case 13:
    //     {
    //         // Prism
    //         line >> TopolPrismQ[0]; //node 1
    //         line >> TopolPrismQ[1]; //node 2
    //         line >> TopolPrismQ[2]; //node 3
    //         line >> TopolPrismQ[3]; //node 4
    //         line >> TopolPrismQ[4]; //node 5
    //         line >> TopolPrismQ[5]; //node 6
            
    //         line >> TopolPrismQ[6];  //node 6 -> 6
    //         line >> TopolPrismQ[8];  //node 8 -> 7
    //         line >> TopolPrismQ[9];  //node 9 -> 8
    //         line >> TopolPrismQ[7];  //node 7 -> 9
    //         line >> TopolPrismQ[10]; //node 11 -> 10
    //         line >> TopolPrismQ[11]; //node 12 -> 11
    //         line >> TopolPrismQ[12]; //node 13 -> 12
    //         line >> TopolPrismQ[14]; //node 14 -> 13
    //         line >> TopolPrismQ[13]; //node 15 -> 14
            
    //         element_id--;
    //         TopolPrismQ[0]--;
    //         TopolPrismQ[1]--;
    //         TopolPrismQ[2]--;
    //         TopolPrismQ[3]--;
    //         TopolPrismQ[4]--;
    //         TopolPrismQ[5]--;
            
    //         TopolPrismQ[6]--;
    //         TopolPrismQ[7]--;
    //         TopolPrismQ[8]--;
    //         TopolPrismQ[9]--;
    //         TopolPrismQ[10]--;
    //         TopolPrismQ[11]--;
    //         TopolPrismQ[12]--;
    //         TopolPrismQ[13]--;
    //         TopolPrismQ[14]--;
            
    //         PanicButton();
    //     }
    //         break;
    //     case 15:
    //         line >> TopolPoint[0];
    //         TopolPoint[0]--;
    //         element_id--;
    //         std::cout << "PLEASE IMPLEMENT ME\n";
    //         PanicButton();
    //         break;
    //     default:
    //     {
    //         std::cout << "Element not impelemented." << std::endl;
    //         PanicButton();
    //     }
    //         break;
    // }
}

void GmshTools::Read(CompMesh& gmesh, const std::string& file_name){
    
    std::string format_version = GetFileVersion(file_name);

    switch(format_version[0]){
        case '3': Read3(gmesh,file_name); break;
        case '4': Read4(gmesh,file_name); break;
        default: 
            std::cout << "GmshTools:: Latest version supported 4.1 \n"
                      << "GmshTools:: Reader not available for the msh file version = " << format_version << '\n'
                      << "GmshTools:: Gmsh can probably export meshes in different legacy versions. Check their documentation for an up-to-date tutorial." << std::endl;
            PanicButton();
    }

    //Delete all null pointers, i.e., elements without a material
    // gmesh.ElementVec().erase(
    //     std::remove(gmesh.ElementVec().begin(), gmesh.ElementVec().end(), nullptr),
    //     gmesh.ElementVec().end()
    // );
    // gmesh.ElementVec().shrink_to_fit();

    RenumberConnectivity(&gmesh);

    gmesh.part_elem= new int[gmesh.NElements()]();
}





Element* InsertElement(CompMesh * gmesh, int & physical_identifier, int & el_type, int  el_identifier, VecInt & node_identifiers){
    
    VecInt Topology;
    int n_nodes = node_identifiers.size();
    Topology.resize(n_nodes);
    Topology.setConstant(-1);
    for (int k_node = 0; k_node<n_nodes; k_node++) {
        Topology[k_node] = node_identifiers[k_node]-GMSH_SHIFT;
    }
    Element* gel = nullptr;
//    el_identifier -= GMSH_SHIFT;

    for (int k = 0; k < node_identifiers.size(); k++){
        gmesh->NodeVec()[Topology[k]] -> pushInverseIncidence(el_identifier);
    };

    //The element should be ignored if there is no material inserted with its physical identifier
    if (!gmesh->Material(physical_identifier)) return gel;

    switch (el_type) {
        case 1:
        {   // Ligelne
        
            gel = new ElementT<ShapeOneDLin>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
            gel->PrintType() = 3;
        }
            break;
        case 2:
        {
            // Linear Triangle
            gel = new ElementT<ShapeTriangleLin>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
            gel->PrintType() = 5;           
            break;
        }
        case 8:
        {
            // Quadratic Line
            gel = new ElementT<ShapeOneDQua>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
            gel->PrintType() = 21;           
            break;
        }
        case 9:
        {
            // Quadratic Triangle
            gel = new ElementT<ShapeTriangleQua>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
            gel->PrintType() = 22;           
            break;
        }
        case 15:
        {
            // Point
            gel = new ElementT<ShapePoint>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
            gel->PrintType() = 1;           
            break;
        }
        case 21:
        {
            // Cubic Triangle
            gel = new ElementT<ShapeTriangleCub>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
            gel->PrintType() = 69;           
            break;
        }
        case 26:
        {
            // Cubic Line
            gel = new ElementT<ShapeOneDCub>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
            gel->PrintType() = 35;           
            break;
        }
        
        default:
        {
            std::cout << "Element not implemented." << std::endl;
            PanicButton();
        }
            break;
    }
    // ElementT<tshape> *gel = new ElementT<tshape>();
    // gmesh->SetElement(el_identifier, gel);
    return gel;
};
