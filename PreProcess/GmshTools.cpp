#include "GmshTools.h"
#include "ShapeHexahedron.h"
#include "ShapeOneDLin.h"
#include "ShapeOneDQua.h"
#include "ShapeOneDCub.h"
#include "ShapeQuadrilateralLin.h"
#include "ShapeQuadrilateralQua.h"
#include "ShapePoint.h"
#include "ShapeTetrahedronLin.h"
#include "ShapeTetrahedronQua.h"
#include "ShapeTetrahedronCub.h"
#include "ShapeTriangleLin.h"
#include "ShapeTriangleQua.h"
#include "ShapeTriangleCub.h"
#include "HierarquicalOneD.h"
#include "HierarquicalQuad.h"
#include "HierarquicalTriangle.h"
#ifdef HAS_METIS
#include <metis.h>
#endif
#include<cstdlib>
#include<fstream>
#include<iostream>
#include "ElementT.h"
#include "ElementWithMem.h"
#include "ElementTransient.h"
#include "PlasticityModel.h"
#include "TransientWeakForm.h"
#include "GeoElement.h"
#include "GeoElementT.h"
#include <vector>

void GmshTools::RenumberConnectivity(GeoMesh *cmesh){
    // Renumber nodes - start
#ifdef HAS_METIS
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
            VecInt connec = cmesh->ElementVec()[elem] -> getGeometricNodes();

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
        VecInt connect = cmesh->ElementVec()[i] -> getGeometricNodes();

        //Reorder connectivity
        for (int k = 0; k < connect.size(); k++) connect[k] = iperm[connect[k]];
        cmesh->ElementVec()[i] -> setGeometricNodes(connect);
    }
    
    for (int i = 0; i < cmesh->NNodes(); i++) cmesh->NodeVec()[i] -> clearInverseIncidence();

    for (int i = 0; i < cmesh->NElements(); i++){
        if (!cmesh->ElementVec()[i]) continue;
        VecInt connect = cmesh->ElementVec()[i] -> getGeometricNodes();

        for (int k = 0; k < connect.size(); k++) cmesh->NodeVec()[connect[k]] -> pushInverseIncidence(i);
    }

#ifdef HAS_PETSC
    PetscLogDouble bytes = 0;
    PetscMemoryGetCurrentUsage(&bytes);
    PetscPrintf(PETSC_COMM_WORLD,"Memory used-1 %g M\n",bytes/(1024*1024));
#endif

    for (int i = 0; i < cmesh->NNodes(); i++){
        for (int j = 0; j < cmesh->NodeVec()[i] -> getNumberOfElements(); j++){
            int elJ = cmesh->NodeVec()[i] -> getInverseIncidenceElement(j);
            for (int k = 0; k <cmesh-> NodeVec()[i] -> getNumberOfElements(); k++)
                cmesh->ElementVec()[elJ] -> pushNeighborElement(cmesh->NodeVec()[i] -> getInverseIncidenceElement(k)); 
        }
    }

#ifdef HAS_PETSC
    PetscMemoryGetCurrentUsage(&bytes);
    PetscPrintf(PETSC_COMM_WORLD,"Memory used00 %g M\n",bytes/(1024*1024));
#endif

    for (int i = 0; i < cmesh->NElements(); i++) cmesh->ElementVec()[i] -> sortEraseNeighborElements();

#ifdef HAS_PETSC
    PetscMemoryGetCurrentUsage(&bytes);
    PetscPrintf(PETSC_COMM_WORLD,"Memory used11 %g M\n",bytes/(1024*1024));        
#endif
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
#endif

}

void GmshTools::BuildNeighbourInformation(GeoMesh *cmesh){
    // Renumber nodes - start

    std::vector<int > neighborNodes;
    int numNd = cmesh->NNodes();
    
    for (int iNode = 0; iNode < numNd; iNode++){
        neighborNodes.reserve(cmesh->NodeVec()[iNode] -> getNumberOfElements()*3);
        neighborNodes.push_back(iNode);

        for (int j = 0; j < cmesh->NodeVec()[iNode] -> getNumberOfElements(); j++){
            int elem = cmesh->NodeVec()[iNode] -> getInverseIncidenceElement(j);
            if (!cmesh->ElementVec()[elem]) continue;
            VecInt connec = cmesh->ElementVec()[elem] -> getGeometricNodes();

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
        neighborNodes.clear();
        neighborNodes.shrink_to_fit();
    }
    
    for (int i = 0; i < cmesh->NNodes(); i++) cmesh->NodeVec()[i] -> clearInverseIncidence();

    for (int i = 0; i < cmesh->NElements(); i++){
        if (!cmesh->ElementVec()[i]) continue;
        VecInt connect = cmesh->ElementVec()[i] -> getGeometricNodes();

        for (int k = 0; k < connect.size(); k++) cmesh->NodeVec()[connect[k]] -> pushInverseIncidence(i);
    }

    for (int i = 0; i < cmesh->NNodes(); i++){
        for (int j = 0; j < cmesh->NodeVec()[i] -> getNumberOfElements(); j++){
            int elJ = cmesh->NodeVec()[i] -> getInverseIncidenceElement(j);
            for (int k = 0; k <cmesh-> NodeVec()[i] -> getNumberOfElements(); k++)
                cmesh->ElementVec()[elJ] -> pushNeighborElement(cmesh->NodeVec()[i] -> getInverseIncidenceElement(k)); 
        }
    }

    for (int i = 0; i < cmesh->NElements(); i++) cmesh->ElementVec()[i] -> sortEraseNeighborElements();

}

static std::string GetFileVersion(const std::string& file_name);
GeoElement* InsertElement(GeoMesh * gmesh, int & physical_identifier, int & el_type, int  el_identifier, VecInt & node_identifiers);

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
        case 36:{
            // Cubic Quadrilateral
            n_nodes = 16;  
        }
            break;
        case 29:{
            // Cubic Tetrahedron
            n_nodes = 20;  
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



void GmshTools::Read4(GeoMesh &gmesh, const std::string &file_name){
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


            for (int i_name = 0; i_name < n_physical_names; i_name++) {
                
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
            REAL x_min, y_min, z_min;
            REAL x_max, y_max, z_max;
            std::vector<int> n_entities = {m_n_points,m_n_curves,m_n_surfaces,m_n_volumes};
            std::vector<int> n_entities_with_physical_tag = {0,0,0,0};
            
            
            for (int i_dim = 0; i_dim <4; i_dim++) {
                for (int i_entity = 0; i_entity < n_entities[i_dim]; i_entity++) {
                    
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
            
            int n_entity_blocks, n_nodes, min_node_tag, max_node_tag;
            read >> n_entity_blocks;
            read >> n_nodes;
            read >> min_node_tag;
            read >> max_node_tag;
            
            int node_id;
            VecDouble coord(3);
            gmesh.SetNumNodes(max_node_tag);
            
            int entity_tag, entity_dim, entity_parametric, entity_nodes;
            for (int i_block = 0; i_block < n_entity_blocks; i_block++)
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
                
                std::vector<int> nodeids(entity_nodes,-1);
                for (int inode = 0; inode < entity_nodes; inode++) {
                    read >> nodeids[inode];
                }
                for (int inode = 0; inode < entity_nodes; inode++) {
                    read >> coord[0];
                    read >> coord[1];
                    read >> coord[2];
                    
                    Node *node = new Node(coord,nodeids[inode]-GMSH_SHIFT);
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
            
            int n_entity_blocks, n_elements, min_element_tag, max_element_tag;
            read >> n_entity_blocks;
            read >> n_elements;
            read >> min_element_tag;
            read >> max_element_tag;
            gmesh.SetNumElements(n_elements);
            int elcount = 0;
            
            int entity_tag, entity_dim, entity_el_type, entity_elements;
            for (int i_block = 0; i_block < n_entity_blocks; i_block++)
            {
                read.getline(buf, 1024);
                read >> entity_dim;
                read >> entity_tag;
                read >> entity_el_type;
                read >> entity_elements;
                
                if(entity_elements == 0){
                    std::cout << "The entity with tag " << entity_tag << " does not have elements to insert" << std::endl;
                }
                
                for (int iel = 0; iel < entity_elements; iel++) {
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

void GmshTools::Read(GeoMesh& gmesh, const std::string& file_name){
    
    std::string format_version = GetFileVersion(file_name);

    switch(format_version[0]){
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
#ifdef HAS_METIS
    RenumberConnectivity(&gmesh);
#else
    BuildNeighbourInformation(&gmesh);
#endif

    // gmesh.part_elem= new int[gmesh.NElements()]();
}

//Inserts elements onto the mesh
GeoElement* InsertElement(GeoMesh * gmesh, int & physical_identifier, int & el_type, int  el_identifier, VecInt & node_identifiers){
    
    VecInt Topology;
    int n_nodes = node_identifiers.size();
    Topology.resize(n_nodes);
    Topology.setConstant(-1);
    for (int k_node = 0; k_node<n_nodes; k_node++) {
        Topology[k_node] = node_identifiers[k_node]-GMSH_SHIFT;
    }
    GeoElement* gel = nullptr;
//    el_identifier -= GMSH_SHIFT;

    for (int k = 0; k < node_identifiers.size(); k++){
        gmesh->NodeVec()[Topology[k]] -> pushInverseIncidence(el_identifier);
    };

    switch (el_type) {
        case 1:
        {   // Line
            gel = new GeoElementT<ShapeOneDLin>(el_identifier,Topology,gmesh,physical_identifier);
            gel->PrintType() = 3;
        }
            break;
        case 2:
        {   // Linear Triangle
            gel = new GeoElementT<ShapeTriangleLin>(el_identifier,Topology,gmesh,physical_identifier);
            gel->PrintType() = 5;           
            break;
        }
        case 3:
        {   // Linear Quadrilateral
            gel = new GeoElementT<ShapeQuadrilateralLin>(el_identifier,Topology,gmesh,physical_identifier);
            gel->PrintType() = 9;
            break;
        }
        case 4:
        {   // Linear Tetrahedron
            gel = new GeoElementT<ShapeTetrahedronLin>(el_identifier,Topology,gmesh,physical_identifier);
            gel->PrintType() = 10;
            break;
        }
        case 5:
        {   // Linear Hexahedron
            gel = new GeoElementT<ShapeHexahedron>(el_identifier,Topology,gmesh,physical_identifier);
            gel->PrintType() = 12;
            break;
        }
        case 8:
        {   // Quadratic Line
            gel = new GeoElementT<ShapeOneDQua>(el_identifier,Topology,gmesh,physical_identifier);
            gel->PrintType() = 21;           
            break;
        }
        case 9:
        {   // Quadratic Triangle
            gel = new GeoElementT<ShapeTriangleQua>(el_identifier,Topology,gmesh,physical_identifier);
            gel->PrintType() = 22;           
            break;
        }
        case 10:
        {   // Quadratic Quad
            gel = new GeoElementT<ShapeQuadrilateralQua>(el_identifier,Topology,gmesh,physical_identifier);
            gel->PrintType() = 28;           
            break;
        }
        case 11:
        {   // Quadratic Tetrahedron
            VecInt Topology2(10);
            Topology2[0] = Topology[0];
            Topology2[1] = Topology[3];
            Topology2[2] = Topology[1];
            Topology2[3] = Topology[2];
            Topology2[4] = Topology[7];
            Topology2[5] = Topology[9];
            Topology2[6] = Topology[4];
            Topology2[7] = Topology[6];
            Topology2[8] = Topology[8];
            Topology2[9] = Topology[5];
            gel = new GeoElementT<ShapeTetrahedronQua>(el_identifier,Topology2,gmesh,physical_identifier);
            gel->PrintType() = 24;
            break;
        }
        case 15:
        {   // Point
            gel = new GeoElementT<ShapePoint>(el_identifier,Topology,gmesh,physical_identifier);
            gel->PrintType() = 1;           
            break;
        }
        case 21:
        {   // Cubic Triangle
            gel = new GeoElementT<ShapeTriangleCub>(el_identifier,Topology,gmesh,physical_identifier);
            gel->PrintType() = 69;           
            break;
        }
        case 26:
        {   // Cubic Line
            gel = new GeoElementT<ShapeOneDCub>(el_identifier,Topology,gmesh,physical_identifier);
            gel->PrintType() = 35;           
            break;
        }
        case 29:
        {   // Cubic Tetrahedron
            VecInt Topology2(20);
            Topology2[0] = Topology[0];
            Topology2[1] = Topology[1];
            Topology2[2] = Topology[2];
            Topology2[3] = Topology[3];
            Topology2[4] = Topology[12];
            Topology2[5] = Topology[13];
            Topology2[6] = Topology[14];
            Topology2[7] = Topology[15];
            Topology2[8] = Topology[10];
            Topology2[9] = Topology[11];
            Topology2[10] = Topology[8];
            Topology2[11] = Topology[9];
            Topology2[12] = Topology[6];
            Topology2[13] = Topology[7];
            Topology2[14] = Topology[4];
            Topology2[15] = Topology[5];
            Topology2[16] = Topology[19];
            Topology2[17] = Topology[18];
            Topology2[18] = Topology[17];
            Topology2[19] = Topology[16];
            gel = new GeoElementT<ShapeTetrahedronCub>(el_identifier,Topology2,gmesh,physical_identifier);
            gel->PrintType() = 71;
            break;
        }

        case 36:
        {
            PanicButton();
            // Cubic Quadrilateral
            // if (plasticmaterial){
            //     gel = new ElementWithMem<ShapeQuadrilateralCub>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
            // } else {
            //     gel = new ElementT<ShapeQuadrilateralCub>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
            // }
            // gel->PrintType() = 70;
            break;
        }
        
        default:
        {
            std::cout << "Element not implemented." << std::endl;
            PanicButton();
        }
            break;
    }
    return gel;
};







///*only for backup */
/*
//Inserts elements onto the mesh
Element* InsertElement(GeoMesh * gmesh, int & physical_identifier, int & el_type, int  el_identifier, VecInt & node_identifiers){
    
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
    PlasticityModel *plasticmaterial = dynamic_cast<PlasticityModel * > (gmesh->Material(physical_identifier));
    TransientWeakForm *transientmaterial = dynamic_cast<TransientWeakForm * > (gmesh->Material(physical_identifier));

    ApproxType approx = gmesh->GetApproxType();

    switch (el_type) {
        case 1:
        {   // Line
            if (plasticmaterial){
                switch (approx){
                case ApproxType::EIsoparametric:
                    gel = new ElementWithMem<ShapeOneDLin,ShapeOneDLin>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                case ApproxType::EHierarquic:
                    gel = new ElementWithMem<ShapeOneDLin,HierarquicalOneD>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                
                default:
                    PanicButton();
                    break;
                }                
            } else if (transientmaterial){
                switch (approx){
                case ApproxType::EIsoparametric:
                    gel = new ElementTransient<ShapeOneDLin,ShapeOneDLin>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                case ApproxType::EHierarquic:
                    gel = new ElementTransient<ShapeOneDLin,HierarquicalOneD>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                
                default:
                    PanicButton();
                    break;
                }                
            } else {
                switch (approx){
                case ApproxType::EIsoparametric:
                    gel = new ElementT<ShapeOneDLin,ShapeOneDLin>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                case ApproxType::EHierarquic:
                    gel = new ElementT<ShapeOneDLin,HierarquicalOneD>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                
                default:
                    PanicButton();
                    break;
                }          
            }
            gel->PrintType() = 3;
        }
            break;
        case 2:
        {
            // Linear Triangle
            if (plasticmaterial){
                switch (approx){
                case ApproxType::EIsoparametric:
                    gel = new ElementWithMem<ShapeTriangleLin,ShapeTriangleLin>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                case ApproxType::EHierarquic:
                    gel = new ElementWithMem<ShapeTriangleLin,HierarquicalTriangle>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                
                default:
                    PanicButton();
                    break;
                }             
                
            } else if (transientmaterial){
                switch (approx){
                case ApproxType::EIsoparametric:
                    gel = new ElementTransient<ShapeTriangleLin,ShapeTriangleLin>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                case ApproxType::EHierarquic:
                    gel = new ElementTransient<ShapeTriangleLin,HierarquicalTriangle>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                
                default:
                    PanicButton();
                    break;
                }       
            } else {
                switch (approx){
                case ApproxType::EIsoparametric:
                    gel = new ElementT<ShapeTriangleLin,ShapeTriangleLin>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                case ApproxType::EHierarquic:
                    gel = new ElementT<ShapeTriangleLin,HierarquicalTriangle>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                
                default:
                    PanicButton();
                    break;
                }          
            }
            gel->PrintType() = 5;           
            break;
        }
        case 3:
        {
            // Linear Quadrilateral
            if (plasticmaterial){
                switch (approx){
                case ApproxType::EIsoparametric:
                    gel = new ElementWithMem<ShapeQuadrilateralLin,ShapeQuadrilateralLin>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                case ApproxType::EHierarquic:
                    gel = new ElementWithMem<ShapeQuadrilateralLin,HierarquicalQuad>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                
                default:
                    PanicButton();
                    break;
                }                
            } else if (transientmaterial){
                switch (approx){
                case ApproxType::EIsoparametric:
                    gel = new ElementTransient<ShapeQuadrilateralLin,ShapeQuadrilateralLin>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                case ApproxType::EHierarquic:
                    gel = new ElementTransient<ShapeQuadrilateralLin,HierarquicalQuad>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                
                default:
                    PanicButton();
                    break;
                }                
            } else {
                switch (approx){
                case ApproxType::EIsoparametric:
                    gel = new ElementT<ShapeQuadrilateralLin,ShapeQuadrilateralLin>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                case ApproxType::EHierarquic:
                    gel = new ElementT<ShapeQuadrilateralLin,HierarquicalQuad>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                
                default:
                    PanicButton();
                    break;
                }          
            }
            gel->PrintType() = 9;
            break;
        }
        case 4:
        {
            // Linear Tetrahedron
            if (plasticmaterial){
                gel = new ElementWithMem<ShapeTetrahedronLin,ShapeTetrahedronLin>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
            } else if (transientmaterial){    
                gel = new ElementTransient<ShapeTetrahedronLin,ShapeTetrahedronLin>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
            } else {
                gel = new ElementT<ShapeTetrahedronLin,ShapeTetrahedronLin>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
            }
            gel->PrintType() = 10;
            break;
        }
        case 5:
        {
            // Linear Hexahedron
            if (plasticmaterial){
                gel = new ElementWithMem<ShapeHexahedron,ShapeHexahedron>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
            } else if (transientmaterial){
                gel = new ElementTransient<ShapeHexahedron,ShapeHexahedron>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
            } else {
                gel = new ElementT<ShapeHexahedron,ShapeHexahedron>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
            }
            gel->PrintType() = 12;
            break;
        }
        case 8:
        {
            // Quadratic Line
            if (plasticmaterial){
                switch (approx){
                case ApproxType::EIsoparametric:
                    gel = new ElementWithMem<ShapeOneDQua,ShapeOneDQua>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                case ApproxType::EHierarquic:
                    gel = new ElementWithMem<ShapeOneDQua,HierarquicalOneD>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                
                default:
                    PanicButton();
                    break;
                }                
            } else if (transientmaterial){
                switch (approx){
                case ApproxType::EIsoparametric:
                    gel = new ElementTransient<ShapeOneDQua,ShapeOneDQua>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                case ApproxType::EHierarquic:
                    gel = new ElementTransient<ShapeOneDQua,HierarquicalOneD>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                
                default:
                    PanicButton();
                    break;
                }                
            } else {
                switch (approx){
                case ApproxType::EIsoparametric:
                    gel = new ElementT<ShapeOneDQua,ShapeOneDQua>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                case ApproxType::EHierarquic:
                    gel = new ElementT<ShapeOneDQua,HierarquicalOneD>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                
                default:
                    PanicButton();
                    break;
                }          
            }
            gel->PrintType() = 21;           
            break;
        }
        case 9:
        {
            // Quadratic Triangle
            if (plasticmaterial){
                switch (approx){
                case ApproxType::EIsoparametric:
                    gel = new ElementWithMem<ShapeTriangleQua,ShapeTriangleQua>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                case ApproxType::EHierarquic:
                    gel = new ElementWithMem<ShapeTriangleQua,HierarquicalTriangle>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                
                default:
                    PanicButton();
                    break;
                }             
                
            } else if (transientmaterial){
                switch (approx){
                case ApproxType::EIsoparametric:
                    gel = new ElementTransient<ShapeTriangleQua,ShapeTriangleQua>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                case ApproxType::EHierarquic:
                    gel = new ElementTransient<ShapeTriangleQua,HierarquicalTriangle>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                
                default:
                    PanicButton();
                    break;
                }       
            } else {
                switch (approx){
                case ApproxType::EIsoparametric:
                    gel = new ElementT<ShapeTriangleQua,ShapeTriangleQua>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                case ApproxType::EHierarquic:
                    gel = new ElementT<ShapeTriangleQua,HierarquicalTriangle>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                
                default:
                    PanicButton();
                    break;
                }          
            }
            gel->PrintType() = 22;           
            break;
        }
        case 10:
        {
            // Quadratic Quad
            if (plasticmaterial){
                switch (approx){
                case ApproxType::EIsoparametric:
                    gel = new ElementWithMem<ShapeQuadrilateralQua,ShapeQuadrilateralQua>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                case ApproxType::EHierarquic:
                    gel = new ElementWithMem<ShapeQuadrilateralQua,HierarquicalQuad>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                
                default:
                    PanicButton();
                    break;
                }                
            } else if (transientmaterial){
                switch (approx){
                case ApproxType::EIsoparametric:
                    gel = new ElementTransient<ShapeQuadrilateralQua,ShapeQuadrilateralQua>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                case ApproxType::EHierarquic:
                    gel = new ElementTransient<ShapeQuadrilateralQua,HierarquicalQuad>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                
                default:
                    PanicButton();
                    break;
                }                
            } else {
                switch (approx){
                case ApproxType::EIsoparametric:
                    gel = new ElementT<ShapeQuadrilateralQua,ShapeQuadrilateralQua>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                case ApproxType::EHierarquic:
                    gel = new ElementT<ShapeQuadrilateralQua,HierarquicalQuad>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                
                default:
                    PanicButton();
                    break;
                }          
            }
            gel->PrintType() = 28;           
            break;
        }
        case 15:
        {
            // Point
            if (plasticmaterial){
                gel = new ElementWithMem<ShapePoint,ShapePoint>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
            } else if (transientmaterial){
                gel = new ElementTransient<ShapePoint,ShapePoint>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
            } else {
                gel = new ElementT<ShapePoint,ShapePoint>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
            }
            gel->PrintType() = 1;           
            break;
        }
        case 21:
        {
            // Cubic Triangle
            if (plasticmaterial){
                switch (approx){
                case ApproxType::EIsoparametric:
                    gel = new ElementWithMem<ShapeTriangleCub,ShapeTriangleCub>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                case ApproxType::EHierarquic:
                    gel = new ElementWithMem<ShapeTriangleCub,HierarquicalTriangle>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                
                default:
                    PanicButton();
                    break;
                }             
                
            } else if (transientmaterial){
                switch (approx){
                case ApproxType::EIsoparametric:
                    gel = new ElementTransient<ShapeTriangleCub,ShapeTriangleCub>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                case ApproxType::EHierarquic:
                    gel = new ElementTransient<ShapeTriangleCub,HierarquicalTriangle>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                
                default:
                    PanicButton();
                    break;
                }       
            } else {
                switch (approx){
                case ApproxType::EIsoparametric:
                    gel = new ElementT<ShapeTriangleCub,ShapeTriangleCub>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                case ApproxType::EHierarquic:
                    gel = new ElementT<ShapeTriangleCub,HierarquicalTriangle>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                
                default:
                    PanicButton();
                    break;
                }          
            }
            gel->PrintType() = 69;           
            break;
        }
        case 26:
        {
            // Cubic Line
            if (plasticmaterial){
                switch (approx){
                case ApproxType::EIsoparametric:
                    gel = new ElementWithMem<ShapeOneDCub,ShapeOneDCub>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                case ApproxType::EHierarquic:
                    gel = new ElementWithMem<ShapeOneDCub,HierarquicalOneD>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                
                default:
                    PanicButton();
                    break;
                }                
            } else if (transientmaterial){
                switch (approx){
                case ApproxType::EIsoparametric:
                    gel = new ElementTransient<ShapeOneDCub,ShapeOneDCub>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                case ApproxType::EHierarquic:
                    gel = new ElementTransient<ShapeOneDCub,HierarquicalOneD>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                
                default:
                    PanicButton();
                    break;
                }                
            } else {
                switch (approx){
                case ApproxType::EIsoparametric:
                    gel = new ElementT<ShapeOneDCub,ShapeOneDCub>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                case ApproxType::EHierarquic:
                    gel = new ElementT<ShapeOneDCub,HierarquicalOneD>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
                    break;
                
                default:
                    PanicButton();
                    break;
                }          
            }
            gel->PrintType() = 35;           
            break;
        }

        case 36:
        {
            PanicButton();
            // Cubic Quadrilateral
            // if (plasticmaterial){
            //     gel = new ElementWithMem<ShapeQuadrilateralCub>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
            // } else {
            //     gel = new ElementT<ShapeQuadrilateralCub>(el_identifier,Topology,gmesh,gmesh->Material(physical_identifier));
            // }
            // gel->PrintType() = 70;
            break;
        }
        
        default:
        {
            std::cout << "Element not implemented." << std::endl;
            PanicButton();
        }
            break;
    }
    return gel;
};
*/