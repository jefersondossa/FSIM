#ifndef GRAPHMESH_H
#define GRAPHMESH_H

#include<fstream>
#include<iostream>
#include<string>
#include "CompMesh.h"

class CompMesh;

/// @brief Implements the graphical mesh, used to print results in a .vtu or .xdmf file
// It differs from the computational mesh because all elements are disconnected
class GraphMesh
{
private:
    // vector of element connects
    std::vector<VecInt>    fElementConnects;
    // vector of mesh nodes
    std::vector<VecDouble> fNodes;
    // vector of element types
    std::vector<int>       fElementTypes;
    // a map relating a graph node to the comp mesh node 
    std::map<int64_t,int64_t> fGNodeToMNode;
    // a map relating a comp mesh node to the graph node
    std::map<int64_t,int64_t> fMNodeToGNode;
    // a map relating a graph element to the comp mesh element
    std::map<int64_t,int64_t> fGElementToMElement;

public:
    // graph mesh constructor
    GraphMesh(CompMesh *cmesh);

    // Default destructor
    ~GraphMesh(){ };

    /// @brief Returns the number of elements in the graph mesh
    /// @return number of elements
    int64_t NElements(){
        return fElementConnects.size();
    }

    /// @brief Returns the number of nodes in the graph mesh 
    /// @return number of nodes
    int64_t NNodes(){
        return fNodes.size();
    }

    /// @brief Returns node coordinates
    /// @param index node index
    /// @return node coordinates
    VecDouble &Node(int64_t index){
        return fNodes[index];
    }

    /// @brief Returns element connectivity
    /// @param index element index
    /// @return element connectivity
    VecInt &Connect(int64_t index){return fElementConnects[index];}

    /// @brief Returns the element type
    /// @param index element index
    /// @return element type
    int &ElType(int64_t index){return fElementTypes[index];}

    /// @brief Returns the correspondence between a graph node to the comp mesh node
    /// @param graphnode graph node
    /// @return comp mesh node
    int64_t GraphNodeToMeshNode(int64_t graphnode){
        return fGNodeToMNode[graphnode];
    };
    
    /// @brief Returns the correspondence between a comp mesh node to the graph node
    /// @param meshnode comp mesh node
    /// @return graph node
    int64_t MeshNodeToGraphNode(int64_t meshnode){
        return fMNodeToGNode[meshnode];
    };

    int64_t GraphElementToMeshElement(int64_t graphel){
        return fGElementToMElement[graphel];
    }
};

#endif