#ifndef GRAPHMESH_H
#define GRAPHMESH_H

#include<fstream>
#include<iostream>
#include<string>
#include "CompMesh.h"

class GraphMesh
{
private:
    std::vector<VecInt>    fElementConnects;
    std::vector<VecDouble> fNodes;
    std::vector<int>       fElementTypes;
    std::map<int64_t,int64_t> fGNodeToMNode;
    std::map<int64_t,int64_t> fMNodeToGNode;

public:
    GraphMesh(CompMesh *cmesh);

    ~GraphMesh(){ };

    int64_t NElements(){return fElementConnects.size();}
    int64_t NNodes(){return fNodes.size();}

    VecDouble &Node(int64_t index){return fNodes[index];}
    VecInt &Connect(int64_t index){return fElementConnects[index];}
    int &ElType(int64_t index){return fElementTypes[index];}
    int64_t GraphNodeToMeshNode(int64_t graphnode){return fGNodeToMNode[graphnode];};
    int64_t MeshNodeToGraphNode(int64_t meshnode){return fMNodeToGNode[meshnode];};
};

#endif