
// #include <iostream>
// #include <math.h>

// #include "DenseEigen.h"
// #include "GeometricNode.h"
// #include "GeometricLine.h"
// #include "GeometricElementT.h"
// #include "GeometricElement.h"

// #define CATCH_CONFIG_MAIN
// #include <catch2/catch.hpp>
// using namespace Catch::literals;

// using std::cout;
// using std::endl;
// using std::cin;

// void CreateNode()
// {
//     //Test default constructor
//     GeometricNode *node1 = new GeometricNode(); 
//     REQUIRE(node1);
//     // node1->Print(std::cout);

//     //Test real constructor
//     VecDouble coord(3);
//     coord[0] = 1.;
//     coord[1] = 2.;
//     coord[2] = 3.;
//     int64_t index = 1000;

//     GeometricNode *node2 = new GeometricNode(coord,index);
//     VecDouble n2coord = node2->GetCoordinates();
//     bool vCoord = false;
//     if (n2coord == coord) vCoord = true;
    
//     // node2->Print(std::cout);

//     REQUIRE(node2);
//     REQUIRE(vCoord);
// }

// template<class T>
// void CreateElement(){

//     //Test default constructor
//     GeometricElementT<T> *el1 = new GeometricElementT<T>(); 
//     REQUIRE(el1);

//     //Test real constructor
//     GeometricMesh * gmesh = new GeometricMesh();
//     gmesh->SetDimension(1);



//     int iEl = 0;
//     int matid = 1;
//     VecInt nodes(2);
//     nodes[0] = iEl;
//     nodes[1] = iEl + 1;
//     auto *gel = new GeometricElementT<T>(nodes, matid, gmesh, iEl);
//     gmesh->SetElement(iEl, gel);
//     REQUIRE(gel);

//     VecDouble xi(1), x(1);
//     MatDouble NodeCo(2,3);
//     xi[0] = 0.;
//     x.setZero();
//     NodeCo.setZero();
//     // NodeCo(0,0) = -1.;
//     // NodeCo(1,0) = -1.;

//     // gmesh->Element(0)->X(xi,x);
//     // gel->X(xi,x);
    
// }

// template<int Dim>
// void CreateMesh(){

//     GeometricMesh * gmesh = new GeometricMesh(Dim);

//     int nElements = 3;
//     int nNodes = pow(nElements+1,Dim);

//     gmesh->SetNumberOfNodes(nNodes);
//     gmesh->SetNumberOfElements(nElements);
//     gmesh->SetGeometricMap(GHierarchical);

//     VecDouble coord(3);
//     coord.setZero();

//     double deltaX = 1./nElements;
//     double deltaY = 0.1;

//     int64_t index = 0;
//     GeometricNode *nodeAux = new GeometricNode(coord,index);  
//     gmesh->SetNode(nodeAux->GetIndex(),*nodeAux);

//     int matid = 1;

//     for (int iEl = 0; iEl < nElements; iEl++)
//     {
//         coord[0] = deltaX*(iEl+1);
//         coord[1] = deltaY*(iEl+1);
        
//         GeometricNode *nodeAux = new GeometricNode(coord,++index);  
//         gmesh->SetNode(nodeAux->GetIndex(),*nodeAux);

//         VecInt nodes(2);
//         nodes[0] = iEl;
//         nodes[1] = iEl + 1;
//         GeometricElement *gel = new GeometricElementT<GeometricLine>(nodes, matid, gmesh, iEl);
//     }
//     gmesh->Print(std::cout);

// }


// TEST_CASE("Geom_Test")
// {
    
//     SECTION("Node Creation")
//     {
//         CreateNode();   
//     }

//     SECTION("Element Creation")
//     {
//         CreateElement<GeometricLine>();   
//     }

//     SECTION("Mesh Creation")
//     {
//         CreateMesh<1>();   
//     }
    
// }
