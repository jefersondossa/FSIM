//------------------------------------------------------------------------------
//--------------------------Universidade de Sao Paulo---------------------------
//----------------------Escola de Engenharia de Sao Carlos----------------------
//----------------Departamento de Engenharia de Estruturas - SET----------------
//------------------------------Sao Carlos - 2017-------------------------------
//------------------------------------------------------------------------------
 
///-----------------------------------------------------------------------------
///----Software developed for analysis of Fluid-Structure Interaction (FSI)-----
///---problems with shell structures and incompressible fluid. The fluid is-----
///---modeled in Arlequin framework and Arbitrary Lagrangian-Eulerian (ALE)-----
///---description in mixed formulation, with triangular/tetrahedral elements----
///with quadratic (six/ten nodes) and linear (four/three nodes) approximations--
///--for velocity and pressure fields, respectively. The structure is modeled---
///-with shell elements free rotations as degrees of freedom in the nonlinear---
///-Positional Finite Element Method framework. Triangular elements with cubic--
///approximation (ten nodes) and seven degrees of freedom per node (three nodal-
///--positions, three components of unconstrained vector and strain rate. The---
///-fluid-shell coupling is performed by a partitioned loosely coupled scheme.--
///-----------------------------------------------------------------------------
  
//------------------------------------------------------------------------------
//---------------------------------Developed by---------------------------------
//-------Jeferson Wilian Dossa Fernandes and Rodolfo Andre Kuche Sanches--------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------STARTS MAIN PROGRAM-----------------------------
//------------------------------------------------------------------------------
static char help[] = "Solves the Incompressible flow problem";

// C++ standard libraries
#include <fstream> 
  
// Developed Header Files
#include "LinearAnalysis.h"
#include "NonLinearAnalysis.h"
#include "IncrementalAnalysis.h"
#include "TransientAnalysis.h"
#include "Arlequin.h"
#include "GmshTools.h"
#include "VTUGenerator.h"
#include "XDMFGenerator.h"
#include "Poisson.h"
#include "L2Projection.h"
#include "Elasticity2D.h"
#include "ElasticityPositional2D.h"
#include "GraphMesh.h"
#include "ElasticTruss.h"
#include "PositionalTruss.h"
#include "MohrCoulomb.h"
#include "VonMises.h"
#include "Elasticity3D.h"
#include "TransientElasticity2D.h"
#include "TransientPoisson.h"
#include "TransientElasticTruss.h"
#include "TransientPositionalTruss.h"
#include "TransientPositional2D.h"
#include "NavierStokes.h"
#include "TransientNavierStokes.h"

int main(int argc, char **args) {

#ifdef HAS_PETSC
    // Starts main program invoking PETSc
    PetscInitialize(&argc, &args, (char*)0, help);
#endif

    // #include "examples/turbine.h"
    // #include "examples/cylinder.h"
    // #include "examples/rectangle.h"
    // #include "examples/rectangledyn.h"
    // #include "examples/poisson.h"
    // #include "examples/arlequin.h"
    // #include "examples/cube.h"
    // #include "examples/truss.h"
    // #include "examples/mohrcoulombtest.h"
    #include "examples/notched.h"
    // #include "examples/geogrelha.h"
    // #include "examples/plasticitytest.h"
     
#ifdef HAS_PETSC
    //Finalize main program   
    PetscFinalize();
#endif

    return 0; 
}
 
 
  




 
