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
#include "src/FSInteraction.hpp"  
#include "src/fluidDomain.h"


int main(int argc, char **args) {

    // Starts main program invoking PETSc
    PetscInitialize(&argc, &args, (char*)0, help);

    int rank, size;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
 
    if (rank == 0){
        std::cout << "2D Incompressible Flows Numerical Analysis" << std::endl;
        std::cout << "Starting.." << std::endl;
        std::cout << "Type the input file name:" << std::endl;
    };

    // Defines the problem dimension
    const int dimension = 2;

    //Type definition
    typedef Fluid<dimension>         FluidModel;
    typedef Arlequin<dimension>      Arlequin;
    typedef FSInteraction<dimension> FSI;
 
//  Create problem variables 
    FluidModel coarseModel, fineModel, control;  
    Arlequin   arlequinProblem; 
    FSI        coupledProblem;  
    
//==========================================================================
//==============================PROBLEM MESHES==============================
//==========================================================================
    //COARSE 
    Geometry* fluid1 = new Geometry(0);

    Point* p000 = fluid1 -> addPoint({0.0,0.0},0.5,false);
    Point* p001 = fluid1 -> addPoint({19.5,0.0},0.5,false);
    Point* p002 = fluid1 -> addPoint({19.5,12.0},0.5,false);
    Point* p003 = fluid1 -> addPoint({0.0,12.0},0.5,false); 
     
  

    Line* l000 = fluid1 -> addLine({p000,p001});
    Line* l001 = fluid1 -> addLine({p001,p002});
    Line* l002 = fluid1 -> addLine({p002,p003});
    Line* l003 = fluid1 -> addLine({p003,p000}); 
     
    LineLoop* ll001 = fluid1 -> addLineLoop({ l000, l001, l002, l003 });
    
    double prog1 = 1.07;
    double prog2 = 1.1;

    fluid1 -> addBoundaryCondition("NEUMANN", l001, {}, {}, "GLOBAL");
    
    fluid1 -> addBoundaryCondition("DIRICHLET", l000, {}, {0.0}, "GLOBAL");
    fluid1 -> addBoundaryCondition("DIRICHLET", l002, {}, {0.0}, "GLOBAL");
    fluid1 -> addBoundaryCondition("DIRICHLET", l003, {31.5}, {0.0}, "GLOBAL");
    
    PlaneSurface* s001 = fluid1 -> addPlaneSurface({ll001});


    //FINE
    Geometry* fluid2 = new Geometry(0);
    double d1 = 0.1;
    double d2 = 0.24;
    double gz = 1.0;

    Point* p010 = fluid2 -> addPoint({4.5,5.5},d1,false);
    Point* p011 = fluid2 -> addPoint({5.5,5.5},d1,false);
    Point* p015 = fluid2 -> addPoint({5.5,5.99999},d1,false);
    Point* p016 = fluid2 -> addPoint({9.5,6.0},d1,false);
    Point* p014 = fluid2 -> addPoint({5.5,6.00001},d1,false);
    Point* p013 = fluid2 -> addPoint({5.5,6.5},d1,false);
    Point* p012 = fluid2 -> addPoint({4.5,6.5},d1,false);
    Point* p017 = fluid2 -> addPoint({3.0,4.5},d2,false);
    Point* p018 = fluid2 -> addPoint({10.5,4.5},d2,false);
    Point* p019 = fluid2 -> addPoint({3.0,7.5},d2,false);
    Point* p020 = fluid2 -> addPoint({10.5,7.5},d2,false);
    Point* p021 = fluid2 -> addPoint({3.0,4.5-gz},d2,false);
    Point* p022 = fluid2 -> addPoint({10.5,4.5-gz},d2,false);
    Point* p023 = fluid2 -> addPoint({10.5+gz,4.5},d2,false);
    Point* p024 = fluid2 -> addPoint({10.5+gz,7.5},d2,false);
    Point* p025 = fluid2 -> addPoint({10.5,7.5+gz},d2,false);
    Point* p026 = fluid2 -> addPoint({3.0,7.5+gz},d2,false);
    Point* p027 = fluid2 -> addPoint({3.0-gz,7.5},d2,false);
    Point* p028 = fluid2 -> addPoint({3.0-gz,4.5},d2,false);

    Line* l010 = fluid2 -> addLine({p010,p011});
    Line* l011 = fluid2 -> addLine({p011,p015});
    Line* l012 = fluid2 -> addLine({p015,p016});
    Line* l013 = fluid2 -> addLine({p016,p014}); 
    Line* l014 = fluid2 -> addLine({p014,p013}); 
    Line* l015 = fluid2 -> addLine({p013,p012}); 
    Line* l016 = fluid2 -> addLine({p012,p010}); 
    Line* l017 = fluid2 -> addLine({p017,p018}); 
    Line* l018 = fluid2 -> addLine({p018,p020}); 
    Line* l019 = fluid2 -> addLine({p020,p019}); 
    Line* l020 = fluid2 -> addLine({p019,p017}); 
    Line* l021 = fluid2 -> addLine({p021,p022}); 
    Line* l022 = fluid2 -> addLine({p022,p018,p023}); 
    Line* l023 = fluid2 -> addLine({p023,p024}); 
    Line* l024 = fluid2 -> addLine({p024,p020,p025}); 
    Line* l025 = fluid2 -> addLine({p025,p026}); 
    Line* l026 = fluid2 -> addLine({p026,p019,p027}); 
    Line* l027 = fluid2 -> addLine({p027,p028}); 
    Line* l028 = fluid2 -> addLine({p028,p017,p021}); 
   
    LineLoop* ll011 = fluid2 -> addLineLoop({ l017, l018, l019, l020, 
                                              l016 -> operator-(), l015 -> operator-(),
                                              l014 -> operator-(), l013 -> operator-(),
                                              l012 -> operator-(), l011 -> operator-(),
                                              l010 -> operator-() });
    LineLoop* ll012 = fluid2 -> addLineLoop({ l021, l022, l023, l024, l025, l026, l027, l028,
                                              l020 -> operator-(), l019 -> operator-(),
                                              l018 -> operator-(), l017 -> operator-() });

    fluid2 -> addBoundaryCondition("DIRICHLET", l010, {0.0}, {0.0}, "GLOBAL");
    fluid2 -> addBoundaryCondition("DIRICHLET", l011, {0.0}, {0.0}, "GLOBAL");
    fluid2 -> addBoundaryCondition("DIRICHLET", l014, {0.0}, {0.0}, "GLOBAL");
    fluid2 -> addBoundaryCondition("DIRICHLET", l015, {0.0}, {0.0}, "GLOBAL");
    fluid2 -> addBoundaryCondition("DIRICHLET", l016, {0.0}, {0.0}, "GLOBAL");

    //FSINTERFACE
    fluid2 -> addBoundaryCondition("FSINTERFACE", l012, {0.0}, {0.0}, "GLOBAL");
    fluid2 -> addBoundaryCondition("FSINTERFACE", l013, {0.0}, {0.0}, "GLOBAL");

    fluid2 -> addBoundaryCondition("NEUMANN", l017, {}, {}, "GLOBAL");
    fluid2 -> addBoundaryCondition("NEUMANN", l018, {}, {}, "GLOBAL");
    fluid2 -> addBoundaryCondition("NEUMANN", l019, {}, {}, "GLOBAL");
    fluid2 -> addBoundaryCondition("NEUMANN", l020, {}, {}, "GLOBAL");

    fluid2 -> addBoundaryCondition("GLUE", l021, {0}, {0}, "GLOBAL");
    fluid2 -> addBoundaryCondition("GLUE", l022, {0}, {0}, "GLOBAL");
    fluid2 -> addBoundaryCondition("GLUE", l023, {0}, {0}, "GLOBAL");
    fluid2 -> addBoundaryCondition("GLUE", l024, {0}, {0}, "GLOBAL");
    fluid2 -> addBoundaryCondition("GLUE", l025, {0}, {0}, "GLOBAL");
    fluid2 -> addBoundaryCondition("GLUE", l026, {0}, {0}, "GLOBAL");
    fluid2 -> addBoundaryCondition("GLUE", l027, {0}, {0}, "GLOBAL");
    fluid2 -> addBoundaryCondition("GLUE", l028, {0}, {0}, "GLOBAL");

    PlaneSurface* s002 = fluid2 -> addPlaneSurface({ll011});
    PlaneSurface* s003 = fluid2 -> addPlaneSurface({ll012});

    if (rank == 0){
  
        FluidDomain* problem = new FluidDomain(fluid1);
        problem -> addSurfaceMaterial({ s001 }, 1.0, 1.0, 1.0, "PLANE_STRESS");
        problem -> generateMesh("T6", "DELAUNAY", "coarse", "", false, false);

        FluidDomain* problem2 = new FluidDomain(fluid2);
        problem2 -> addSurfaceMaterial({ s002, s003 }, 1.0, 1.0, 1.0, "PLANE_STRESS");
        problem2 -> generateMesh("T6", "DELAUNAY", "fine", "", false, false);

        //problem -> readInput("exemplo.msh",0);
	};

//==========================================================================
//==============================PROBLEM MESHES==============================
//==========================================================================

	MPI_Barrier(PETSC_COMM_WORLD);   

    coarseModel.meshReading(fluid1,"flutter.txt","coarse.msh","mirror_coarse.txt",0);
    fineModel.meshReading(fluid2,"flutter.txt","fine.msh","mirror_fine.txt",0);

    arlequinProblem.setFluidModels(coarseModel, fineModel); 
    // arlequinProblem.solveArlequinProblem(3, 1.e-7, 2, 1); 

    char in_solid[32] = "cantilever.txt";
    coupledProblem.setArlequinAndSolidModels(arlequinProblem,in_solid);
    coupledProblem.solveFSIProblemGaussSeidelArlequin(20000);

	// coarseModel.solveTransientProblem(3,1.e-7,true);    
    // 
    // coupledProblem.setFluidAndSolidModels(coarseModel,in_solid);
    // coupledProblem.solveFSIProblemGaussSeidel(10000);

       
     // 
     // 
     // arlequinProblem.solveArlequinProblemMoving(4, 1.e-7, 2, 1); 
          
     // 
    //Finalize main program   
    PetscFinalize();
 
    return 0; 
}
 
 
  




 