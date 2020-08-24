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
        //Coarse
    //Coarse
        Geometry* fluid1 = new Geometry(0);
        double e = 5;

        double a1 = 5;
        double z = 0.;
        double a = 50.;
        double b = 150.;
        double as = 6.;

        Point* p000 = fluid1 -> addPoint({-a,-a},a1,false);
        Point* p001 = fluid1 -> addPoint({-4-e,-a},a1,false);
        Point* p002 = fluid1 -> addPoint({4+e,-a},a1,false);
        Point* p003 = fluid1 -> addPoint({b, -a},a1,false);
        Point* p004 = fluid1 -> addPoint({  -a,-4-e},a1,false);
        Point* p005 = fluid1 -> addPoint({-4-e,-4-e},a1/as,false);
        Point* p006 = fluid1 -> addPoint({ 4+e,-4-e},a1/as,false);
        Point* p007 = fluid1 -> addPoint({   b,-4-e},a1,false);
        Point* p008 = fluid1 -> addPoint({  -a, 4+e},a1,false);
        Point* p009 = fluid1 -> addPoint({-4-e, 4+e},a1/as,false);
        Point* p010 = fluid1 -> addPoint({ 4+e, 4+e},a1/as,false);
        Point* p011 = fluid1 -> addPoint({   b, 4+e},a1,false);
        Point* p012 = fluid1 -> addPoint({  -a,a},a1,false);
        Point* p013 = fluid1 -> addPoint({-4-e,a},a1,false);
        Point* p014 = fluid1 -> addPoint({ 4+e,a},a1,false);
        Point* p015 = fluid1 -> addPoint({   b,a},a1,false);
        

        Line* l001 = fluid1 -> addLine({p000,p001});
        Line* l002 = fluid1 -> addLine({p001,p002});
        Line* l003 = fluid1 -> addLine({p002,p003});
        Line* l004 = fluid1 -> addLine({p000,p004});
        Line* l005 = fluid1 -> addLine({p001,p005});        
        Line* l006 = fluid1 -> addLine({p002,p006});
        Line* l007 = fluid1 -> addLine({p003,p007});
        Line* l008 = fluid1 -> addLine({p004,p005});
        Line* l009 = fluid1 -> addLine({p005,p006});
        Line* l010 = fluid1 -> addLine({p006,p007});
        Line* l011 = fluid1 -> addLine({p004,p008});
        Line* l012 = fluid1 -> addLine({p005,p009});
        Line* l013 = fluid1 -> addLine({p006,p010});
        Line* l014 = fluid1 -> addLine({p007,p011});
        Line* l015 = fluid1 -> addLine({p008,p009});
        Line* l016 = fluid1 -> addLine({p009,p010});
        Line* l017 = fluid1 -> addLine({p010,p011});        
        Line* l018 = fluid1 -> addLine({p008,p012});
        Line* l019 = fluid1 -> addLine({p009,p013});
        Line* l020 = fluid1 -> addLine({p010,p014});
        Line* l021 = fluid1 -> addLine({p011,p015});
        Line* l022 = fluid1 -> addLine({p012,p013});
        Line* l023 = fluid1 -> addLine({p013,p014});
        Line* l024 = fluid1 -> addLine({p014,p015});
       
        LineLoop* ll001 = fluid1 -> addLineLoop({ l001, l005, l008 -> operator-(), l004 -> operator-() });
        LineLoop* ll002 = fluid1 -> addLineLoop({ l002, l006, l009 -> operator-(), l005 -> operator-() });
        LineLoop* ll003 = fluid1 -> addLineLoop({ l003, l007, l010 -> operator-(), l006 -> operator-() });
        LineLoop* ll004 = fluid1 -> addLineLoop({ l008, l012, l015 -> operator-(), l011 -> operator-() });
        LineLoop* ll005 = fluid1 -> addLineLoop({ l009, l013, l016 -> operator-(), l012 -> operator-() });
        LineLoop* ll006 = fluid1 -> addLineLoop({ l010, l014, l017 -> operator-(), l013 -> operator-() });
        LineLoop* ll007 = fluid1 -> addLineLoop({ l015, l019, l022 -> operator-(), l018 -> operator-() });
        LineLoop* ll008 = fluid1 -> addLineLoop({ l016, l020, l023 -> operator-(), l019 -> operator-() });
        LineLoop* ll009 = fluid1 -> addLineLoop({ l017, l021, l024 -> operator-(), l020 -> operator-() });
        
        //std::vector<LineLoop*> lin = {ll0, ll1};
        PlaneSurface* s001 = fluid1 -> addPlaneSurface({ll001});
        PlaneSurface* s002 = fluid1 -> addPlaneSurface({ll002});
        PlaneSurface* s003 = fluid1 -> addPlaneSurface({ll003});
        PlaneSurface* s004 = fluid1 -> addPlaneSurface({ll004});
        PlaneSurface* s005 = fluid1 -> addPlaneSurface({ll005});
        PlaneSurface* s006 = fluid1 -> addPlaneSurface({ll006});
        PlaneSurface* s007 = fluid1 -> addPlaneSurface({ll007});
        PlaneSurface* s008 = fluid1 -> addPlaneSurface({ll008});
        PlaneSurface* s009 = fluid1 -> addPlaneSurface({ll009});

        
        fluid1 -> addBoundaryCondition("NEUMANN", l005, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l006, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l007, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l012, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l013, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l014, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l019, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l020, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l021, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l008, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l009, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l010, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l015, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l016, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l017, {}, {}, "GLOBAL");
        
        fluid1 -> addBoundaryCondition("DIRICHLET", l001, {}, {0.0}, "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l002, {}, {0.0}, "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l003, {}, {0.0}, "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l022, {}, {0.0}, "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l023, {}, {0.0}, "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l024, {}, {0.0}, "GLOBAL");

        fluid1 -> addBoundaryCondition("DIRICHLET", l004, {20.0}, {0.0}, "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l011, {20.0}, {0.0}, "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l018, {20.0}, {0.0}, "GLOBAL");
        
        
        // Geometry* fluid1 = new Geometry(0);
        // double a1 = 1.8;
        // double z = 0.;
        // double a = 40.;
        // double b = 120.;

        // Point* p000 = fluid1 -> addPoint({-a,-a},a1,false);
        // Point* p001 = fluid1 -> addPoint({ z,-a},a1,false);
        // Point* p002 = fluid1 -> addPoint({ b,-a},a1,false);
        // Point* p003 = fluid1 -> addPoint({-a, z},a1,false);
        // Point* p004 = fluid1 -> addPoint({ z, z},a1/3,false);
        // Point* p005 = fluid1 -> addPoint({ b, z},a1,false);
        // Point* p006 = fluid1 -> addPoint({-a, a},a1,false);
        // Point* p007 = fluid1 -> addPoint({ z, a},a1,false);
        // Point* p008 = fluid1 -> addPoint({ b, a},a1,false);

        // Line* l001 = fluid1 -> addLine({p000,p001});
        // Line* l002 = fluid1 -> addLine({p001,p002});
        // Line* l003 = fluid1 -> addLine({p000,p003});
        // Line* l004 = fluid1 -> addLine({p001,p004});
        // Line* l005 = fluid1 -> addLine({p002,p005});        
        // Line* l006 = fluid1 -> addLine({p003,p004});
        // Line* l007 = fluid1 -> addLine({p004,p005});
        // Line* l008 = fluid1 -> addLine({p003,p006});
        // Line* l009 = fluid1 -> addLine({p004,p007});
        // Line* l010 = fluid1 -> addLine({p005,p008});
        // Line* l011 = fluid1 -> addLine({p006,p007});
        // Line* l012 = fluid1 -> addLine({p007,p008});
       
        // LineLoop* ll001 = fluid1 -> addLineLoop({ l001, l004, l006 -> operator-(), l003 -> operator-() });
        // LineLoop* ll002 = fluid1 -> addLineLoop({ l002, l005, l007 -> operator-(), l004 -> operator-() });
        // LineLoop* ll003 = fluid1 -> addLineLoop({ l006, l009, l011 -> operator-(), l008 -> operator-() });
        // LineLoop* ll004 = fluid1 -> addLineLoop({ l007, l010, l012 -> operator-(), l009 -> operator-() });
        
        // //std::vector<LineLoop*> lin = {ll0, ll1};
        // PlaneSurface* s001 = fluid1 -> addPlaneSurface({ll001});
        // PlaneSurface* s002 = fluid1 -> addPlaneSurface({ll002});
        // PlaneSurface* s003 = fluid1 -> addPlaneSurface({ll003});
        // PlaneSurface* s004 = fluid1 -> addPlaneSurface({ll004});
        
        // fluid1 -> addBoundaryCondition("NEUMANN", l004, {}, {}, "GLOBAL");
        // fluid1 -> addBoundaryCondition("NEUMANN", l005, {}, {}, "GLOBAL");
        // fluid1 -> addBoundaryCondition("NEUMANN", l006, {}, {}, "GLOBAL");
        // fluid1 -> addBoundaryCondition("NEUMANN", l007, {}, {}, "GLOBAL");
        // fluid1 -> addBoundaryCondition("NEUMANN", l009, {}, {}, "GLOBAL");
        // fluid1 -> addBoundaryCondition("NEUMANN", l010, {}, {}, "GLOBAL");
        
        // fluid1 -> addBoundaryCondition("DIRICHLET", l001, {}, {0.0}, "GLOBAL");
        // fluid1 -> addBoundaryCondition("DIRICHLET", l002, {}, {0.0}, "GLOBAL");
        // fluid1 -> addBoundaryCondition("DIRICHLET", l011, {}, {0.0}, "GLOBAL");
        // fluid1 -> addBoundaryCondition("DIRICHLET", l012, {}, {0.0}, "GLOBAL");

        // fluid1 -> addBoundaryCondition("DIRICHLET", l003, {350.0}, {0.0}, "GLOBAL");
        // fluid1 -> addBoundaryCondition("DIRICHLET", l008, {350.0}, {0.0}, "GLOBAL");
        
        


        // Fine
        Geometry* fluid2 = new Geometry(1);

        double a2 = 0.15;
        double a4 = 0.1;
        double a3 = 0.5;

		Point* p100 = fluid2 -> addPoint({ 0.000, 0.000},a2,false);
        Point* p101 = fluid2 -> addPoint({ 1.850, 0.000},a2,false);
        Point* p102 = fluid2 -> addPoint({-1.850, 0.000},a2,false);
        Point* p103 = fluid2 -> addPoint({-3.850, 0.000},a2,false);
        Point* p104 = fluid2 -> addPoint({ 0.150, 0.000},a4,false);
        Point* p105 = fluid2 -> addPoint({-0.150, 0.000},a4,false);
        Point* p106 = fluid2 -> addPoint({ 3.850, 0.000},a2,false);
        Point* p107 = fluid2 -> addPoint({ 0.000,-7.00},a3,false);
        Point* p108 = fluid2 -> addPoint({ 7.00, 0.000},a3,false);
        Point* p109 = fluid2 -> addPoint({ 0.000, 7.00},a3,false);
        Point* p110 = fluid2 -> addPoint({-7.00, 0.000},a3,false);
        Point* p111 = fluid2 -> addPoint({ 0.000,-9.00},a3,false);
        Point* p112 = fluid2 -> addPoint({ 9.00, 0.000},a3,false);
        Point* p113 = fluid2 -> addPoint({ 0.000, 9.00},a3,false);
        Point* p114 = fluid2 -> addPoint({-9.00, 0.000},a3,false);
        
        Line* l100 = fluid2 -> addLine({p104,p102,p103});
        Line* l101 = fluid2 -> addLine({p105,p101,p106});
        Line* l102 = fluid2 -> addLine({p107,p100,p108});
        Line* l103 = fluid2 -> addLine({p108,p100,p109});
        Line* l104 = fluid2 -> addLine({p109,p100,p110});
        Line* l105 = fluid2 -> addLine({p110,p100,p107});
        Line* l106 = fluid2 -> addLine({p111,p100,p112});
        Line* l107 = fluid2 -> addLine({p112,p100,p113});
        Line* l108 = fluid2 -> addLine({p113,p100,p114});
        Line* l109 = fluid2 -> addLine({p114,p100,p111});
        Line* l110 = fluid2 -> addLine({p104,p102,p103});
        Line* l111 = fluid2 -> addLine({p105,p101,p106});
        Line* l112 = fluid2 -> addLine({p110,p103});
        // Line* l113 = fluid2 -> addLine({p105,p104});
        Line* l114 = fluid2 -> addLine({p106,p108});
        Line* l115 = fluid2 -> addLine({p104,p100,p105});
        Line* l116 = fluid2 -> addLine({p105,p100,p104});
         
        LineLoop* ll101 = fluid2->addLineLoop({ l112 -> operator-(), l110, l115 -> operator-(), l101 -> operator-(), l114 -> operator-(), l102, l105 });
        LineLoop* ll102 = fluid2->addLineLoop({ l112, l100 -> operator-(), l116 -> operator-(), l111, l114, l103, l104 });
        LineLoop* ll103 = fluid2->addLineLoop({ l106, l107, l108, l109, l102 -> operator-(), l103 -> operator-(), l104 -> operator-(), l105 -> operator-() });
        
        //std::vector<LineLoop*> lin = {ll2, ll1};
        //std::vector<LineLoop*> lin = {ll2};
        //PlaneSurface* s101 = fluid2 -> addPlaneSurface(lin);
        PlaneSurface* s101 = fluid2 -> addPlaneSurface({ll101});
        PlaneSurface* s102 = fluid2 -> addPlaneSurface({ll102});
        PlaneSurface* s103 = fluid2 -> addPlaneSurface({ll103}); 
        
        fluid2 -> addBoundaryCondition("NEUMANN", l102, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l103, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l104, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l105, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l112, {}, {}, "GLOBAL");
        // fluid2 -> addBoundaryCondition("NEUMANN", l113, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l114, {}, {}, "GLOBAL");
        
        fluid2 -> addBoundaryCondition("GLUE", l106, {0}, {0}, "GLOBAL");
        fluid2 -> addBoundaryCondition("GLUE", l107, {0}, {0}, "GLOBAL");
        fluid2 -> addBoundaryCondition("GLUE", l108, {0}, {0}, "GLOBAL");
        fluid2 -> addBoundaryCondition("GLUE", l109, {0}, {0}, "GLOBAL");
        
        fluid2 -> addBoundaryCondition("FSINTERFACE", l100, {0}, {0}, "GLOBAL"); // FSI
        fluid2 -> addBoundaryCondition("FSINTERFACE", l101, {0}, {0}, "GLOBAL"); // FSI
        fluid2 -> addBoundaryCondition("FSINTERFACE", l110, {0}, {0}, "GLOBAL"); // FSI
        fluid2 -> addBoundaryCondition("FSINTERFACE", l111, {0}, {0}, "GLOBAL"); // FSI
        fluid2 -> addBoundaryCondition("DIRICHLET", l115, {0}, {0}, "GLOBAL");
        fluid2 -> addBoundaryCondition("DIRICHLET", l116, {0}, {0}, "GLOBAL");       

    if (rank == 0){
  
        FluidDomain* problem = new FluidDomain(fluid1);
        problem -> addSurfaceMaterial({ s001,s002,s003,s004,s005,s006,s007,s008,s009 }, 1.0, 1.0, 1.0, "PLANE_STRESS");
        problem -> generateMesh("T6", "DELAUNAY", "coarse", "", false, true);

        FluidDomain* problem2 = new FluidDomain(fluid2);
        problem2 -> addSurfaceMaterial({ s101,s102,s103 }, 1.0, 1.0, 1.0, "PLANE_STRESS");
        problem2 -> generateMesh("T6", "DELAUNAY", "fine", "", false, true);

        //problem -> readInput("exemplo.msh",0);
	};
 
//==========================================================================
//==============================PROBLEM MESHES==============================
//==========================================================================

	MPI_Barrier(PETSC_COMM_WORLD);   

    coarseModel.meshReading(fluid1,"problem_data.txt","coarse.msh","mirror.txt",0);
    fineModel.meshReading(fluid2,"problem_data.txt","fine.msh","mirror_fine.txt",0);
   // }    
	MPI_Barrier(PETSC_COMM_WORLD);
  
    char in_solid[32] = "turbine.txt";
         
  
    arlequinProblem.setFluidModels(coarseModel, fineModel);
    coupledProblem.setArlequinAndSolidModels(arlequinProblem,in_solid);
 
    // coarseModel.readInitialValues("saidaVelCoarse102000.txt","saidaVelCoarse102001.txt");
    // fineModel.readInitialValues("saidaVelFine102000.txt","saidaVelFine102001.txt");
    
    arlequinProblem.solveArlequinProblem(4, 1.e-7, 2, 1);

    coupledProblem.solveFSIProblemGaussSeidelArlequin(100000);

 
        
    // arlequinProblem.setFluidModels(coarseModel, fineModel) ;  

    // // arlequinProblem.solveArlequinProblem(3, 1.e-7, 2, 1); 
    
    // arlequinProblem.solveArlequinProblem(3, 1.e-7, 2, 1);
     
    //Finalize main program   
    PetscFinalize();
 
    return 0; 
}
 
 
  




 
