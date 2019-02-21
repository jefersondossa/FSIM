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
   //FSI        coupledProblem;  
    
//==========================================================================
//==============================PROBLEM MESHES==============================
//==========================================================================
        //Coarse
        Geometry* fluid1 = new Geometry(0);

        double x1 = -18.5;
        double x2 = -2.5;
        double x3 = 2.5;
        double x4 = 38.5;
        double x5 = -2;
        double x6 = -0.5;
        double x7 = 0.0;
        double x8 = 0.5;
        double x9 = 2;

		double y1 = -25.5;
        double y2 = -2.5;
        double y3 = 2.5;
        double y4 = 25.5;
        double y5 = -2;
        double y6 = -0.5;
        double y7 = 0.0;
        double y8 = 0.5;
        double y9 = 2;

        double h1 = 5; double h2 = 3; double h3 = 10;
        double v1 = 5; double v2 = 3;

        Point* p001 = fluid1 -> addPoint({x1,y1},0.5,false);
        Point* p002 = fluid1 -> addPoint({x2,y1},0.5,false);
        Point* p003 = fluid1 -> addPoint({x3,y1},0.5,false);
        Point* p004 = fluid1 -> addPoint({x4,y1},0.5,false);
        Point* p005 = fluid1 -> addPoint({x1,y2},0.5,false);
        Point* p006 = fluid1 -> addPoint({x2,y2},0.5,false);
        Point* p007 = fluid1 -> addPoint({x3,y2},0.5,false);
        Point* p008 = fluid1 -> addPoint({x4,y2},0.5,false);
        Point* p009 = fluid1 -> addPoint({x1,y3},0.5,false);
        Point* p010 = fluid1 -> addPoint({x2,y3},0.5,false);
        Point* p011 = fluid1 -> addPoint({x3,y3},0.5,false);
        Point* p012 = fluid1 -> addPoint({x4,y3},0.5,false);
        Point* p013 = fluid1 -> addPoint({x1,y4},0.5,false);
        Point* p014 = fluid1 -> addPoint({x2,y4},0.5,false);
        Point* p015 = fluid1 -> addPoint({x3,y4},0.5,false);
        Point* p016 = fluid1 -> addPoint({x4,y4},0.5,false);

        Line* l001 = fluid1 -> addLine({p001,p002});
        Line* l002 = fluid1 -> addLine({p002,p003});
        Line* l003 = fluid1 -> addLine({p003,p004});
        Line* l004 = fluid1 -> addLine({p001,p005}); 
        Line* l005 = fluid1 -> addLine({p002,p006});
        Line* l006 = fluid1 -> addLine({p003,p007});
        Line* l007 = fluid1 -> addLine({p004,p008});
        Line* l008 = fluid1 -> addLine({p005,p006}); 
        Line* l009 = fluid1 -> addLine({p006,p007});
        Line* l010 = fluid1 -> addLine({p007,p008});
        Line* l011 = fluid1 -> addLine({p005,p009});
        Line* l012 = fluid1 -> addLine({p006,p010}); 
        Line* l013 = fluid1 -> addLine({p007,p011});
        Line* l014 = fluid1 -> addLine({p008,p012});
        Line* l015 = fluid1 -> addLine({p009,p010});
        Line* l016 = fluid1 -> addLine({p010,p011}); 
        Line* l017 = fluid1 -> addLine({p011,p012});
        Line* l018 = fluid1 -> addLine({p009,p013});
        Line* l019 = fluid1 -> addLine({p010,p014});
        Line* l020 = fluid1 -> addLine({p011,p015}); 
        Line* l021 = fluid1 -> addLine({p012,p016});
        Line* l022 = fluid1 -> addLine({p013,p014});
        Line* l023 = fluid1 -> addLine({p014,p015});
        Line* l024 = fluid1 -> addLine({p015,p016}); 
         
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

        double prog1 = 1.07;
        double prog2 = 1.1;

        fluid1 -> transfiniteLine({l001}, h1, 1/prog2);
        fluid1 -> transfiniteLine({l002}, h2);
        fluid1 -> transfiniteLine({l003}, h3, prog1);
        fluid1 -> transfiniteLine({l004}, v1, 1/prog1);
        fluid1 -> transfiniteLine({l005}, v1, 1/prog1);
        fluid1 -> transfiniteLine({l006}, v1, 1/prog1);
        fluid1 -> transfiniteLine({l007}, v1, 1/prog1);
        fluid1 -> transfiniteLine({l008}, h1, 1/prog2);
        fluid1 -> transfiniteLine({l009}, h2);
        fluid1 -> transfiniteLine({l010}, h3, prog1);
        fluid1 -> transfiniteLine({l011}, v2);
        fluid1 -> transfiniteLine({l012}, v2);
        fluid1 -> transfiniteLine({l013}, v2);
        fluid1 -> transfiniteLine({l014}, v2);
        fluid1 -> transfiniteLine({l015}, h1, 1/prog2);
        fluid1 -> transfiniteLine({l016}, h2);
        fluid1 -> transfiniteLine({l017}, h3, prog1);
        fluid1 -> transfiniteLine({l018}, v1, prog1);
        fluid1 -> transfiniteLine({l019}, v1, prog1);
        fluid1 -> transfiniteLine({l020}, v1, prog1);
        fluid1 -> transfiniteLine({l021}, v1, prog1);
        fluid1 -> transfiniteLine({l022}, h1, 1/prog2);
        fluid1 -> transfiniteLine({l023}, h2);
        fluid1 -> transfiniteLine({l024}, h3, prog1);

        fluid1 -> transfiniteSurface({s001}, "Alternated", {p001,p002,p006,p005});
        fluid1 -> transfiniteSurface({s002}, "Alternated", {p002,p003,p007,p006});
        fluid1 -> transfiniteSurface({s003}, "Alternated", {p003,p004,p008,p007});
        fluid1 -> transfiniteSurface({s004}, "Alternated", {p005,p006,p010,p009});
        fluid1 -> transfiniteSurface({s005}, "Alternated", {p006,p007,p011,p010});
        fluid1 -> transfiniteSurface({s006}, "Alternated", {p007,p008,p012,p011});
        fluid1 -> transfiniteSurface({s007}, "Alternated", {p009,p010,p014,p013});
        fluid1 -> transfiniteSurface({s008}, "Alternated", {p010,p011,p015,p014});
        fluid1 -> transfiniteSurface({s009}, "Alternated", {p011,p012,p016,p015});

        fluid1 -> addBoundaryCondition("NEUMANN", l005, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l006, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l008, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l009, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l010, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l012, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l013, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l015, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l016, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l017, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l019, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l020, {}, {}, "GLOBAL");

        fluid1 -> addBoundaryCondition("DIRICHLET", l001, {1.0}, {0.0}, "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l002, {1.0}, {0.0}, "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l003, {1.0}, {0.0}, "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l022, {1.0}, {0.0}, "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l023, {1.0}, {0.0}, "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l024, {1.0}, {0.0}, "GLOBAL");

        fluid1 -> addBoundaryCondition("NEUMANN", l007, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l014, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l021, {}, {}, "GLOBAL");

        fluid1 -> addBoundaryCondition("DIRICHLET", l004, {1.0}, {0.0}, "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l011, {1.0}, {0.0}, "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l018, {1.0}, {0.0}, "GLOBAL");
        


        // Fine
        Geometry* fluid2 = new Geometry(1);

        Point* p101 = fluid2 -> addPoint({x7,y7},0.1,false);
        Point* p102 = fluid2 -> addPoint({x7,y6},0.1,false);
        Point* p103 = fluid2 -> addPoint({x8,y7},0.1,false);
        Point* p104 = fluid2 -> addPoint({x7,y8},0.1,false);
        Point* p105 = fluid2 -> addPoint({x6,y7},0.1,false);
        Point* p106 = fluid2 -> addPoint({x7,y5},0.1,false);
        Point* p107 = fluid2 -> addPoint({x9,y7},0.1,false);
        Point* p108 = fluid2 -> addPoint({x7,y9},0.1,false);
        Point* p109 = fluid2 -> addPoint({x5,y7},0.1,false);
        Point* p110 = fluid2 -> addPoint({x7,y2},0.1,false);
        Point* p111 = fluid2 -> addPoint({x3,y7},0.1,false);
        Point* p112 = fluid2 -> addPoint({x7,y3},0.1,false);
        Point* p113 = fluid2 -> addPoint({x2,y7},0.1,false);
        
        Line* l101 = fluid2 -> addLine({p102,p001,p103});
        Line* l102 = fluid2 -> addLine({p103,p001,p104});
        Line* l103 = fluid2 -> addLine({p104,p001,p105});
        Line* l104 = fluid2 -> addLine({p105,p001,p102});
        Line* l105 = fluid2 -> addLine({p106,p001,p107});
        Line* l106 = fluid2 -> addLine({p107,p001,p108});
        Line* l107 = fluid2 -> addLine({p108,p001,p109});
        Line* l108 = fluid2 -> addLine({p109,p001,p106});
        Line* l109 = fluid2 -> addLine({p110,p001,p111});
        Line* l110 = fluid2 -> addLine({p111,p001,p112});
        Line* l111 = fluid2 -> addLine({p112,p001,p113});
        Line* l112 = fluid2 -> addLine({p113,p001,p110});
        Line* l113 = fluid2 -> addLine({p106,p002});
        Line* l114 = fluid2 -> addLine({p103,p107});
        Line* l115 = fluid2 -> addLine({p104,p108});
        Line* l116 = fluid2 -> addLine({p109,p105});
        Line* l117 = fluid2 -> addLine({p110,p106});
        Line* l118 = fluid2 -> addLine({p107,p111});
        Line* l119 = fluid2 -> addLine({p108,p112});
        Line* l120 = fluid2 -> addLine({p113,p109});
        

        LineLoop* ll101 = fluid2->addLineLoop({ l105, l114 -> operator-(), l101 -> operator-(), l113 -> operator-() });
        LineLoop* ll102 = fluid2->addLineLoop({ l114, l106, l115 -> operator-(), l102 -> operator-() });
        LineLoop* ll103 = fluid2->addLineLoop({ l107, l116, l103 -> operator-(), l115 });
        LineLoop* ll104 = fluid2->addLineLoop({ l108, l113, l104 -> operator-(), l116 -> operator-() });
        LineLoop* ll105 = fluid2->addLineLoop({ l109, l118 -> operator-(), l105 -> operator-(), l117 -> operator-() });
        LineLoop* ll106 = fluid2->addLineLoop({ l118, l110, l119 -> operator-(), l106 -> operator-() });
        LineLoop* ll107 = fluid2->addLineLoop({ l111, l120, l107 -> operator-(), l119 });
        LineLoop* ll108 = fluid2->addLineLoop({ l112, l117, l108 -> operator-(), l120 -> operator-() });

        //std::vector<LineLoop*> lin = {ll2, ll1};
        //std::vector<LineLoop*> lin = {ll2};
        //PlaneSurface* s101 = fluid2 -> addPlaneSurface(lin);
        PlaneSurface* s101 = fluid2 -> addPlaneSurface({ll101});
        PlaneSurface* s102 = fluid2 -> addPlaneSurface({ll102});
        PlaneSurface* s103 = fluid2 -> addPlaneSurface({ll103});
        PlaneSurface* s104 = fluid2 -> addPlaneSurface({ll104});
        PlaneSurface* s105 = fluid2 -> addPlaneSurface({ll105});
        PlaneSurface* s106 = fluid2 -> addPlaneSurface({ll106});
        PlaneSurface* s107 = fluid2 -> addPlaneSurface({ll107});
        PlaneSurface* s108 = fluid2 -> addPlaneSurface({ll108});
       

        fluid2 -> addBoundaryCondition("NEUMANN", l105, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l106, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l107, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l108, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l113, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l114, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l115, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l116, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l117, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l118, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l119, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l120, {}, {}, "GLOBAL");

        fluid2 -> addBoundaryCondition("GLUE", l109, {0}, {0}, "GLOBAL");
        fluid2 -> addBoundaryCondition("GLUE", l110, {0}, {0}, "GLOBAL");
        fluid2 -> addBoundaryCondition("GLUE", l111, {0}, {0}, "GLOBAL");
        fluid2 -> addBoundaryCondition("GLUE", l112, {0}, {0}, "GLOBAL");
        
        fluid2 -> addBoundaryCondition("DIRICHLET", l101, {0}, {0}, "GLOBAL");
        fluid2 -> addBoundaryCondition("DIRICHLET", l102, {0}, {0}, "GLOBAL");
        fluid2 -> addBoundaryCondition("DIRICHLET", l103, {0}, {0}, "GLOBAL");
        fluid2 -> addBoundaryCondition("DIRICHLET", l104, {0}, {0}, "GLOBAL");

        //Transfinite lines 
        int circ = 20; int r1 = 20; int r2 = 3;
        //corners
        fluid2 -> transfiniteLine({ l101 }, circ);
        fluid2 -> transfiniteLine({ l102 }, circ);
        fluid2 -> transfiniteLine({ l103 }, circ);
        fluid2 -> transfiniteLine({ l104 }, circ);
        fluid2 -> transfiniteLine({ l105 }, circ);
        fluid2 -> transfiniteLine({ l106 }, circ);
        fluid2 -> transfiniteLine({ l107 }, circ);
        fluid2 -> transfiniteLine({ l108 }, circ);
        fluid2 -> transfiniteLine({ l109 }, circ);
        fluid2 -> transfiniteLine({ l110 }, circ);
        fluid2 -> transfiniteLine({ l111 }, circ);
        fluid2 -> transfiniteLine({ l112 }, circ);

        //r1
        fluid2 -> transfiniteLine({ l113 }, r1);
        fluid2 -> transfiniteLine({ l114 }, r1);
        fluid2 -> transfiniteLine({ l115 }, r1);
        fluid2 -> transfiniteLine({ l116 }, r1);
        //2r
        fluid2 -> transfiniteLine({ l117 }, r2);
        fluid2 -> transfiniteLine({ l118 }, r2);
        fluid2 -> transfiniteLine({ l119 }, r2);
        fluid2 -> transfiniteLine({ l120 }, r2);
       
        //Surfaces
        fluid2 -> transfiniteSurface({s101}, "Alternated", {p102,p106,p107,p103});
        fluid2 -> transfiniteSurface({s102}, "Alternated", {p103,p107,p108,p104});
        fluid2 -> transfiniteSurface({s103}, "Alternated", {p104,p108,p109,p105});
        fluid2 -> transfiniteSurface({s104}, "Alternated", {p105,p109,p106,p102});
        fluid2 -> transfiniteSurface({s105}, "Alternated", {p106,p110,p111,p107});
        fluid2 -> transfiniteSurface({s106}, "Alternated", {p107,p111,p112,p108});
        fluid2 -> transfiniteSurface({s107}, "Alternated", {p108,p112,p113,p109});
        fluid2 -> transfiniteSurface({s108}, "Alternated", {p109,p113,p110,p106});
        

    if (rank == 0){
  
        FluidDomain* problem = new FluidDomain(fluid1);
        problem -> addSurfaceMaterial({ s001,s002,s003,s004,s005,s006,s007,s008,s009 }, 1.0, 1.0, 1.0, "PLANE_STRESS");
        problem -> generateMesh("T6", "DELAUNAY", "coarse", "", false, true);

        FluidDomain* problem2 = new FluidDomain(fluid2);
        problem2 -> addSurfaceMaterial({ s101,s102,s103,s104,s105,s106,s107,s108 }, 1.0, 1.0, 1.0, "PLANE_STRESS");
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


 	//control.dataReading("cylinder_coarse_ref.txt","mirror.txt");
    //control.solveTransientProblem(2, 1.e-6, 2); 

    //Data reading     
    // Cylinder
    // fineModel.dataReading("cyl_fine_str3.txt","mirror_fine.txt");  
    // coarseModel.dataReading("cyl_coarse2.txt","mirror_coarse.txt");
    // fineModel.dataReading("cylinder_fine_ref2.txt","mirror_fine.txt");  
    // coarseModel.dataReading("cylinder_coarse_ref.txt","mirror_coarse.txt");
    // fineModel.dataReading("cyl_fine_novo3.txt","mirror_fine.txt");  
    // coarseModel.dataReading("cyl_coarse_novo2.txt","mirror_coarse.txt");
  
    // Cavity   
    // fineModel.dataReading("fine4.txt","mirror_fine.txt"); 
    // coarseModel.dataReading("coarse.txt","mirror_coarse.txt");
  
    // Helice 
    // fineModel.dataReading("cavpropfcir.txt","mirror_fine.txt"); 
    // coarseModel.dataReading("cavity1.txt","mirror_coarse.txt");

    // Flutter
    // fineModel.dataReading("flutter_fine2.txt","mirror_fine.txt"); 
    // coarseModel.dataReading("flutter_coarse.txt","mirror_coarse.txt");
     
    //control.dataReading("32x32.txt","mirror_control.txt");
 
    //char in_solid[32] = "cantilever.txt";
        
  
    // arlequinProblem.setFluidModels(coarseModel, fineModel);
    // coupledProblem.setArlequinAndSolidModels(arlequinProblem,in_solid);
 
    // coupledProblem.solveFSIProblemGaussSeidelArlequin(1000);

    //control.dataReading(".txt","mirror_coarse.txt");
  
    //Problem solving  
         
    //Solve Problem function needs three parameters:  
    //1- The maximum number of iterations in the Newton-Raphson process
    //2- The maximum relative error in the Newton-Raphson process (DU)
    //3- The type of problem to be solved:  
    //   1- Steady Stokes problem;   
    //   2- Steady Navier-Stokes problem (solves the steady Stokes
    //      in the first step to perform the initial guess);  
     
    //coarseModel.solveTransientProblem(1, 1.e-15, 2);
 
    //fluidCoarseMesh.solveTransientProblem(4, 1.e-16, 2); 
  
    
    //control.solveTransientProblem(2, 1.e-6, 2); 
     
    //fineModel.solveTransientProblem(2, 1.e-16, 2);  
 
    // fineModel.solveSteadyLaplaceProblem(8, 1.e-16, 2);
         
  
    //1- The maximum number of iterations in the Newton-Raphson process
    //2- The maximum relative error in the Newton-Raphson process (DU)
    //3- The type of problem to be solved:
    //   1- Stokes problem;   
    //   2- Navier-Stokes problem (solves the steady Stokes
    //      in the first step to perform the initial guess);
    //4- 0 - Steady problem 
    //   1 - Transient problem  
      
    coarseModel.readInitialValues("initialVelocityCoarse2.txt","initialPressureCoarse2.txt");
    fineModel.readInitialValues("initialVelocityFine2.txt","initialPressureFine2.txt");


       
     arlequinProblem.setFluidModels(coarseModel, fineModel) ; 

     arlequinProblem.solveArlequinProblem(3, 1.e-7, 2, 1); 
          
     
    //Finalize main program   
    PetscFinalize();
 
    return 0; 
}
 
 
  




 
