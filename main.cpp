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
#include "FSInteraction.hpp"  


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
    
    //Data reading     
    // Cylinder
    // fineModel.dataReading("cyl_fine.txt","mirror_fine.txt");  
    // coarseModel.dataReading("cyl_coarse2.txt","mirror_coarse.txt");
    // fineModel.dataReading("cyl_fine_str4.txt","mirror_fine.txt");  
    // coarseModel.dataReading("cyl_coarse_str.txt","mirror_coarse.txt");
    fineModel.dataReading("various_fine.txt","mirror_fine.txt");  
    coarseModel.dataReading("various_coarse.txt","mirror_coarse.txt");
  
    // Cavity  
    // fineModel.dataReading("fine2.txt","mirror_fine.txt"); 
    // coarseModel.dataReading("coarse.txt","mirror_coarse.txt");
  
    // Helice 
    // fineModel.dataReading("cavpropfcir.txt","mirror_fine.txt"); 
    // coarseModel.dataReading("cavity1.txt","mirror_coarse.txt");

    // Flutter
    // fineModel.dataReading("flutter_fine2.txt","mirror_fine.txt"); 
    // coarseModel.dataReading("flutter_coarse.txt","mirror_coarse.txt");
     
    control.dataReading("32x32.txt","mirror_control.txt");

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
      
       
     arlequinProblem.setFluidModels(coarseModel, fineModel) ; 
       
     arlequinProblem.solveArlequinProblem(1, 1.e-7, 2, 1);
          
     
    //Finalize main program   
    PetscFinalize();
 
    return 0; 
}
 
 
  




 
