#include <fstream>
#include <LinearAnalysis.h>
#include <GmshTools.h>
#include <VTUGenerator.h>
#include <L2Projection.h>
#include <Elasticity2D.h>
#include <memory>
#include <catch2/catch.hpp>
#define fTolerance 1.e-10
#define CATCH_CONFIG_MAIN
using namespace Catch::literals;

const double Elas = 1.;
const double nu = 0.;

auto exactSol2D = [](const VecDouble &coord, VecDouble &u, MatrixDouble &gradU){
    const auto &x=coord[0];
    const auto &y=coord[1];

    u[0] = x*(x-1.)+y*(y-1.);
    u[1] = x*(x-1.)+y*(y-1.);
    gradU(0,0) = 2.*x-1.;
    gradU(1,0) = 2.*x-1.;
    gradU(0,1) = 2.*y-1.;
    gradU(1,1) = 2.*y-1.;
};

auto forcing2D = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];

    force[0] = -(Elas*(-5.+nu))/(2.*(-1.+nu*nu));
    force[1] = -(Elas*(-5.+nu))/(2.*(-1.+nu*nu));
};

void SolveProblem2D(GeoMesh *gmesh, int order){
    CompMesh * cmesh = new CompMesh(gmesh,ApproxType::EIsoparametric);
    cmesh->Dimension() = 2;

    Elasticity2D * mat = new Elasticity2D(6,Elas,nu,true,1.);

    mat->SetForcingFunction(forcing2D);
    mat->SetExactSolution(exactSol2D);
    cmesh->InsertMaterial(mat);

    //BC
    MatrixDouble val1(2,2);
    val1.setZero();
    VecDouble val2(2);
    val2.setZero();
    
    // Homogeneous Dirichlet
    L2Projection * matbc1 = new L2Projection(5,1,BoundaryConditionType::kDirichlet,val1,val2);
    matbc1->SetExactSolution(exactSol2D);

    cmesh->InsertMaterial(matbc1);
 
    cmesh->AutoBuild();
    cmesh->Print("cmesh.txt");

    LinearAnalysis an(cmesh,SolverType::ELDLt);
   
    std::vector<std::string> ScalarNames, VectorNames;
    // ScalarNames = {"Solution","ExactSolution"};
    VectorNames = {"Displacement","ExactDisplacement","ExactForce"};
    an.Run();

    VTUGenerator::PrintResults(cmesh,"resultPoisson",ScalarNames,VectorNames); 

    VecDouble errors(4);
    an.PostProcessError(errors);
    REQUIRE(errors[0]<fTolerance);
    REQUIRE(errors[1]<fTolerance);
    REQUIRE(errors[2]<fTolerance);
};

TEST_CASE("Poisson_test","[Poisson]")
{
    SECTION("Check Elasticity 2D"){
        //Geometric Mesh
        GeoMesh * gmesh1 = new GeoMesh();
        GmshTools::Read(*gmesh1,"../../UnitTest/placa.msh");
        //gmesh1->Print("gmesh.txt");
        SolveProblem2D(gmesh1, 1);
    }

}