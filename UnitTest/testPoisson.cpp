#include <fstream>
#include <LinearAnalysis.h>
#include <GmshTools.h>
#include <VTUGenerator.h>
#include <L2Projection.h>
#include <Poisson.h>
#include <memory>
#include <catch2/catch.hpp>
#define fTolerance 1.e-10
#define CATCH_CONFIG_MAIN
using namespace Catch::literals;

auto forcing = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    force[0] = std::sin(M_PI*x/1);
};

auto exactSol = [](const VecDouble &coord, VecDouble &u, MatrixDouble &gradU){
    const auto &x=coord[0];
    const auto &y=coord[1];

    u[0] = std::sin(M_PI*x/1)/ (M_PI * M_PI);
    gradU(0,0) = std::cos(M_PI*x/1)*1/(M_PI);
};

auto forcing1 = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    force[0] = 0;
};

auto exactSol1 = [](const VecDouble &coord, VecDouble &u, MatrixDouble &gradU){
    const auto &x=coord[0];
    const auto &y=coord[1];

    u[0] = x;
    gradU(0,0) = 1.;
};

auto forcing2 = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    force[0] = -1.;
};

auto exactSol2 = [](const VecDouble &coord, VecDouble &u, MatrixDouble &gradU){
    const auto &x=coord[0];
    const auto &y=coord[1];

    u[0] = (x*x-x)/2.;
    gradU(0,0) = (x-0.5);
};

auto forcing3 = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    force[0] = x;
};

auto exactSol3 = [](const VecDouble &coord, VecDouble &u, MatrixDouble &gradU){
    const auto &x=coord[0];
    const auto &y=coord[1];

    u[0] = ((-x*x*x/6.)+(x/6.));
    gradU(0,0) = (1./6.*(1-3*x*x));
};

auto exactSol2D = [](const VecDouble &coord, VecDouble &u, MatrixDouble &gradU){
    const auto &x=coord[0];
    const auto &y=coord[1];

    u[0] = (x*(x-1.)*y*(y-1.));
    gradU(0,0) = ((-1.+2.*x)*(-1.+y)*y);
    gradU(1,0) = ((-1.+x)*x*(-1.+2.*y));
    // vs: u[0], solução do problema em questão;
    // vs: gradU(0,0), derivada de U em relação a x;
    // vs: gradU(0,1), derivada de U em relação a y;
};

auto forcing2D = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    force[0] = -2.*(-1.+x)*x+-2.*(-1.+y)*y;
    // -2.*(-1.+x)*x+-2.*(-1.+y)*y
    // vs: Laplaciano de U[0] = grad(0,0) + grad(0,1)
};

auto exactSol3D = [](const VecDouble &coord, VecDouble &u, MatrixDouble &gradU){
    const auto &x=coord[0];
    const auto &y=coord[1];

    u[0] = (x*(x-1.)*y*(y-1.));
    gradU(0,0) = ((-1.+2.*x)*(-1.+y)*y);
    gradU(1,0) = ((-1.+x)*x*(-1.+2.*y));
    // vs: u[0], solução do problema em questão;
    // vs: gradU(0,0), derivada de U em relação a x;
    // vs: gradU(0,1), derivada de U em relação a y;
};

auto forcing3D = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    force[0] = -2.*(-1.+x)*x+-2.*(-1.+y)*y;
    // -2.*(-1.+x)*x+-2.*(-1.+y)*y
    // vs: Laplaciano de U[0] = grad(0,0) + grad(0,1)
};

void SolveProblemHarmonic(GeoMesh *gmesh){
    CompMesh* cmesh = new CompMesh(gmesh,ApproxType::EIsoparametric);
    cmesh->Dimension() = 1;

    Poisson * mat = new Poisson(6,1,1);
    mat->SetForcingFunction(forcing);
    mat->SetExactSolution(exactSol);
    cmesh->InsertMaterial(mat);

    //BC
    MatrixDouble val1(1,1);
    val1.setZero();
    VecDouble val2(1);
    val2.setZero();
    
    //Apoio fixo
    L2Projection * matbc1 = new L2Projection(3,0,BoundaryConditionType::kDirichlet,val1,val2);

    //Apoio fixo
    L2Projection * matbc2 = new L2Projection(4,0,BoundaryConditionType::kDirichlet,val1,val2);

    cmesh->InsertMaterial(matbc1);
    cmesh->InsertMaterial(matbc2);
 
    cmesh->AutoBuild();
    cmesh->Print("cmesh.txt");

    LinearAnalysis an(cmesh,SolverType::ELDLt);
   
    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"Solution","ExactSolution"};
    VectorNames = {"Derivative"};
    an.Run();

    //Put a check criterion here.

    VTUGenerator::PrintResults(cmesh,"resultPoisson",ScalarNames,VectorNames); 

    VecDouble errors(3);
    an.PostProcessError(errors);
};

void SolveProblem(GeoMesh *gmesh, int order){
    CompMesh* cmesh = new CompMesh(gmesh,ApproxType::EIsoparametric);
    cmesh->Dimension() = 1;

    Poisson * mat = new Poisson(6,1,1);
    switch (order)
    {
    case 1:
        mat->SetForcingFunction(forcing1);
        mat->SetExactSolution(exactSol1);
        break;

    case 2:
        mat->SetForcingFunction(forcing2);
        mat->SetExactSolution(exactSol2);
        break;

    case 3:
        mat->SetForcingFunction(forcing3);
        mat->SetExactSolution(exactSol3);
        break;
    
    default:
        std::cout<<"Order not implemented\n";
        break;
    }
    cmesh->InsertMaterial(mat);

    //BC
    MatrixDouble val1(1,1);
    val1.setZero();
    VecDouble val2(1);
    val2.setZero();
    
    //Apoio fixo
    L2Projection * matbc1 = new L2Projection(3,0,BoundaryConditionType::kDirichlet,val1,val2);

    //Apoio fixo
    if (order == 1){
        val2[0]=1;
    }
    L2Projection * matbc2 = new L2Projection(4,0,BoundaryConditionType::kDirichlet,val1,val2);

    cmesh->InsertMaterial(matbc1);
    cmesh->InsertMaterial(matbc2);
 
    cmesh->AutoBuild();
    cmesh->Print("cmesh.txt");

    LinearAnalysis an(cmesh,SolverType::ELDLt);
   
    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"Solution","ExactSolution"};
    VectorNames = {"Derivative"};
    an.Run();

    // VTUGenerator::PrintResults(cmesh,"resultPoisson",ScalarNames,VectorNames); 

    VecDouble errors(3);
    an.PostProcessError(errors);
    REQUIRE(errors[0]<fTolerance);
    REQUIRE(errors[1]<fTolerance);
    REQUIRE(errors[2]<fTolerance);
};

void SolveProblemHierarquic(GeoMesh *gmesh, int order){
    CompMesh* cmesh = new CompMesh(gmesh,ApproxType::EHierarquic);
    cmesh->Dimension() = 1;
    cmesh->SetDefaultOrder(order);

    Poisson * mat = new Poisson(6,1,1);
    switch (order)
    {
    case 1:
        mat->SetForcingFunction(forcing1);
        mat->SetExactSolution(exactSol1);
        break;

    case 2:
        mat->SetForcingFunction(forcing2);
        mat->SetExactSolution(exactSol2);
        break;

    case 3:
        mat->SetForcingFunction(forcing3);
        mat->SetExactSolution(exactSol3);
        break;
    
    default:
        std::cout<<"Order not implemented\n";
        break;
    }

    cmesh->InsertMaterial(mat);

    //BC
    MatrixDouble val1(1,1);
    val1.setZero();
    VecDouble val2(1);
    val2.setZero();
    
    //Apoio fixo
    L2Projection * matbc1 = new L2Projection(3,0,BoundaryConditionType::kDirichlet,val1,val2);

    //Apoio fixo
    if (order == 1){
        val2[0]=1;
    }
    L2Projection * matbc2 = new L2Projection(4,0,BoundaryConditionType::kDirichlet,val1,val2);

    cmesh->InsertMaterial(matbc1);
    cmesh->InsertMaterial(matbc2);
 
    cmesh->AutoBuild();
    cmesh->Print("cmesh.txt");

    LinearAnalysis an(cmesh,SolverType::ELDLt);
   
    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"Solution","ExactSolution"};
    VectorNames = {"Derivative"};
    an.Run();

    //Put a check criterion here.

    // VTUGenerator::PrintResults(cmesh,"resultPoisson",ScalarNames,VectorNames); 

    VecDouble errors(3);
    an.PostProcessError(errors);
    REQUIRE(errors[0]<fTolerance);
    REQUIRE(errors[1]<fTolerance);
    REQUIRE(errors[2]<fTolerance);
};

void SolveProblem2D(GeoMesh *gmesh, int order){
    CompMesh * cmesh = new CompMesh(gmesh,ApproxType::EIsoparametric);
    cmesh->Dimension() = 2;

    Poisson * mat = new Poisson(6,2,1);

    mat->SetForcingFunction(forcing2D);
    mat->SetExactSolution(exactSol2D);
    cmesh->InsertMaterial(mat);

    //BC
    MatrixDouble val1(1,1);
    val1.setZero();
    VecDouble val2(1);
    val2.setZero();
    
    // Homogeneous Dirichlet
    L2Projection * matbc1 = new L2Projection(5,1,BoundaryConditionType::kDirichlet,val1,val2);

    cmesh->InsertMaterial(matbc1);
 
    cmesh->AutoBuild();
    cmesh->Print("cmesh.txt");

    LinearAnalysis an(cmesh,SolverType::ELDLt);
   
    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"Solution","ExactSolution"};
    VectorNames = {"Derivative"};
    an.Run();

    VTUGenerator::PrintResults(cmesh,"resultPoisson",ScalarNames,VectorNames); 

    VecDouble errors(3);
    an.PostProcessError(errors);
    REQUIRE(errors[0]<fTolerance);
    REQUIRE(errors[1]<fTolerance);
    REQUIRE(errors[2]<fTolerance);
};

void SolveProblem2DHierarquic(GeoMesh *gmesh, int order){
    CompMesh * cmesh = new CompMesh(gmesh,ApproxType::EHierarquic);
    cmesh->SetDefaultOrder(2);
    cmesh->Dimension() = 2;

    Poisson * mat = new Poisson(6,2,1);

    mat->SetForcingFunction(forcing2D);
    mat->SetExactSolution(exactSol2D);
    cmesh->InsertMaterial(mat);

    //BC
    MatrixDouble val1(1,1);
    val1.setZero();
    VecDouble val2(1);
    val2.setZero();
    
    // Homogeneous Dirichlet
    L2Projection * matbc1 = new L2Projection(5,1,BoundaryConditionType::kDirichlet,val1,val2);

    cmesh->InsertMaterial(matbc1);
 
    cmesh->AutoBuild();
    cmesh->Print("cmesh.txt");

    LinearAnalysis an(cmesh,SolverType::ELDLt);
   
    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"Solution","ExactSolution"};
    VectorNames = {"Derivative"};
    an.Run();

    VTUGenerator::PrintResults(cmesh,"resultPoisson",ScalarNames,VectorNames); 

    VecDouble errors(3);
    an.PostProcessError(errors);
    REQUIRE(errors[0]<fTolerance);
    REQUIRE(errors[1]<fTolerance);
    REQUIRE(errors[2]<fTolerance);
};

void SolveProblem3D(GeoMesh *gmesh, int order){
    //TODO Implement 3d test
    CompMesh * cmesh = new CompMesh(gmesh,ApproxType::EIsoparametric);
    cmesh->Dimension() = 3;

    Poisson * mat = new Poisson(14,3,1);

    mat->SetForcingFunction(forcing3D);
    mat->SetExactSolution(exactSol3D);
    cmesh->InsertMaterial(mat);

    //BC
    MatrixDouble val1(1,1);
    val1.setZero();
    VecDouble val2(1);
    val2.setZero();
    
    // Homogeneous Dirichlet
    L2Projection * matbc1 = new L2Projection(13,2,BoundaryConditionType::kDirichlet,val1,val2);

    cmesh->InsertMaterial(matbc1);
 
    cmesh->AutoBuild();
    cmesh->Print("cmesh.txt");

    LinearAnalysis an(cmesh,SolverType::ELDLt);
   
    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"Solution","ExactSolution"};
    VectorNames = {"Derivative"};
    an.Run();
    an.PrintGlobalMatrix();
    an.PrintGlobalRhs();

    VTUGenerator::PrintResults(cmesh,"resultPoisson",ScalarNames,VectorNames); 

    VecDouble errors(3);
    an.PostProcessError(errors);
    REQUIRE(errors[0]<fTolerance);
    REQUIRE(errors[1]<fTolerance);
    REQUIRE(errors[2]<fTolerance);
};

TEST_CASE("Poisson_test","[Poisson]")
{
    SECTION("Check Isoparametric"){
        //Geometric Mesh
        GeoMesh * gmesh1 = new GeoMesh();
        GmshTools::Read(*gmesh1,"../../UnitTest/poisson1d-1.msh");
        //gmesh1->Print("gmesh.txt");

        SolveProblem(gmesh1, 1);

        //Geometric Mesh
        GeoMesh * gmesh2 = new GeoMesh();
        GmshTools::Read(*gmesh2,"../../UnitTest/poisson1d-2.msh");
        //gmesh1->Print("gmesh.txt");

        SolveProblem(gmesh2, 2);

        //Geometric Mesh
        GeoMesh * gmesh3 = new GeoMesh();
        GmshTools::Read(*gmesh3,"../../UnitTest/poisson1d-3.msh");
        //gmesh1->Print("gmesh.txt");

        SolveProblem(gmesh3, 3);
    }

    SECTION("Check Hierarquic"){
        //Geometric Mesh
        GeoMesh * gmesh1 = new GeoMesh();
        GmshTools::Read(*gmesh1,"../../UnitTest/poisson1d-1.msh");
        //gmesh1->Print("gmesh.txt");

        for (int iorder = 1; iorder < 3; iorder++){
            SolveProblemHierarquic(gmesh1, iorder);    
        }
    }

    SECTION("Check 2D"){
        //Geometric Mesh
        GeoMesh * gmesh1 = new GeoMesh();
        GmshTools::Read(*gmesh1,"../../UnitTest/placa.msh");
        //gmesh1->Print("gmesh.txt");
        SolveProblem2D(gmesh1, 1);
    }

     SECTION("Check 2DHierarquic"){
        //Geometric Mesh
        GeoMesh * gmesh1 = new GeoMesh();
        GmshTools::Read(*gmesh1,"../../UnitTest/placa.msh");
        //gmesh1->Print("gmesh.txt");
        SolveProblem2DHierarquic(gmesh1, 1);
    }

}

/*  
int main()
{

Geometric Mesh
GeoMesh * gmesh = new GeoMesh();
GmshTools::Read(*gmesh,"../../UnitTest/poisson1d.msh");
gmesh->Print("gmesh.txt");

SolveProblem1(gmesh);
SolveProblem2(gmesh);
SolveProblem3(gmesh);
SolveProblemHarmonic(gmesh);

return 0;

}
*/