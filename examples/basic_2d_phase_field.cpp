#include <fstream>
#include <NonLinearAnalysis.h>
#include <GmshTools.h>
#include <VTUGenerator.h>
#include <L2Projection.h>
#include <PhaseField.h>
#include <memory>

const auto exactSolPoisson = [](const VecDouble &coord, VecDouble &u, MatrixDouble &gradU){
    const auto &x=coord[0];
    const auto &y=coord[1];

    // u[0] = x * x * (x-1.) * y * y * (y-1.);
    // gradU(0,0) = x*(3*x-2.)*(y-1.)*y*y;
    // gradU(1,0) = x*x*(x-1.)*y*(3.*y-2.);
    // u[0] = x;
    // gradU(0,0) = 0.;
    // gradU(1,0) = 0.;

    u[0] = y - sinh(y)/sinh(1);
};

auto forcingFunctionPoisson = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    force[0] = y;
    // force[0] = 1.-x*x;
    // force[0] = -2. * (x*x*(1.-3.*y) - (y-1.)*y*y + 3.*x*(y-1.)*y*y + x*x*x*(3.*y-1.));
};

int main()
{
    std::unique_ptr<CompMesh> model = std::make_unique<CompMesh>();

    constexpr auto kInternalMatId = 15;

    auto *govEquation = new PhaseField(kInternalMatId, 2);
    govEquation->SetForcingFunction([](const VecDouble &, VecDouble &){});
    // govEquation->SetExactSolution(exactSolPoisson);
    model->InsertMaterial(govEquation);

    MatrixDouble val1(1, 1);
    VecDouble val2(1);
    val1.setZero();
    val2.setZero();

    constexpr auto kLeftMatId = 16;
    auto *leftBC = new L2Projection(kLeftMatId, 2, BoundaryConditionType::kNeumann, val1, val2);
    model->InsertMaterial(leftBC);

    val1.setZero();
    val2.setZero();
    val2[0] = 0;

    constexpr auto kRightMatId = 17;
    auto *rightBC = new L2Projection(kRightMatId, 2, BoundaryConditionType::kNeumann, val1, val2);
    model->InsertMaterial(rightBC);

    val2.setZero();
    val1.setZero();
    val2[0] = 0;

    constexpr auto kTopMatId = 18;
    auto *topBC = new L2Projection(kTopMatId, 2, BoundaryConditionType::kDirichlet, val1, val2);
    model->InsertMaterial(topBC);

    val2.setZero();
    val2[0] = 1;

    constexpr auto kBotttomMatId = 19;
    auto *bottomBC = new L2Projection(kBotttomMatId, 2, BoundaryConditionType::kDirichlet, val1, val2);
    model->InsertMaterial(bottomBC);

    GmshTools::Read(*model, "../rectangle.msh");

    NonLinearAnalysis an(model.get(), SolverType::AMGCLBiCGStab, 1e-6, 1000);

    an.Run();

    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"Solution", "ExactSolution"};
    VectorNames = {"Derivative"};
    VTUGenerator::PrintResults(model.get(), "basic_2d_phase_field", ScalarNames, VectorNames, {}, 1000);

    return 0;
}