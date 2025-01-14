#include <fstream>
#include <NonLinearAnalysis.h>
#include <GmshTools.h>
#include <VTUGenerator.h>
#include <L2Projection.h>
#include <PhaseField.h>
#include <memory>

int main()
{
    std::unique_ptr<CompMesh> model = std::make_unique<CompMesh>();

    constexpr auto kInternalMatId = 15;

    auto *govEquation = new PhaseField(kInternalMatId, 2);
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

    NonLinearAnalysis an(model.get(), SolverType::ELU, 1e-3, 1000);

    an.Run();

    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"Solution"};
    VectorNames = {"Derivative"};
    VTUGenerator::PrintResults(model.get(), "basic_2d_phase_field", ScalarNames, VectorNames);

    return 0;
}