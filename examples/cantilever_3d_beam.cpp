#include <fstream>
#include <LinearAnalysis.h>
#include <GmshTools.h>
#include <VTUGenerator.h>
#include <L2Projection.h>
#include <Elasticity3D.h>
#include <memory>

int main()
{
    std::unique_ptr<CompMesh> model = std::make_unique<CompMesh>();

    constexpr auto kVolumeMatId = 15;

    auto *govEquation = new Elasticity3D(kVolumeMatId, 1e3, 0.3);
    model->InsertMaterial(govEquation);

    MatrixDouble val1(3, 3);
    val1.setZero();
    VecDouble val2(3);
    val2.setZero();
    val2[0] = 0.;

    constexpr auto kEngasteMatId = 17;

    auto *engasteBC = new L2Projection(kEngasteMatId, 3, BoundaryConditionType::kDirichlet, val1, val2);
    model->InsertMaterial(engasteBC);

    val1.setZero();
    val2.setZero();
    val2[0] = 0.;

    constexpr auto kFreeMatId = 18;

    auto *freeBC = new L2Projection(kFreeMatId, 3, BoundaryConditionType::kNeumann, val1, val2);
    model->InsertMaterial(freeBC);

    val2.setZero();
    val2[1] = -10;

    constexpr auto kLoadMatId = 16;
    auto *loadBC = new L2Projection(kLoadMatId, 3, BoundaryConditionType::kNeumann, val1, val2);
    model->InsertMaterial(loadBC);

    GmshTools::Read(*model, "../cantilever_3d_beam.msh");

    LinearAnalysis an(model.get(), SolverType::ECholmod);

    an.Run();

    std::vector<std::string> ScalarNames, VectorNames;
    VectorNames = {"Displacement", "Stress"};
    VTUGenerator::PrintResults(model.get(), "cantilever_3d_beam", ScalarNames, VectorNames);

    return 0;
}