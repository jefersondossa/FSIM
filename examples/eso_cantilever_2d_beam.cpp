#include <fstream>
#include <LinearAnalysis.h>
#include <GmshTools.h>
#include <VTUGenerator.h>
#include <L2Projection.h>
#include <Elasticity2D.h>
#include <memory>
#include <eso.hpp>

int main()
{
    std::unique_ptr<CompMesh> model = std::make_unique<CompMesh>();

    constexpr auto kVolumeMatId = 15;

    auto *govEquation = new Elasticity2D(kVolumeMatId, 1e3, 0.3);
    model->InsertMaterial(govEquation);

    MatrixDouble val1(2, 2);
    val1.setZero();
    VecDouble val2(2);
    val2.setZero();
    val2[0] = 0.;

    constexpr auto kEngasteMatId = 17;

    auto *engasteBC = new L2Projection(kEngasteMatId, 2, BoundaryConditionType::kDirichlet, val1, val2);
    model->InsertMaterial(engasteBC);

    val1.setZero();
    val2.setZero();
    val2[0] = 0.;

    constexpr auto kFreeMatId = 18;

    auto *freeBC = new L2Projection(kFreeMatId, 2, BoundaryConditionType::kNeumann, val1, val2);
    model->InsertMaterial(freeBC);

    val2.setZero();
    val2[1] = -10;

    constexpr auto kLoadMatId = 16;
    auto *loadBC = new L2Projection(kLoadMatId, 2, BoundaryConditionType::kNeumann, val1, val2);
    model->InsertMaterial(loadBC);

    GmshTools::Read(*model, "../cantilever_2d_beam.msh");

    auto an = std::make_unique<LinearAnalysis>(model.get(), SolverType::AMGCLBiCGStab);

    RunEso(*model, *an, 0.7);

    return 0;
}