#include <fstream>
#include <TransientAnalysis.h>
#include <GmshTools.h>
#include <VTUGenerator.h>
#include <VTUReader.h>
#include <L2Projection.h>
#include <PhaseField.h>
#include <Elasticity2D.h>
#include <ElementTransient.h>
#include <ShapePoint.h>
#include <Node.h>
#include <memory>
#include <LagrangeMultiplier.h>
#include <ElementLagrangeMultiplier.h>
#include <ShapeQuadrilateralLin.h>
#include <ShapeTriangleLin.h>
#include <NullWeakForm.h>
#include <algorithm>
#include <functional>
#include <stdexcept>
#include <string>

static constexpr double V_F = 0.55;

namespace
{
bool file_exists(const std::string &filename)
{
    std::ifstream input(filename);
    return static_cast<bool>(input);
}

std::string resolve_vtu_filename(const std::string &filename_or_prefix, int initial_step)
{
    if (file_exists(filename_or_prefix))
    {
        return filename_or_prefix;
    }

    const auto indexed_filename = filename_or_prefix + std::to_string(initial_step) + ".vtu";
    if (file_exists(indexed_filename))
    {
        return indexed_filename;
    }

    throw std::runtime_error("could not find VTU file from \"" + filename_or_prefix + "\"");
}

double sample_row_wise_grid(const VTUReader::Grid &grid, double x, double y, double lx, double ly)
{
    if (grid.empty() || grid.front().empty())
    {
        throw std::runtime_error("initial VTU field is empty");
    }
    if (lx <= 0.0 || ly <= 0.0)
    {
        throw std::runtime_error("LX and LY must be positive");
    }

    const std::size_t nrows = grid.size();
    const std::size_t ncols = grid.front().size();
    for (const auto &row : grid)
    {
        if (row.size() != ncols)
        {
            throw std::runtime_error("initial VTU field is not rectangular");
        }
    }

    const double x_clamped = std::clamp(x, 0.0, lx);
    const double y_clamped = std::clamp(y, 0.0, ly);

    const double col_pos = (ncols == 1) ? 0.0 : x_clamped * static_cast<double>(ncols - 1) / lx;
    const double row_pos = (nrows == 1) ? 0.0 : y_clamped * static_cast<double>(nrows - 1) / ly;

    const std::size_t col0 = static_cast<std::size_t>(std::floor(col_pos));
    const std::size_t row0 = static_cast<std::size_t>(std::floor(row_pos));
    const std::size_t col1 = std::min(col0 + 1, ncols - 1);
    const std::size_t row1 = std::min(row0 + 1, nrows - 1);

    const double tx = col_pos - static_cast<double>(col0);
    const double ty = row_pos - static_cast<double>(row0);

    const double v00 = grid[row0][col0];
    const double v01 = grid[row0][col1];
    const double v10 = grid[row1][col0];
    const double v11 = grid[row1][col1];

    const double top = (1.0 - tx) * v00 + tx * v01;
    const double bottom = (1.0 - tx) * v10 + tx * v11;
    return (1.0 - ty) * top + ty * bottom;
}
}

void SetupBoundaryConditionsElasticity2D2(CompMesh &modelElasticity2D)
{
    MatrixDouble val1(2, 2);
    VecDouble val2(2);
    constexpr auto kVolumeMatId = 15;
    auto *govEquationElasticity2D = new Elasticity2D(kVolumeMatId, 1e3, 0.3);
    modelElasticity2D.InsertMaterial(govEquationElasticity2D);

    val1.setZero();
    val2.setZero();
    val2[0] = 0.;

    constexpr auto kEngasteMatId = 16;

    auto *engasteBC = new L2Projection(kEngasteMatId, 2, BoundaryConditionType::kDirichlet, val1, val2);
    modelElasticity2D.InsertMaterial(engasteBC);

    val1.setZero();
    val2.setZero();
    val2[0] = 0.;

    constexpr auto kFreeMatTopId = 18;

    auto *freeBC = new L2Projection(kFreeMatTopId, 2, BoundaryConditionType::kNeumann, val1, val2);
    modelElasticity2D.InsertMaterial(freeBC);

    val1.setZero();
    val2.setZero();
    val2[0] = 0.;

    constexpr auto kFreeMatBottomId = 19;

    auto *freeBCEl2D = new L2Projection(kFreeMatBottomId, 2, BoundaryConditionType::kNeumann, val1, val2);
    modelElasticity2D.InsertMaterial(freeBCEl2D);

    val2.setZero();
    val2[0] = 10;
    constexpr auto kLoadMatId = 45;
    auto *El2D = new L2Projection(kLoadMatId, 2, BoundaryConditionType::kNeumann, val1, val2);
    modelElasticity2D.InsertMaterial(El2D);

    GmshTools::Read(modelElasticity2D, "../../rectangle.msh");

    // TODO: Disabled because now we recalculate the stiffness matrix contribution per node (based on phi)
    // Disables memory on elements
    // (makes sure elemental stiffness is not recalculated)
    // for(auto& [_, pWeakForm] : modelElasticity2D.MaterialVector()) {
    // govEquationElasticity2D->SetHasMemory(false);
    // }
}

void SetupBoundaryConditionsElasticity2DCantileverRightBottom(CompMesh& modelElasticity2D)
{
    MatrixDouble val1(2, 2);
    VecDouble val2(2);
    constexpr auto kVolumeMatId = 15;
    auto *govEquationElasticity2D = new Elasticity2D(kVolumeMatId, 1e3, 0.3);
    modelElasticity2D.InsertMaterial(govEquationElasticity2D);
    val1.setZero();
    val2.setZero();
    constexpr auto kEngasteMatId = 18;
    auto *engasteBC = new L2Projection(kEngasteMatId, 2, BoundaryConditionType::kDirichlet, val1, val2);
    modelElasticity2D.InsertMaterial(engasteBC);
    val1.setZero();
    val2.setZero();
    val2[1] = -10;
    constexpr auto kLoadMatId = 17;
    auto *El2D2 = new L2Projection(kLoadMatId, 2, BoundaryConditionType::kNeumann, val1, val2);
    modelElasticity2D.InsertMaterial(El2D2);
    GmshTools::Read(modelElasticity2D, "../../rectangle.msh");

    // TODO: Disabled because now we recalculate the stiffness matrix contribution per node (based on phi)
    // Disables memory on elements
    // (makes sure elemental stiffness is not recalculated)
    // for(auto& [_, pWeakForm] : modelElasticity2D.MaterialVector()) {
        // govEquationElasticity2D->SetHasMemory(false);
    // }
}

void SetupBoundaryConditionsElasticity2DMBB(CompMesh &modelElasticity2D)
{
    MatrixDouble val1(2, 2);
    VecDouble val2(2);
    constexpr auto kVolumeMatId = 15;
    auto *govEquationElasticity2D = new Elasticity2D(kVolumeMatId, 1e3, 0.3);
    modelElasticity2D.InsertMaterial(govEquationElasticity2D);

    val1.setZero();
    val2.setZero();
    val2[0] = 1.0;
    constexpr auto kEngasteMatId = 16;
    auto *engasteBC = new L2Projection(kEngasteMatId, 2, BoundaryConditionType::kDirectionalHomogeneousDirichlet, val1, val2);
    modelElasticity2D.InsertMaterial(engasteBC);

    val1.setZero();
    val2.setZero();
    val2[1] = 1.0;
    constexpr auto kBottomRight = 17;
    auto *bottomRightBC = new L2Projection(kBottomRight, 2, BoundaryConditionType::kDirectionalHomogeneousDirichlet, val1, val2);
    modelElasticity2D.InsertMaterial(bottomRightBC);

    val2.setZero();
    val1.setZero();
    val2[0] = 0.;
    constexpr auto kFreeMatTopId = 18;
    auto *freeBC = new L2Projection(kFreeMatTopId, 2, BoundaryConditionType::kNeumann, val1, val2);
    modelElasticity2D.InsertMaterial(freeBC);

    val1.setZero();
    val2.setZero();
    val2[0] = 0.;
    constexpr auto kFreeMatBottomId = 19;
    auto *freeBCEl2D = new L2Projection(kFreeMatBottomId, 2, BoundaryConditionType::kNeumann, val1, val2);
    modelElasticity2D.InsertMaterial(freeBCEl2D);

    val1.setZero();
    val2.setZero();
    val2[1] = -1;
    constexpr auto kLoadMatId = 23;
    auto *El2D2 = new L2Projection(kLoadMatId, 2, BoundaryConditionType::kNeumann, val1, val2);
    modelElasticity2D.InsertMaterial(El2D2);
    GmshTools::Read(modelElasticity2D, "../../rectangle_mbb.msh");
}

void SetupBoundaryConditionsElasticity2DBridge(CompMesh &modelElasticity2D)
{
    MatrixDouble val1(2, 2);
    VecDouble val2(2);
    constexpr auto kVolumeMatId = 15;
    auto *govEquationElasticity2D = new Elasticity2D(kVolumeMatId, 1e3, 0.3);
    modelElasticity2D.InsertMaterial(govEquationElasticity2D);

    // val1.setZero();
    // val2.setZero();
    // constexpr auto kLeftMatId = 16;
    // auto *leftBC = new L2Projection(kLeftMatId, 2, BoundaryConditionType::kNeumann, val1, val2);
    // modelElasticity2D.InsertMaterial(leftBC);

    val1.setZero();
    val2.setZero();
    val2[1] = 1.0;
    constexpr auto kBottomRight = 17;
    auto *bottomRightBC = new L2Projection(kBottomRight, 2, BoundaryConditionType::kDirectionalHomogeneousDirichlet, val1, val2);
    modelElasticity2D.InsertMaterial(bottomRightBC);

    // val2.setZero();
    // val1.setZero();
    // val2[0] = 0.;
    // constexpr auto kFreeMatTopId = 18;
    // auto *freeBC = new L2Projection(kFreeMatTopId, 2, BoundaryConditionType::kNeumann, val1, val2);
    // modelElasticity2D.InsertMaterial(freeBC);

    // val1.setZero();
    // val2.setZero();
    // constexpr auto kFreeMatBottomId = 19;
    // auto *freeBCEl2D = new L2Projection(kFreeMatBottomId, 2, BoundaryConditionType::kNeumann, val1, val2);
    // modelElasticity2D.InsertMaterial(freeBCEl2D);

    val1.setZero();
    val2.setZero();
    val2[1] = -40;
    constexpr auto kMiddleBottom = 23;
    auto *El2D2 = new L2Projection(kMiddleBottom, 2, BoundaryConditionType::kNeumann, val1, val2);
    modelElasticity2D.InsertMaterial(El2D2);

    val1.setZero();
    val2.setZero();
    constexpr auto kLeftBottom = 24;
    auto *bottomLeftBC = new L2Projection(kLeftBottom, 2, BoundaryConditionType::kDirichlet, val1, val2);
    modelElasticity2D.InsertMaterial(bottomLeftBC);

    GmshTools::Read(modelElasticity2D, "../../rectangle_bridge.msh");
}

void SetupBoundaryConditionsPhaseFieldBridge(CompMesh &modelPhaseField)
{
    MatrixDouble val1(1, 1);
    VecDouble val2(1);

    val1.setZero();
    val2.setZero();
    constexpr auto kPhaseFieldInternalMatId = 15;
    auto *govEquationPF = new PhaseField(kPhaseFieldInternalMatId, 2);
    modelPhaseField.InsertMaterial(govEquationPF);

    // val1.setZero();
    // val2.setZero();
    // constexpr auto kEngasteMatId = 16;
    // auto *engasteBC = new L2Projection(kEngasteMatId, 2, BoundaryConditionType::kNeumann, val1, val2);
    // modelPhaseField.InsertMaterial(engasteBC);

    val1.setZero();
    val2.setZero();
    constexpr auto kBottomRight = 17;
    auto *bottomRightBC = new L2Projection(kBottomRight, 2, BoundaryConditionType::kNeumann, val1, val2);
    modelPhaseField.InsertMaterial(bottomRightBC);

    // val2.setZero();
    // val1.setZero();
    // constexpr auto kFreeMatTopId = 18;
    // auto *freeBC = new L2Projection(kFreeMatTopId, 2, BoundaryConditionType::kNeumann, val1, val2);
    // modelPhaseField.InsertMaterial(freeBC);

    // val1.setZero();
    // val2.setZero();
    // constexpr auto kFreeMatBottomId = 19;
    // auto *freeBCEl2D = new L2Projection(kFreeMatBottomId, 2, BoundaryConditionType::kNeumann, val1, val2);
    // modelPhaseField.InsertMaterial(freeBCEl2D);

    val1.setZero();
    val2.setZero();
    constexpr auto kLoadMatId = 23;
    auto *El2D2 = new L2Projection(kLoadMatId, 2, BoundaryConditionType::kNeumann, val1, val2);
    modelPhaseField.InsertMaterial(El2D2);

    val1.setZero();
    val2.setZero();
    constexpr auto kLeftBottom = 24;
    auto *leftBottomBC = new L2Projection(kLeftBottom, 2, BoundaryConditionType::kNeumann, val1, val2);
    modelPhaseField.InsertMaterial(leftBottomBC);

    do
    {
        break;
        constexpr auto kPhaseFieldLeftMatId = 16;
        auto *leftBCPF = new L2Projection(kPhaseFieldLeftMatId, 2, BoundaryConditionType::kNeumann, val1, val2);
        modelPhaseField.InsertMaterial(leftBCPF);
    
        constexpr auto kPhaseFieldRightMatId = 45;
        auto *rightBCPF = new L2Projection(kPhaseFieldRightMatId, 2, BoundaryConditionType::kNeumann, val1, val2);
        modelPhaseField.InsertMaterial(rightBCPF);
    
        constexpr auto kPhaseFieldTopMatId = 18;
        auto *topBCPF = new L2Projection(kPhaseFieldTopMatId, 2, BoundaryConditionType::kNeumann, val1, val2);
        modelPhaseField.InsertMaterial(topBCPF);
    
        // val2[0] = 1.0;
        constexpr auto kPhaseFieldBotttomMatId = 19;
        auto *bottomBCPF = new L2Projection(kPhaseFieldBotttomMatId, 2, BoundaryConditionType::kNeumann, val1, val2);
        modelPhaseField.InsertMaterial(bottomBCPF);
    
        // {
        //     constexpr auto kLoadMatId = 17;
        //     auto *loaddd = new L2Projection(kLoadMatId, 2, BoundaryConditionType::kNeumann, val1, val2);
        //     modelPhaseField.InsertMaterial(loaddd);
        // }
    
        // {
        //     constexpr auto kLoadMatId = 23;
        //     auto *loaddd = new L2Projection(kLoadMatId, 2, BoundaryConditionType::kNeumann, val1, val2);
        //     modelPhaseField.InsertMaterial(loaddd);
    } while (false);
    
    GmshTools::Read(modelPhaseField, "../../rectangle_bridge.msh");
    // return;
    // For the phase field implementation, we also need to impose an additional restrain to the final volume.
    // To do so, we create a new node, which serves to store the additional Lagrange multiplier variable.
    // However, the LM has the contribution of all elements. Thus, the new node need to be set as a new connectivity for all elements.
    // Note that this completely changes the data structure, and additional care is needed to handle it.
    // For instance, changes in the shape functions and assemble algorithm may be needed.
    // First, create a new node and a new element to store the LM.
    VecDouble Coor(3);
    Coor[0] = 0.0;
    Coor[1] = 0.0;
    Coor[2] = 0.0;
    Node *newNode = new Node(Coor, modelPhaseField.NodeVec().size(), 1);
    newNode->AllocateTimeDerivatives();
    modelPhaseField.NodeVec().push_back(newNode);
    int64_t index = newNode->Index();
    VecInt connect(1);
    connect[0] = index;
    // Point
    Element *gel = nullptr;
    int matnull = modelPhaseField.GetNewMaterialId();
    std::cout << "matnul " << matnull << std::endl;
    NullWeakForm *nullwf = new NullWeakForm(matnull, 1);
    gel = new ElementTransient<ShapePoint>(index, connect, &modelPhaseField, nullwf);
    modelPhaseField.InsertElement(gel);

    // Creates the Lagrange multiplier material
    int matlagmult = modelPhaseField.GetNewMaterialId();
    LagrangeMultiplier *lagmult = new LagrangeMultiplier(matlagmult, 1);
    // Sets the initial/final volume constraint as a forcing function
    //  const auto diam = 2*100.0/700;
    //  double finalVol = 2.-12.*M_PI*diam*diam/4.0;
    static constexpr auto diam = 2 * 100.0 / 700;
    static constexpr auto holesVol = 12. * M_PI * diam * diam / 4.0;
    static constexpr auto finalVol = 2. - holesVol;
    static constexpr auto finalVolRel = finalVol / 2.0;
    auto forcingFunction = [](const VecDouble &coord, VecDouble &force)
    {
        force[0] = V_F;
        // force[0] = 0.0;
    };
    lagmult->SetForcingFunction(forcingFunction);
    modelPhaseField.InsertMaterial(lagmult);
    // Now, create new elements of the type LagrangeMultiplier, which will be used to impose the additional constrain.
    //  Here, the point element will always be the second and the volumetric element, the first.
    //  std::cout << "NElements before LagrangeMultiplier: " << modelPhaseField.NElements() << '\n';
    for (auto &element : modelPhaseField.ElementVec())
    {
        if (element->Dimension() != modelPhaseField.Dimension())
            continue;

        Element *gelmult = nullptr;
        auto *wf = modelPhaseField.Material(matlagmult);
        auto newindex = modelPhaseField.NElements();
        auto type = element->Type();

        switch (type)
        {
        case ElementType::EQuadrilateral:
            if (modelPhaseField.GetDefaultOrder() == 1)
            {
                gelmult = new ElementLagrangeMultiplier<ShapeQuadrilateralLin>(newindex, element, gel, &modelPhaseField, wf);
            }
            else
            {
                PanicButton();
            }
            break;
        case ElementType::ETriangle:
            if (modelPhaseField.GetDefaultOrder() == 1)
            {
                gelmult = new ElementLagrangeMultiplier<ShapeTriangleLin>(newindex, element, gel, &modelPhaseField, wf);
            }
            else
            {
                PanicButton();
            }
            break;

        default:
            PanicButton();
            break;
        }
        modelPhaseField.InsertElement(gelmult);
    }
    // std::cout << "NElements after LagrangeMultiplier: " << modelPhaseField.NElements() << '\n';
    delete[] modelPhaseField.part_elem;
    modelPhaseField.part_elem = new int[modelPhaseField.NElements()]();
}

void SetupBoundaryConditionsPhaseFieldMBB(CompMesh &modelPhaseField)
{
    MatrixDouble val1(1, 1);
    VecDouble val2(1);

    val1.setZero();
    val2.setZero();
    constexpr auto kPhaseFieldInternalMatId = 15;
    auto *govEquationPF = new PhaseField(kPhaseFieldInternalMatId, 2);
    modelPhaseField.InsertMaterial(govEquationPF);


    val1.setZero();
    val2.setZero();
    constexpr auto kEngasteMatId = 16;
    auto *engasteBC = new L2Projection(kEngasteMatId, 2, BoundaryConditionType::kNeumann, val1, val2);
    modelPhaseField.InsertMaterial(engasteBC);

    val1.setZero();
    val2.setZero();
    constexpr auto kBottomRight = 17;
    auto *bottomRightBC = new L2Projection(kBottomRight, 2, BoundaryConditionType::kNeumann, val1, val2);
    modelPhaseField.InsertMaterial(bottomRightBC);

    val2.setZero();
    val1.setZero();
    constexpr auto kFreeMatTopId = 18;
    auto *freeBC = new L2Projection(kFreeMatTopId, 2, BoundaryConditionType::kNeumann, val1, val2);
    modelPhaseField.InsertMaterial(freeBC);

    val1.setZero();
    val2.setZero();
    constexpr auto kFreeMatBottomId = 19;
    auto *freeBCEl2D = new L2Projection(kFreeMatBottomId, 2, BoundaryConditionType::kNeumann, val1, val2);
    modelPhaseField.InsertMaterial(freeBCEl2D);

    val1.setZero();
    val2.setZero();
    constexpr auto kLoadMatId = 23;
    auto *El2D2 = new L2Projection(kLoadMatId, 2, BoundaryConditionType::kNeumann, val1, val2);
    modelPhaseField.InsertMaterial(El2D2);

    
    GmshTools::Read(modelPhaseField, "../../rectangle_mbb.msh");
    // return;
    // For the phase field implementation, we also need to impose an additional restrain to the final volume.
    // To do so, we create a new node, which serves to store the additional Lagrange multiplier variable.
    // However, the LM has the contribution of all elements. Thus, the new node need to be set as a new connectivity for all elements.
    // Note that this completely changes the data structure, and additional care is needed to handle it.
    // For instance, changes in the shape functions and assemble algorithm may be needed.
    // First, create a new node and a new element to store the LM.
    VecDouble Coor(3);
    Coor[0] = 0.0;
    Coor[1] = 0.0;
    Coor[2] = 0.0;
    Node *newNode = new Node(Coor, modelPhaseField.NodeVec().size(), 1);
    newNode->AllocateTimeDerivatives();
    modelPhaseField.NodeVec().push_back(newNode);
    int64_t index = newNode->Index();
    VecInt connect(1);
    connect[0] = index;
    // Point
    Element *gel = nullptr;
    int matnull = modelPhaseField.GetNewMaterialId();
    std::cout << "matnul " << matnull << std::endl;
    NullWeakForm *nullwf = new NullWeakForm(matnull, 1);
    gel = new ElementTransient<ShapePoint>(index, connect, &modelPhaseField, nullwf);
    modelPhaseField.InsertElement(gel);

    // Creates the Lagrange multiplier material
    int matlagmult = modelPhaseField.GetNewMaterialId();
    LagrangeMultiplier *lagmult = new LagrangeMultiplier(matlagmult, 1);
    // Sets the initial/final volume constraint as a forcing function
    //  const auto diam = 2*100.0/700;
    //  double finalVol = 2.-12.*M_PI*diam*diam/4.0;
    static constexpr auto diam = 2 * 100.0 / 700;
    static constexpr auto holesVol = 12. * M_PI * diam * diam / 4.0;
    static constexpr auto finalVol = 2. - holesVol;
    static constexpr auto finalVolRel = finalVol / 2.0;
    auto forcingFunction = [](const VecDouble &coord, VecDouble &force)
    {
        force[0] = V_F;
        // force[0] = 0.0;
    };
    lagmult->SetForcingFunction(forcingFunction);
    modelPhaseField.InsertMaterial(lagmult);
    // Now, create new elements of the type LagrangeMultiplier, which will be used to impose the additional constrain.
    //  Here, the point element will always be the second and the volumetric element, the first.
    //  std::cout << "NElements before LagrangeMultiplier: " << modelPhaseField.NElements() << '\n';
    for (auto &element : modelPhaseField.ElementVec())
    {
        if (element->Dimension() != modelPhaseField.Dimension())
            continue;

        Element *gelmult = nullptr;
        auto *wf = modelPhaseField.Material(matlagmult);
        auto newindex = modelPhaseField.NElements();
        auto type = element->Type();

        switch (type)
        {
        case ElementType::EQuadrilateral:
            if (modelPhaseField.GetDefaultOrder() == 1)
            {
                gelmult = new ElementLagrangeMultiplier<ShapeQuadrilateralLin>(newindex, element, gel, &modelPhaseField, wf);
            }
            else
            {
                PanicButton();
            }
            break;
        case ElementType::ETriangle:
            if (modelPhaseField.GetDefaultOrder() == 1)
            {
                gelmult = new ElementLagrangeMultiplier<ShapeTriangleLin>(newindex, element, gel, &modelPhaseField, wf);
            }
            else
            {
                PanicButton();
            }
            break;

        default:
            PanicButton();
            break;
        }
        modelPhaseField.InsertElement(gelmult);
    }
    // std::cout << "NElements after LagrangeMultiplier: " << modelPhaseField.NElements() << '\n';
    delete[] modelPhaseField.part_elem;
    modelPhaseField.part_elem = new int[modelPhaseField.NElements()]();
}

void SetupBoundaryConditionsPhaseFieldShear(CompMesh &modelPhaseField)
{
    MatrixDouble val1(1, 1);
    VecDouble val2(1);

    val1.setZero();
    val2.setZero();
    constexpr auto kPhaseFieldInternalMatId = 15;
    auto *govEquationPF = new PhaseField(kPhaseFieldInternalMatId, 2);
    modelPhaseField.InsertMaterial(govEquationPF);

    val2.setZero();
    val1.setZero();
    constexpr auto kFreeMatTopId = 18;
    auto *freeBC = new L2Projection(kFreeMatTopId, 2, BoundaryConditionType::kNeumann, val1, val2);
    modelPhaseField.InsertMaterial(freeBC);

    val1.setZero();
    val2.setZero();
    constexpr auto kFreeMatBottomId = 17;
    auto *freeBCEl2D = new L2Projection(kFreeMatBottomId, 2, BoundaryConditionType::kNeumann, val1, val2);
    modelPhaseField.InsertMaterial(freeBCEl2D);

    
    GmshTools::Read(modelPhaseField, "../../rectangle.msh");
    // return;
    // For the phase field implementation, we also need to impose an additional restrain to the final volume.
    // To do so, we create a new node, which serves to store the additional Lagrange multiplier variable.
    // However, the LM has the contribution of all elements. Thus, the new node need to be set as a new connectivity for all elements.
    // Note that this completely changes the data structure, and additional care is needed to handle it.
    // For instance, changes in the shape functions and assemble algorithm may be needed.
    // First, create a new node and a new element to store the LM.
    VecDouble Coor(3);
    Coor[0] = 0.0;
    Coor[1] = 0.0;
    Coor[2] = 0.0;
    Node *newNode = new Node(Coor, modelPhaseField.NodeVec().size(), 1);
    newNode->AllocateTimeDerivatives();
    modelPhaseField.NodeVec().push_back(newNode);
    int64_t index = newNode->Index();
    VecInt connect(1);
    connect[0] = index;
    // Point
    Element *gel = nullptr;
    int matnull = modelPhaseField.GetNewMaterialId();
    std::cout << "matnul " << matnull << std::endl;
    NullWeakForm *nullwf = new NullWeakForm(matnull, 1);
    gel = new ElementTransient<ShapePoint>(index, connect, &modelPhaseField, nullwf);
    modelPhaseField.InsertElement(gel);

    // Creates the Lagrange multiplier material
    int matlagmult = modelPhaseField.GetNewMaterialId();
    LagrangeMultiplier *lagmult = new LagrangeMultiplier(matlagmult, 1);
    // Sets the initial/final volume constraint as a forcing function
    //  const auto diam = 2*100.0/700;
    //  double finalVol = 2.-12.*M_PI*diam*diam/4.0;
    static constexpr auto diam = 2 * 100.0 / 700;
    static constexpr auto holesVol = 12. * M_PI * diam * diam / 4.0;
    static constexpr auto finalVol = 2. - holesVol;
    static constexpr auto finalVolRel = finalVol / 2.0;
    auto forcingFunction = [](const VecDouble &coord, VecDouble &force)
    {
        force[0] = V_F;
        // force[0] = 0.0;
    };
    lagmult->SetForcingFunction(forcingFunction);
    modelPhaseField.InsertMaterial(lagmult);
    // Now, create new elements of the type LagrangeMultiplier, which will be used to impose the additional constrain.
    //  Here, the point element will always be the second and the volumetric element, the first.
    //  std::cout << "NElements before LagrangeMultiplier: " << modelPhaseField.NElements() << '\n';
    for (auto &element : modelPhaseField.ElementVec())
    {
        if (element->Dimension() != modelPhaseField.Dimension())
            continue;

        Element *gelmult = nullptr;
        auto *wf = modelPhaseField.Material(matlagmult);
        auto newindex = modelPhaseField.NElements();
        auto type = element->Type();

        switch (type)
        {
        case ElementType::EQuadrilateral:
            if (modelPhaseField.GetDefaultOrder() == 1)
            {
                gelmult = new ElementLagrangeMultiplier<ShapeQuadrilateralLin>(newindex, element, gel, &modelPhaseField, wf);
            }
            else
            {
                PanicButton();
            }
            break;
        case ElementType::ETriangle:
            if (modelPhaseField.GetDefaultOrder() == 1)
            {
                gelmult = new ElementLagrangeMultiplier<ShapeTriangleLin>(newindex, element, gel, &modelPhaseField, wf);
            }
            else
            {
                PanicButton();
            }
            break;

        default:
            PanicButton();
            break;
        }
        modelPhaseField.InsertElement(gelmult);
    }
    // std::cout << "NElements after LagrangeMultiplier: " << modelPhaseField.NElements() << '\n';
    delete[] modelPhaseField.part_elem;
    modelPhaseField.part_elem = new int[modelPhaseField.NElements()]();
}

int main(int argc, char *argv[])
{
    bool use_vtu_initial_condition = false;
    int initial_step = 0;
    std::string initial_vtu_filename;
    std::string initial_field_name;
    double initial_lx = 0.0;
    double initial_ly = 0.0;
    VTUReader::Grid initial_field;

    if (argc == 6)
    {
        try
        {
            initial_vtu_filename = resolve_vtu_filename(argv[1], std::stoi(argv[5]));
            initial_field_name = argv[2];
            initial_lx = std::stod(argv[3]);
            initial_ly = std::stod(argv[4]);
            initial_step = std::stoi(argv[5]);
            initial_field = VTUReader::read_field_from_vtu(initial_vtu_filename, initial_field_name);
            if (initial_field.empty() || initial_field.front().empty())
            {
                throw std::runtime_error("failed to load initial field from VTU");
            }
            use_vtu_initial_condition = true;
            std::cout << "Initializing u0 from " << initial_vtu_filename << " using field \"" << initial_field_name << "\"\n";
        }
        catch (const std::exception &error)
        {
            std::cerr << "Invalid VTU initialization arguments: " << error.what() << '\n';
            return 1;
        }
    }
    else if (argc != 1)
    {
        std::cerr << "Usage: " << argv[0] << " [filename field_name lx ly initial_step]\n";
        return 1;
    }

    std::unique_ptr<CompMesh> modelElasticity2D = std::make_unique<CompMesh>();
    // SetupBoundaryConditionsElasticity2D2(*modelElasticity2D);
    // SetupBoundaryConditionsElasticity2DBridge(*modelElasticity2D);
    SetupBoundaryConditionsElasticity2DMBB(*modelElasticity2D);
    LinearAnalysis anElasticity2D(modelElasticity2D.get(), SolverType::ELU);

    std::unique_ptr<CompMesh> modelPhaseField = std::make_unique<CompMesh>();
    SetupBoundaryConditionsPhaseFieldMBB(*modelPhaseField);
    TransientAnalysis anPhaseField(modelPhaseField.get(), SolverType::ELU, true);
    // NonLinearAnalysis anPhaseField(modelPhaseField.get(), SolverType::ELU);
    anPhaseField.SetMaxIter(10);

    constexpr auto min_val = 1e-3;
    constexpr auto max_val = 1.0;

    std::function<double(double, double)> u0;
    if (use_vtu_initial_condition)
    {
        u0 = [initial_field = std::move(initial_field), initial_lx, initial_ly](double x, double y) -> double
        {
            return sample_row_wise_grid(initial_field, x, y, initial_lx, initial_ly);
        };
    }
    else
    {
        u0 = [](double x, double y) -> double
        {
        constexpr auto A = 0.2;
        const double V = V_F;
        return V + A * sin(6 * M_PI * x) * sin(3 * M_PI * y);
        };
    }

    for (size_t i = 0; i < anPhaseField.MeshVector()[0]->NNodes(); i++)
    {
        auto &node = anPhaseField.MeshVector()[0]->NodeVec()[i];

        const auto &node_pos = node->getCoordinates();

        auto &fSol = node->Solution();
        auto &fSolPrev = node->PrevSolution();
        fSol.resize(1);
        fSolPrev.resize(1);

        fSol[0] = u0(node_pos[0], node_pos[1]);
        fSolPrev[0] = fSol[0];
    }

    std::vector<std::string> ScalarNamesPhaseField, VectorNamesPhaseField;
    ScalarNamesPhaseField = {"Solution", "TimeDerivative", "WeightFunction"};
    VectorNamesPhaseField = {"Derivative"};

    anPhaseField.PrintVariables("basic_2d_phase_field_dt", ScalarNamesPhaseField, VectorNamesPhaseField);

    if (modelElasticity2D->NElements() != modelPhaseField->NElements())
    {
        // We assume that we are using the same mesh for phase field and elasticity!
        // Yes, but the check need to be different since now we have elements corresponding to the Lagrange multiplier
        // PanicButton();
    }

    for (size_t inode = 0; inode < modelElasticity2D->NNodes(); inode++)
    {
        modelElasticity2D->NodeVec()[inode]->setWeightFunction(1.0);
    }
    for (int64_t i_el = 0; i_el < modelElasticity2D->NElements(); i_el++){
        auto &elemElas2D = modelElasticity2D->ElementVec()[i_el];            
        if (elemElas2D->Dimension() != modelElasticity2D->Dimension()) continue;
        elemElas2D->setIntegPointWeightFunction();
    }

    std::vector<std::string> ScalarNamesElasticity2D, VectorNamesElasticity2D;
    ScalarNamesElasticity2D = {"Compliance", "ComplianceSensibility", "WeightFunction"};
    VectorNamesElasticity2D = {"Displacement", "Stress", "Strain"};

    int i = initial_step;
    
    for (int64_t i_el = 0; i_el < modelElasticity2D->NElements(); i_el++)
    {
        auto elemPhaseField = modelPhaseField->ElementVec()[i_el];

        elemPhaseField->IntegrationData().fJ = 1;
    }

    int step = initial_step;
    // Solving
    while (true)
    {
        // std::cout << "********** OPTIMIZATION STEP " << step << " **********\n";
        // For the elasticity problem, transfer phase field solution to the weightFunction data structure.
        for (size_t inode = 0; inode < modelElasticity2D->NNodes(); inode++)
        {
            if (step == 0){
                modelElasticity2D->NodeVec()[inode]->setWeightFunction(1.0);
            } else {
                modelElasticity2D->NodeVec()[inode]->setWeightFunction(std::min(std::max(std::pow(modelPhaseField->NodeVec()[inode]->Solution()[0], 3), 1e-2), 1.0));
                // modelElasticity2D->NodeVec()[inode]->setWeightFunction(std::max(min_val, pow(modelPhaseField->NodeVec()[inode]->Solution()[0], 3)));
            }
        }
        for (int64_t i_el = 0; i_el < modelElasticity2D->NElements(); i_el++){
            auto &elemElas2D = modelElasticity2D->ElementVec()[i_el];            
            if (elemElas2D->Dimension() != modelElasticity2D->Dimension()) continue;
            elemElas2D->setIntegPointWeightFunction();
        }
        step++;
        // std::cout << "Solving Solid Mechanics\n";
        anElasticity2D.Run();
        VTUGenerator::PrintResults(modelElasticity2D.get(), "cantilever_2d_beam", ScalarNamesElasticity2D, VectorNamesElasticity2D, {}, 0);
        if(i % 200 == 0) {
            VTUGenerator::PrintResults(modelElasticity2D.get(), "cantilever_2d_beam", ScalarNamesElasticity2D, VectorNamesElasticity2D, {}, i);
        }
        std::vector<double> compliances;
        // For each elasticity element, compute the compliance and then transfer the information
        // to the phase field corresponding element. The value of J need to be computed for each integration point.
        // Thus the //J data structure need to be a vector with lenght equal to the number of integration points.
        // Otherwhise, the value of J can also be stored in the weight function data structure.
        for (int64_t i_el = 0; i_el < modelElasticity2D->NElements() && i_el < modelPhaseField->NElements(); i_el++)
        {
            auto elemElas2D = modelElasticity2D->ElementVec()[i_el];
            auto elemPhaseField = modelPhaseField->ElementVec()[i_el];

            if (elemElas2D->Dimension() != modelElasticity2D->Dimension() || elemElas2D->Dimension() != modelElasticity2D->Dimension())
            {
                continue;
            }

            for (int inode = 0; inode < elemElas2D->NElNodes(); inode++)
            {
                auto xparametric = elemElas2D->NodeCoord(inode);
                auto connect = modelElasticity2D->ElementVec()[i_el]->getConnectivity();
                elemElas2D->IntegrationData().fAdimCoord = xparametric;
                if (!elemElas2D->IntegrationData().fNeedsSol || !elemElas2D->IntegrationData().fNeedsDSol)
                {
                    elemElas2D->IntegrationData().fNeedsSol = true;
                    elemElas2D->IntegrationData().fNeedsDSol = true;
                    elemElas2D->IntegrationData().fSol.resize(elemElas2D->GetWeakForm()->NState());
                    elemElas2D->IntegrationData().fDSolDx.resize(elemElas2D->GetWeakForm()->NState(), elemElas2D->Dimension());
                }
                elemElas2D->ComputeJacobian();
                elemElas2D->ComputeSpatialDerivatives();
                elemElas2D->interpolateSolution();
                const auto compliance_var_idx = elemElas2D->GetWeakForm()->VariableIndex("ComplianceSensibility");
                const auto nvar = elemElas2D->GetWeakForm()->NSolutionVariables(compliance_var_idx);
                VecDouble compl_vec(nvar);
                elemElas2D->Solution(compliance_var_idx, compl_vec);
                compliances.push_back(compl_vec[0]);
            }
        }

        const auto maxCompliance = *std::max_element(compliances.begin(), compliances.end());

        int64_t complianceIndex = 0;
        for (int64_t i_el = 0; i_el < modelElasticity2D->NElements() && i_el < modelPhaseField->NElements(); i_el++)
        {
            auto elemElas2D = modelElasticity2D->ElementVec()[i_el];
            auto elemPhaseField = modelPhaseField->ElementVec()[i_el];

            if (elemElas2D->Dimension() != modelElasticity2D->Dimension() || elemPhaseField->Dimension() != modelPhaseField->Dimension())
            {
                continue;
            }

            for (int inode = 0; inode < elemElas2D->NElNodes(); inode++)
            {
                auto connect = modelElasticity2D->ElementVec()[i_el]->getConnectivity();
                modelPhaseField->NodeVec()[connect[inode]]->setWeightFunction(compliances[complianceIndex] / maxCompliance);
                complianceIndex++;
            }
        }

        for (int64_t i_el = 0; i_el < modelPhaseField->NElements() && i_el < modelElasticity2D->NElements(); i_el++){
            auto &elemPhaseField = modelPhaseField->ElementVec()[i_el];            
            if (elemPhaseField->Dimension() != modelPhaseField->Dimension()) continue;
            elemPhaseField->setIntegPointWeightFunction();
        }

        std::cout << "Solving Phase-Field\n";
        for(size_t i = 0; i < 40; i++) {
            anPhaseField.Run(1);
        }
        if(i % 10 == 0) {
            VTUGenerator::PrintResults(modelPhaseField.get(), "phase_field_2d_", ScalarNamesPhaseField, VectorNamesPhaseField, {}, i);
        }

        std::cout << "LAGRANGE MULTIPLIER = " << modelPhaseField->NodeVec().back()->Solution()[0] << '\n';
        i++;
    }

    return 0;
}
